#!/usr/bin/env python3
"""
Export PostgreSQL 'fluxstd' to Hive-partitioned Parquet for DuckDB.

- Reads database connection info and default filters from a TOML file.
- Allows overriding key settings via command-line arguments (argparse).
- Reads Postgres via psycopg2 with server-side cursors for efficient chunking.
- Computes HEALPix columns using astropy-healpix (default HEALPix order = ring).
- Writes Parquet partitioned by hpx32_group in Hive style: hpx32_group=<value>/... [web:188]
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import signal
import sys
import time
import tomllib
from typing import Any

import duckdb
import numpy as np
import pandas as pd
import psycopg2
from astropy import units as u
from astropy.coordinates import Latitude, Longitude
from astropy_healpix import lonlat_to_healpix
from loguru import logger


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export Postgres table to Hive-partitioned Parquet for DuckDB."
    )

    parser.add_argument(
        "-c",
        "--config",
        required=True,
        help="Path to TOML config file (with [targetdb.db] and [targetdb.fluxstd]).",
    )
    parser.add_argument(
        "-d",
        "--outdir",
        "--out",
        dest="out",
        required=True,
        help="Output directory for partitioned Parquet dataset.",
    )

    # Source selection
    parser.add_argument("--schema", default="public", help="Postgres schema name.")
    parser.add_argument("--table", default="fluxstd", help="Postgres table name.")

    # Filtering (override TOML)
    parser.add_argument(
        "--input-catalog-id",
        type=int,
        action="append",
        help="Filter: input_catalog_id (repeatable). Overrides TOML if provided.",
    )
    parser.add_argument(
        "--version",
        action="append",
        help="Filter: version string (repeatable). Overrides TOML if provided.",
    )

    # Chunking
    parser.add_argument(
        "--order-key",
        default="fluxstd_id",
        help="Monotonic numeric key for chunking (default: fluxstd_id).",
    )
    parser.add_argument(
        "--chunk-rows",
        type=int,
        default=2_000_000,
        help="Rows per chunk read from Postgres.",
    )

    # Sky columns
    parser.add_argument("--ra-col", default="ra", help="RA column name (deg, ICRS).")
    parser.add_argument("--dec-col", default="dec", help="Dec column name (deg, ICRS).")
    parser.add_argument(
        "--columns",
        help="Comma-separated list of columns to export (default: predefined subset).",
    )
    parser.add_argument(
        "--all-columns",
        action="store_true",
        help="Export all columns (SELECT *). Overrides --columns.",
    )

    # HEALPix
    parser.add_argument("--nside", type=int, default=32, help="HEALPix NSIDE.")
    parser.add_argument(
        "--healpix-order",
        choices=["ring", "nested"],
        default="ring",
        help="HEALPix ordering scheme (default: ring).",
    )
    parser.add_argument(
        "--group-size",
        type=int,
        default=64,
        help="Grouping factor: hpx_group = hpx // group_size.",
    )

    # Performance tuning
    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Number of DuckDB threads for parallel processing (default: 16).",
    )

    # Cleanup
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Clean up temporary directory after compaction.",
    )

    # Compaction-only mode
    parser.add_argument(
        "--compact-only",
        action="store_true",
        help="Skip data fetch and only perform compaction from existing tmp directory.",
    )

    return parser.parse_args()


def load_toml(path: str) -> dict[str, dict[str, Any]]:
    with open(path, "rb") as f:
        return tomllib.load(f)  # [web:455]


def validate_sql_identifier(name: str, param_name: str) -> str:
    """
    Validate SQL identifier to prevent injection attacks.

    Args:
        name: The identifier to validate
        param_name: Name of the parameter (for error messages)

    Returns
    -------
        The validated identifier

    Raises
    ------
        ValueError: If the identifier contains invalid characters
    """
    # PostgreSQL/DuckDB identifier rules: start with letter or underscore,
    # followed by letters, digits, underscores, or dollar signs
    if not re.match(r"^[a-zA-Z_][a-zA-Z0-9_$]*$", name):
        raise ValueError(
            f"Invalid {param_name}: '{name}'. "
            f"SQL identifiers must start with a letter or underscore, "
            f"and contain only letters, numbers, underscores, or dollar signs."
        )

    # Additional safety: reject SQL keywords that could be dangerous
    sql_keywords = {
        "SELECT",
        "DROP",
        "DELETE",
        "INSERT",
        "UPDATE",
        "CREATE",
        "ALTER",
        "GRANT",
        "REVOKE",
        "UNION",
        "EXEC",
        "EXECUTE",
    }
    if name.upper() in sql_keywords:
        raise ValueError(f"Invalid {param_name}: '{name}' is a reserved SQL keyword")

    return name


def validate_output_path(path: str) -> str:
    """
    Validate output directory path to prevent path traversal attacks.

    Args:
        path: The output path to validate

    Returns
    -------
        The absolute path

    Raises
    ------
        ValueError: If the path contains suspicious patterns
    """
    # Check for path traversal attempts
    if ".." in path:
        raise ValueError(f"Path traversal detected in output path: '{path}'")

    # Convert to absolute path
    abs_path = os.path.abspath(path)

    return abs_path


def build_pg_connstr(db_cfg: dict[str, str | int]) -> str:
    # libpq connection string for psycopg2
    return (
        f"host={db_cfg['host']} "
        f"port={db_cfg.get('port', 5432)} "
        f"dbname={db_cfg['dbname']} "
        f"user={db_cfg['user']} "
        f"password={db_cfg['password']}"
    )


def build_pg_dict(db_cfg: dict[str, str | int]) -> dict[str, Any]:
    # Dictionary format for psycopg2.connect()
    return {
        "host": db_cfg["host"],
        "port": db_cfg.get("port", 5432),
        "dbname": db_cfg["dbname"],
        "user": db_cfg["user"],
        "password": db_cfg["password"],
    }


def list_or_none[T](xs: list[T] | None) -> list[T] | None:
    return xs if (xs is not None and len(xs) > 0) else None


def get_default_columns() -> list[str]:
    """
    Get the default list of columns to export.

    Returns
    -------
        List of default column names
    """
    return [
        "fluxstd_id",
        "obj_id",
        "ra",
        "dec",
        "epoch",
        "parallax",
        "parallax_error",
        "pmra",
        "pmra_error",
        "pmdec",
        "pmdec_error",
        "input_catalog_id",
        "version",
        "psf_flux_g",
        "psf_flux_r",
        "psf_flux_i",
        "psf_flux_z",
        "psf_flux_y",
        "psf_flux_error_g",
        "psf_flux_error_r",
        "psf_flux_error_i",
        "psf_flux_error_z",
        "psf_flux_error_y",
        "filter_g",
        "filter_r",
        "filter_i",
        "filter_z",
        "filter_y",
        "prob_f_star",
        "teff_brutus",
        "teff_brutus_low",
        "teff_brutus_high",
        "logg_brutus",
        "logg_brutus_low",
        "logg_brutus_high",
        "teff_gspphot",
        "teff_gspphot_lower",
        "teff_gspphot_upper",
        "is_fstar_gaia",
    ]


def build_column_list(
    columns_arg: str | None,
    all_columns: bool,
    order_key: str,
    ra_col: str,
    dec_col: str,
) -> str:
    """
    Build SQL column list, ensuring required columns are included.

    Args:
        columns_arg: Comma-separated column list from user (or None for default)
        all_columns: If True, return "*" to export all columns
        order_key: Order key column (must be included)
        ra_col: RA column name (must be included for HEALPix)
        dec_col: Dec column name (must be included for HEALPix)

    Returns
    -------
        SQL column list string (e.g., "col1, col2, col3" or "*")
    """
    # If --all-columns flag is set, return "*"
    if all_columns:
        logger.info("Using --all-columns: exporting all columns")
        return "*"

    # Determine which columns to use
    if columns_arg is None:
        # Use default column list
        cols = get_default_columns()
        logger.info("Using default column list")
    else:
        # Parse user-provided columns
        cols = [c.strip() for c in columns_arg.split(",") if c.strip()]

    # Ensure required columns are included
    required = {order_key, ra_col, dec_col}
    cols_set = set(cols)

    missing = required - cols_set
    if missing:
        logger.warning(f"Adding required columns: {', '.join(missing)}")
        cols.extend(missing)

    return ", ".join(cols)


def build_filters(
    args: argparse.Namespace, toml_cfg: dict[str, Any]
) -> tuple[list[int] | None, list[str] | None]:
    flux_cfg = toml_cfg["targetdb"]["fluxstd"]

    input_ids = list_or_none(args.input_catalog_id)
    versions = list_or_none(args.version)

    if input_ids is None:
        input_ids = flux_cfg.get("input_catalog_id", None)
    if versions is None:
        versions = flux_cfg.get("version", None)

    return input_ids, versions


def make_where_clause(input_ids: list[int] | None, versions: list[str] | None) -> str:
    clauses = []
    if input_ids:
        ids_sql = ", ".join(str(int(x)) for x in input_ids)
        clauses.append(f"input_catalog_id IN ({ids_sql})")
    elif versions:
        ver_sql = ", ".join("'" + str(v).replace("'", "''") + "'" for v in versions)
        clauses.append(f"version IN ({ver_sql})")
    return ("WHERE " + " AND ".join(clauses)) if clauses else ""


def add_healpix_columns(
    df: pd.DataFrame,
    ra_col: str,
    dec_col: str,
    nside: int,
    healpix_order: str,
    group_size: int,
    hpx_col: str,
    hpx_group_col: str,
) -> pd.DataFrame:
    ra = df[ra_col].to_numpy(dtype=np.float64)
    dec = df[dec_col].to_numpy(dtype=np.float64)

    lon = Longitude(ra, u.deg)
    lat = Latitude(dec, u.deg)

    hpx = lonlat_to_healpix(lon, lat, nside=nside, order=healpix_order).astype(np.int64)

    out = df.copy()
    out[hpx_col] = hpx
    out[hpx_group_col] = (out[hpx_col] // group_size).astype(np.int64)
    return out


def compact_partitions(
    tmp_path: str,
    release_path: str,
    hpx_group_col: str,
    threads: int,
) -> int:
    """
    Compact partitioned Parquet files from tmp to release directory.

    Args:
        tmp_path: Path to temporary directory with partitioned data
        release_path: Path to release directory for compacted data
        hpx_group_col: Name of the partition column
        threads: Number of DuckDB threads to use

    Returns
    -------
        Total number of rows in the compacted dataset
    """
    logger.info("Starting Parquet compaction...")
    logger.info(f"Source: {tmp_path}")
    logger.info(f"Destination: {release_path}")

    # Create release directory
    os.makedirs(release_path, exist_ok=True)

    # Connect to DuckDB
    con = duckdb.connect(database=":memory:")
    con.execute(f"SET threads={threads};")

    try:
        # Read all partitions and write compacted output
        # COPY with PARTITION_BY will create one file per partition instead of many
        compact_start = time.time()
        logger.info("Reading and compacting partitions...")

        con.execute(
            f"""
            COPY (SELECT * FROM read_parquet('{tmp_path}/**/*.parquet'))
            TO '{release_path}'
            (FORMAT parquet,
             PARTITION_BY ({hpx_group_col}),
             FILENAME_PATTERN 'data')
        """
        )

        compact_time = time.time() - compact_start
        logger.info(f"Compaction completed in {compact_time:.1f}s")

        # Verify row count
        count_result = con.execute(
            f"SELECT COUNT(*) FROM read_parquet('{release_path}/**/*.parquet')"
        ).fetchone()
        total_rows = count_result[0] if count_result else 0

        logger.success(f"Compacted {total_rows:,} rows to {release_path}")

        return total_rows

    finally:
        con.close()


def signal_handler(sig, frame):
    logger.warning("Interrupt signal received (Ctrl-C)")
    logger.info("Waiting for current database operation to complete...")
    logger.info("This may take a few moments if a query is in progress")
    # Let the exception propagate naturally
    raise KeyboardInterrupt


def main():
    # Set up signal handler
    signal.signal(signal.SIGINT, signal_handler)

    args = parse_args()

    # Validate output path to prevent path traversal
    base_out_path = validate_output_path(args.out)

    # Set up tmp and release directories
    tmp_path = os.path.join(base_out_path, "tmp")
    release_path = os.path.join(base_out_path, "release")

    # Generate HEALPix column names based on nside
    hpx_col = f"hpx{args.nside}"
    hpx_group_col = f"hpx{args.nside}_group"

    # Compact-only mode: skip data fetch and only perform compaction
    if args.compact_only:
        logger.info("=" * 80)
        logger.info("Running in COMPACT-ONLY mode")
        logger.info("=" * 80)
        logger.info(f"Base output directory: {base_out_path}")
        logger.info(f"Temporary directory: {tmp_path}")
        logger.info(f"Release directory: {release_path}")
        logger.info(f"HEALPix group column: {hpx_group_col}")

        # Validate that tmp directory exists and has data
        if not os.path.exists(tmp_path):
            logger.error(f"Temporary directory does not exist: {tmp_path}")
            logger.error("Cannot run --compact-only without existing tmp data.")
            logger.error("Please run without --compact-only first to fetch data.")
            sys.exit(1)

        # Check if tmp has any parquet files
        parquet_files = list(os.walk(tmp_path))
        has_parquet = any(
            f.endswith(".parquet") for _, _, files in parquet_files for f in files
        )
        if not has_parquet:
            logger.error(f"No Parquet files found in: {tmp_path}")
            logger.error("Cannot run --compact-only without existing data.")
            sys.exit(1)

        logger.info(f"Found existing data in {tmp_path}")

        # Perform compaction
        logger.info("=" * 80)
        logger.info("Starting compaction phase")
        logger.info("=" * 80)

        compacted_rows = compact_partitions(
            tmp_path=tmp_path,
            release_path=release_path,
            hpx_group_col=hpx_group_col,
            threads=args.threads,
        )

        logger.success(
            f"Compaction complete! {compacted_rows:,} rows written to {release_path}"
        )

        # Clean up tmp directory if --clean flag is set
        if args.clean:
            logger.info("Cleaning up temporary directory...")
            try:
                shutil.rmtree(tmp_path)
                logger.success(f"Removed temporary directory: {tmp_path}")
            except Exception as e:
                logger.error(f"Failed to remove temporary directory: {e}")
        else:
            logger.info(
                f"Temporary directory preserved at: {tmp_path}\n"
                f"Use --clean flag to automatically remove it after compaction."
            )

        return

    # Normal mode: fetch data and then compact
    logger.info(f"Loading configuration from {args.config}")
    cfg = load_toml(args.config)

    db_cfg = cfg["targetdb"]["db"]
    pg_conn_dict = build_pg_dict(db_cfg)

    input_ids, versions = build_filters(args, cfg)
    base_where = make_where_clause(input_ids, versions)

    # Validate SQL identifiers to prevent injection attacks
    schema = validate_sql_identifier(args.schema, "schema name")
    table = validate_sql_identifier(args.table, "table name")
    order_key = validate_sql_identifier(args.order_key, "order key")
    ra_col = validate_sql_identifier(args.ra_col, "RA column")
    dec_col = validate_sql_identifier(args.dec_col, "Dec column")

    # Build column list, ensuring required columns are included
    column_list = build_column_list(
        args.columns, args.all_columns, order_key, ra_col, dec_col
    )

    logger.info(f"Base output directory: {base_out_path}")
    logger.info(f"Temporary directory: {tmp_path}")
    logger.info(f"Release directory: {release_path}")
    logger.info(f"Columns to export: {column_list}")
    logger.info(f"HEALPix columns: {hpx_col}, {hpx_group_col} (nside={args.nside})")

    # Guard: Check if base output directory already exists with content
    # This prevents accidental overwrite of existing data in "full export to new dir" workflow
    if os.path.exists(base_out_path):
        contents = os.listdir(base_out_path)
        if contents:
            logger.error(
                f"Output directory already exists and is not empty: {base_out_path}\n"
                f"Found {len(contents)} items: {contents[:5]}{'...' if len(contents) > 5 else ''}\n"
                f"For safety, this script will not overwrite existing data.\n"
                f"Please use a new/empty directory for each export (e.g., releases/<tag>/)."
            )
            sys.exit(1)
        logger.info("Output directory exists but is empty - proceeding")
    else:
        os.makedirs(base_out_path)
        logger.info(f"Created output directory: {base_out_path}")

    # Create tmp directory for initial export
    os.makedirs(tmp_path, exist_ok=True)
    logger.info(f"Created temporary directory: {tmp_path}")

    # Connect to PostgreSQL using psycopg2
    pg_conn = None
    cursor = None
    duckdb_con = None

    try:
        logger.info(f"Connecting to PostgreSQL: {schema}.{table}")
        pg_conn = psycopg2.connect(**pg_conn_dict)

        # Use a server-side cursor for efficient chunking
        # Named cursors in psycopg2 are server-side cursors
        cursor_name = "fetchiter"
        cursor = pg_conn.cursor(name=cursor_name)
        cursor.arraysize = args.chunk_rows  # Set fetch size

        # Build the SQL query
        if base_where:
            where_sql = base_where
        else:
            where_sql = ""

        sql = f"""
            SELECT {column_list}
            FROM {schema}.{table}
            {where_sql}
            ORDER BY {order_key}
        """

        logger.info("Executing query with server-side cursor...")
        logger.debug(f"SQL query: {sql.strip()}")
        cursor.execute(sql)

        # Initialize DuckDB for Parquet writing
        duckdb_con = duckdb.connect(database=":memory:")
        duckdb_con.execute(f"SET threads={args.threads};")
        logger.info(f"DuckDB threads set to {args.threads}")

        total = 0
        start_time = time.time()
        chunk_num = 0

        while True:
            chunk_start = time.time()
            chunk_num += 1

            logger.info(f"Fetching chunk {chunk_num} ({args.chunk_rows:,} rows)...")
            fetch_start = time.time()

            # Fetch chunk using server-side cursor
            rows = cursor.fetchmany(args.chunk_rows)
            fetch_time = time.time() - fetch_start

            if not rows:
                logger.info("No more rows to fetch")
                break

            # Get column names from cursor description
            column_names = [desc[0] for desc in cursor.description]

            # Convert to pandas DataFrame
            df = pd.DataFrame(rows, columns=column_names)
            logger.info(f"Fetched {len(df):,} rows in {fetch_time:.1f}s")

            logger.info("Computing HEALPix columns...")
            healpix_start = time.time()
            df = add_healpix_columns(
                df,
                ra_col=ra_col,
                dec_col=dec_col,
                nside=args.nside,
                healpix_order=args.healpix_order,
                group_size=args.group_size,
                hpx_col=hpx_col,
                hpx_group_col=hpx_group_col,
            )
            healpix_time = time.time() - healpix_start
            logger.info(f"HEALPix computed in {healpix_time:.1f}s")

            duckdb_con.register("chunk_df", df)

            # DuckDB COPY with PARTITION_BY writes Hive-partitioned datasets
            # Use APPEND to accumulate chunks into the same partitioned output directory
            # FILENAME_PATTERN with {uuid} prevents filename collisions across chunks
            logger.info("Writing to Parquet...")
            write_start = time.time()
            duckdb_con.execute(
                f"""
                COPY (SELECT * FROM chunk_df)
                TO '{tmp_path}'
                (FORMAT parquet,
                 PARTITION_BY ({hpx_group_col}),
                 APPEND,
                 FILENAME_PATTERN 'part_{{uuid}}')
            """
            )
            write_time = time.time() - write_start
            logger.info(f"Written to Parquet in {write_time:.1f}s")

            last_key = int(df[order_key].max())
            total += len(df)

            # Calculate progress metrics
            chunk_time = time.time() - chunk_start
            elapsed = time.time() - start_time
            rate = total / elapsed if elapsed > 0 else 0

            logger.success(
                f"Chunk {chunk_num} complete in {chunk_time:.1f}s | "
                f"Total: {total:,} rows | "
                f"Overall rate: {rate:,.0f} rows/sec | "
                f"Elapsed: {elapsed:.1f}s | "
                f"Last {order_key}={last_key}"
            )

        elapsed_total = time.time() - start_time
        logger.success(
            f"Export complete! Total: {total:,} rows in {elapsed_total:.1f}s "
            f"({total/elapsed_total:,.0f} rows/sec)"
        )

        # Verify the exported data in tmp
        logger.info("Verifying exported Parquet files in tmp...")
        verify_start = time.time()
        try:
            verify_con = duckdb.connect(database=":memory:")
            count_result = verify_con.execute(
                f"SELECT COUNT(*) FROM read_parquet('{tmp_path}/**/*.parquet')"
            ).fetchone()
            parquet_count = count_result[0] if count_result else 0
            verify_time = time.time() - verify_start

            logger.info(f"Verification completed in {verify_time:.1f}s")
            logger.info(f"Rows exported to PostgreSQL: {total:,}")
            logger.info(f"Rows found in Parquet files: {parquet_count:,}")

            if parquet_count == total:
                logger.success("✓ Row count matches! Data integrity verified.")
            else:
                diff = abs(parquet_count - total)
                logger.error(
                    f"✗ Row count mismatch! Difference: {diff:,} rows "
                    f"({parquet_count:,} in Parquet vs {total:,} exported)"
                )
                logger.warning(
                    "This may indicate data loss during export. "
                    "Please investigate before using the data."
                )
            verify_con.close()
        except Exception as e:
            logger.error(f"Failed to verify Parquet files: {e}")
            logger.warning("Skipping verification, but export was completed")

        # Compact partitions from tmp to release
        logger.info("=" * 80)
        logger.info("Starting compaction phase")
        logger.info("=" * 80)

        compacted_rows = compact_partitions(
            tmp_path=tmp_path,
            release_path=release_path,
            hpx_group_col=hpx_group_col,
            threads=args.threads,
        )

        logger.success(
            f"Compaction complete! {compacted_rows:,} rows written to {release_path}"
        )

        # Clean up tmp directory if --clean flag is set
        if args.clean:
            logger.info("Cleaning up temporary directory...")
            try:
                shutil.rmtree(tmp_path)
                logger.success(f"Removed temporary directory: {tmp_path}")
            except Exception as e:
                logger.error(f"Failed to remove temporary directory: {e}")
        else:
            logger.info(
                f"Temporary directory preserved at: {tmp_path}\n"
                f"Use --clean flag to automatically remove it after compaction."
            )

    except KeyboardInterrupt:
        elapsed = time.time() - start_time
        logger.warning(f"Export interrupted by user after {elapsed:.1f}s")
        logger.info(f"Partial export: {total:,} rows exported before interruption")
        logger.warning(
            "Note: Automatic resume is not implemented. Next run will start from the beginning."
        )
        sys.exit(1)

    finally:
        # Clean up connections
        if cursor:
            cursor.close()
        if pg_conn:
            pg_conn.close()
        if duckdb_con:
            duckdb_con.close()
        logger.info("Database connections closed")


if __name__ == "__main__":
    main()
