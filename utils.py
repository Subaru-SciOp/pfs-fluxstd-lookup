#!/usr/bin/env python3
"""
Flux standard catalog matching utilities using HEALPix spatial indexing.

This module provides functions to match astronomical catalogs with a flux standard
catalog stored in DuckDB-readable Parquet files. The matching process uses HEALPix
indexing for efficient spatial queries and applies proper motion and parallax
corrections to coordinates.

Main Functions
--------------
add_healpix_columns : Add HEALPix indices to a DataFrame
fetch_fluxstd : Query flux standard catalog by HEALPix indices
generate_fluxstd_coords : Generate SkyCoord with proper motion corrections
match_catalog : Perform sky coordinate matching between catalogs

Examples
--------
>>> import pandas as pd
>>> df = pd.read_csv("input_catalog.csv")
>>> df = add_healpix_columns(df, nside=32)
>>> hpx_indices = df["hpx32"].unique().tolist()
>>> df_fluxstd = fetch_fluxstd(hpx_indices, "./data/fluxstd/")
>>> df_matched = match_catalog(df, df_fluxstd)
"""
import math
import os
import time
import warnings
from pathlib import Path
from typing import Any

import duckdb
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import (
    Distance,
    Latitude,
    Longitude,
    SkyCoord,
    match_coordinates_sky,
)
from astropy.time import Time
from astropy_healpix import lonlat_to_healpix, neighbours
from dotenv import load_dotenv
from erfa import ErfaWarning
from loguru import logger

load_dotenv(override=False)


def _env_int(name: str, default: int) -> int:
    value = os.getenv(name)
    if value is None or value.strip() == "":
        return default
    try:
        return int(value)
    except ValueError:
        logger.warning(f"Invalid env var {name}='{value}', using default={default}")
        return default


def _env_float(name: str, default: float) -> float:
    value = os.getenv(name)
    if value is None or value.strip() == "":
        return default
    try:
        return float(value)
    except ValueError:
        logger.warning(f"Invalid env var {name}='{value}', using default={default}")
        return default


def _env_str(name: str, default: str, *, allowed: set[str] | None = None) -> str:
    value = os.getenv(name)
    if value is None or value.strip() == "":
        return default
    value = value.strip()
    if allowed is not None and value not in allowed:
        logger.warning(
            f"Invalid env var {name}='{value}', allowed={sorted(allowed)}; using default={default}"
        )
        return default
    return value


def _env_bool(name: str, default: bool) -> bool:
    """
    Get boolean value from environment variable.

    Parameters
    ----------
    name : str
        Environment variable name
    default : bool
        Default value if not set or invalid

    Returns
    -------
    bool
        Boolean value from environment or default
    """
    value = os.getenv(name)
    if value is None or value.strip() == "":
        return default
    value = value.strip().lower()
    if value in ("true", "1", "yes", "on"):
        return True
    elif value in ("false", "0", "no", "off"):
        return False
    else:
        logger.warning(f"Invalid env var {name}='{value}', using default={default}")
        return default


# Configuration constants (env-overridable)
DEFAULT_DUCKDB_THREADS = _env_int("DUCKDB_THREADS", 16)
"""Default number of threads for DuckDB operations (env: DUCKDB_THREADS)"""

DEFAULT_PARALLAX_FALLBACK = _env_float("FLUXSTD_PARALLAX_FALLBACK_MAS", 1e-7)
"""Default parallax value (mas) for objects with poor parallax measurements.
Env: FLUXSTD_PARALLAX_FALLBACK_MAS.

This very small but non-zero value prevents division by zero while maintaining
the correct sign for distance calculations.
"""

DEFAULT_PM_SNR_THRESHOLD = _env_float("FLUXSTD_PM_SNR_THRESHOLD", 3.0)
"""SNR threshold for accepting proper motion and parallax (env: FLUXSTD_PM_SNR_THRESHOLD)"""

DEFAULT_MATCH_SEPARATION = _env_float("FLUXSTD_MATCH_SEPARATION_ARCSEC", 1.0) * u.arcsec
"""Default maximum separation for catalog matching (env: FLUXSTD_MATCH_SEPARATION_ARCSEC)"""

DEFAULT_NSIDE = _env_int("FLUXSTD_NSIDE", 32)
"""Default HEALPix NSIDE parameter (env: FLUXSTD_NSIDE)"""

DEFAULT_HEALPIX_ORDER = _env_str(
    "FLUXSTD_HEALPIX_ORDER", "ring", allowed={"ring", "nested"}
)
"""Default HEALPix ordering scheme (env: FLUXSTD_HEALPIX_ORDER)"""

DEFAULT_BATCH_SIZE = _env_int("FLUXSTD_BATCH_SIZE", 20)
"""Number of HEALPix pixels to process per batch (env: FLUXSTD_BATCH_SIZE)"""

DEFAULT_INCLUDE_NEIGHBORS = _env_bool("FLUXSTD_INCLUDE_NEIGHBORS", True)
"""Include neighboring HEALPix pixels in search (env: FLUXSTD_INCLUDE_NEIGHBORS)"""

# Flux standard quality selection criteria (env-overridable)
DEFAULT_MIN_PROB_F_STAR = _env_float("FLUXSTD_MIN_PROB_F_STAR", 0.5)
"""Minimum stellar probability (env: FLUXSTD_MIN_PROB_F_STAR)"""

DEFAULT_MAX_PROB_F_STAR = _env_float("FLUXSTD_MAX_PROB_F_STAR", 1.0)
"""Maximum stellar probability (env: FLUXSTD_MAX_PROB_F_STAR)"""

DEFAULT_MIN_PSF_MAG_G = _env_float("FLUXSTD_MIN_PSF_MAG_G", 17.0)
"""Minimum PS1 g-band PSF magnitude (env: FLUXSTD_MIN_PSF_MAG_G)"""

DEFAULT_MAX_PSF_MAG_G = _env_float("FLUXSTD_MAX_PSF_MAG_G", 19.0)
"""Maximum PS1 g-band PSF magnitude (env: FLUXSTD_MAX_PSF_MAG_G)"""

DEFAULT_MIN_PSF_MAG_R = _env_float("FLUXSTD_MIN_PSF_MAG_R", 17.0)
"""Minimum Gaia G-filter PSF magnitude (env: FLUXSTD_MIN_PSF_MAG_R)"""

DEFAULT_MAX_PSF_MAG_R = _env_float("FLUXSTD_MAX_PSF_MAG_R", 19.0)
"""Maximum Gaia G-filter PSF magnitude (env: FLUXSTD_MAX_PSF_MAG_R)"""

DEFAULT_MIN_TEFF_BRUTUS = _env_float("FLUXSTD_MIN_TEFF_BRUTUS", 5000.0)
"""Minimum Teff (K) from Brutus (env: FLUXSTD_MIN_TEFF_BRUTUS)"""

DEFAULT_MAX_TEFF_BRUTUS = _env_float("FLUXSTD_MAX_TEFF_BRUTUS", 8000.0)
"""Maximum Teff (K) from Brutus (env: FLUXSTD_MAX_TEFF_BRUTUS)"""

DEFAULT_MIN_TEFF_GSPPHOT = _env_float("FLUXSTD_MIN_TEFF_GSPPHOT", 5000.0)
"""Minimum Teff (K) from Gaia GSP-Phot (env: FLUXSTD_MIN_TEFF_GSPPHOT)"""

DEFAULT_MAX_TEFF_GSPPHOT = _env_float("FLUXSTD_MAX_TEFF_GSPPHOT", 8000.0)
"""Maximum Teff (K) from Gaia GSP-Phot (env: FLUXSTD_MAX_TEFF_GSPPHOT)"""


def _validate_path(path: str, param_name: str) -> str:
    """
    Validate file system path to prevent path traversal attacks.

    Parameters
    ----------
    path : str
        The path to validate
    param_name : str
        Name of the parameter (for error messages)

    Returns
    -------
    str
        The absolute, normalized path

    Raises
    ------
    ValueError
        If the path contains suspicious patterns
    FileNotFoundError
        If the path does not exist
    """
    # Check for null bytes (common injection vector)
    if "\0" in path:
        raise ValueError(f"Invalid {param_name}: path contains null bytes")

    # Convert to absolute path to normalize
    abs_path = str(Path(path).resolve())

    # Check path exists
    if not Path(abs_path).exists():
        raise FileNotFoundError(f"Path does not exist: {abs_path}")

    return abs_path


def add_healpix_columns(
    df: pd.DataFrame,
    ra_col: str = "ra",
    dec_col: str = "dec",
    nside: int = DEFAULT_NSIDE,
    healpix_order: str = DEFAULT_HEALPIX_ORDER,
) -> pd.DataFrame:
    """
    Add HEALPix index columns to a DataFrame based on RA/Dec coordinates.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with celestial coordinates
    ra_col : str, default="ra"
        Name of RA column in degrees
    dec_col : str, default="dec"
        Name of Dec column in degrees
    nside : int, default=32
        HEALPix NSIDE parameter
    healpix_order : str, default="ring"
        HEALPix ordering scheme, 'ring' or 'nested'

    Returns
    -------
    pd.DataFrame
        Copy of input DataFrame with additional hpx{nside} column

    Raises
    ------
    ValueError
        If df is empty, columns are missing, or parameters are invalid
    """
    # Validate inputs
    if df.empty:
        raise ValueError("Input DataFrame is empty")

    if ra_col not in df.columns:
        raise ValueError(
            f"RA column '{ra_col}' not found in DataFrame. "
            f"Available columns: {list(df.columns)}"
        )

    if dec_col not in df.columns:
        raise ValueError(
            f"Dec column '{dec_col}' not found in DataFrame. "
            f"Available columns: {list(df.columns)}"
        )

    if healpix_order not in ["ring", "nested"]:
        raise ValueError(
            f"healpix_order must be 'ring' or 'nested', got '{healpix_order}'"
        )

    if nside <= 0 or (nside & (nside - 1)) != 0:
        raise ValueError(f"nside must be a positive power of 2, got {nside}")

    # Check for null values
    null_ra = df[ra_col].isna().sum()
    null_dec = df[dec_col].isna().sum()
    if null_ra > 0 or null_dec > 0:
        raise ValueError(
            f"Coordinate columns contain null values: {null_ra} in {ra_col}, {null_dec} in {dec_col}"
        )

    ra = df[ra_col].to_numpy(dtype=np.float64)
    dec = df[dec_col].to_numpy(dtype=np.float64)

    lon = Longitude(ra, u.deg)
    lat = Latitude(dec, u.deg)

    logger.info(
        f"Generating HEALPix indices for {len(ra)} input objects with nside={nside} and order={healpix_order}"
    )

    hpx = lonlat_to_healpix(lon, lat, nside=nside, order=healpix_order).astype(np.int64)

    logger.info(f"Number of unique HEALPix indices: {len(np.unique(hpx))}")

    dfout = df.copy()
    dfout[f"hpx{nside}"] = hpx

    return dfout


def fetch_fluxstd(
    hpx_indices: list[int], duckdb_data_root: str = ".", nside: int = DEFAULT_NSIDE
) -> pd.DataFrame:
    """
    Fetch flux standard catalog data from DuckDB-readable Parquet files.

    Queries Parquet files filtered by HEALPix indices.
    This allows efficient loading of only the relevant subset of the catalog.

    Parameters
    ----------
    hpx_indices : list[int]
        List of unique HEALPix indices to query
    duckdb_data_root : str, default="."
        Root directory containing Parquet files
    nside : int, default=32
        HEALPix NSIDE parameter

    Returns
    -------
    pd.DataFrame
        DataFrame containing flux standard catalog entries matching the HEALPix indices

    Raises
    ------
    ValueError
        If hpx_indices is empty or contains invalid values
    FileNotFoundError
        If duckdb_data_root does not exist
    """
    # Validate inputs
    if not hpx_indices:
        raise ValueError("hpx_indices list is empty")

    # Validate path
    validated_path = _validate_path(duckdb_data_root, "duckdb_data_root")

    # Convert HEALPix indices to integers to prevent SQL injection
    try:
        safe_hpx_indices = [int(x) for x in hpx_indices]
    except (ValueError, TypeError) as e:
        raise ValueError(f"HEALPix indices must be integers: {e}") from e

    # Validate nside
    if not isinstance(nside, int) or nside <= 0:
        raise ValueError(f"nside must be a positive integer, got {nside}")

    with duckdb.connect(database=":memory:") as duckdb_conn:
        logger.info("Loading fluxstd catalog into DuckDB")
        duckdb_conn.execute(f"SET threads={DEFAULT_DUCKDB_THREADS};")

        logger.info(
            f"Querying fluxstd catalog by {len(safe_hpx_indices)} unique HEALPix indices"
        )

        # Build safe SQL query with validated integers
        hpx_list = ", ".join(str(int(x)) for x in safe_hpx_indices)

        try:
            query_result = duckdb_conn.execute(
                f"""
                SELECT * FROM read_parquet('{validated_path}/**/*.parquet', hive_partitioning = true)
                WHERE hpx{int(nside)} IN ({hpx_list})
                """
            )
            df_db = query_result.df()
        except Exception as e:
            logger.error(f"DuckDB query failed: {e}")
            raise

        if df_db.empty:
            logger.warning(
                "No flux standard catalog entries found matching the HEALPix indices"
            )
        else:
            logger.info(f"Retrieved {len(df_db)} flux standard catalog entries")

    return df_db


def generate_fluxstd_coords(
    df_db: pd.DataFrame,
) -> tuple[SkyCoord, pd.Series, pd.Series]:
    """
    Generate SkyCoord objects for flux standard catalog with proper motion and parallax.

    Applies proper motion and parallax corrections, propagating coordinates to epoch 2000.0.
    Bad proper motion and parallax measurements (SNR < 3.0) are set to zero/default values.
    Prioritizes PS1+Gaia flux standards, using Gaia-only as fallback.

    Parameters
    ----------
    df_db : pd.DataFrame
        DataFrame with flux standard catalog data

    Returns
    -------
    tuple[SkyCoord, pd.Series, pd.Series]
        Tuple containing:
        - SkyCoord object with coordinates propagated to epoch 2000.0
        - Boolean pandas Series indicating good flux standard objects
        - String pandas Series indicating catalog source ("PS1+Gaia" or "Gaia-only")

    Raises
    ------
    ValueError
        If required columns are missing or df_db is empty
    """
    # Validate inputs
    if df_db.empty:
        raise ValueError("Input DataFrame df_db is empty")

    # Use module-level constants for quality selection criteria
    min_prob_f_star = DEFAULT_MIN_PROB_F_STAR
    max_prob_f_star = DEFAULT_MAX_PROB_F_STAR

    min_psf_mag_g = DEFAULT_MIN_PSF_MAG_G
    max_psf_mag_g = DEFAULT_MAX_PSF_MAG_G

    min_psf_mag_r = DEFAULT_MIN_PSF_MAG_R
    max_psf_mag_r = DEFAULT_MAX_PSF_MAG_R

    min_teff_brutus = DEFAULT_MIN_TEFF_BRUTUS
    max_teff_brutus = DEFAULT_MAX_TEFF_BRUTUS

    min_teff_gspphot = DEFAULT_MIN_TEFF_GSPPHOT
    max_teff_gspphot = DEFAULT_MAX_TEFF_GSPPHOT

    min_psf_flux_g = (max_psf_mag_g * u.ABmag).to(u.nJy).value
    max_psf_flux_g = (min_psf_mag_g * u.ABmag).to(u.nJy).value

    min_psf_flux_r = (max_psf_mag_r * u.ABmag).to(u.nJy).value
    max_psf_flux_r = (min_psf_mag_r * u.ABmag).to(u.nJy).value

    good_fstar_brutus = (df_db["prob_f_star"] >= min_prob_f_star) & (
        df_db["prob_f_star"] <= max_prob_f_star
    )
    good_flux_ps1 = (df_db["psf_flux_g"] >= min_psf_flux_g) & (
        df_db["psf_flux_g"] <= max_psf_flux_g
    )
    good_teff_brutus = (df_db["teff_brutus"] >= min_teff_brutus) & (
        df_db["teff_brutus"] <= max_teff_brutus
    )
    good_fluxstd_ps1 = good_fstar_brutus & good_flux_ps1 & good_teff_brutus

    # Detect stars with PS1 data available
    has_ps1_data = df_db["prob_f_star"].notna()

    good_fstar_gaia = df_db["is_fstar_gaia"]
    good_teff_gaia = (df_db["teff_gspphot"] >= min_teff_gspphot) & (
        df_db["teff_gspphot"] <= max_teff_gspphot
    )
    good_flux_gaia = (df_db["psf_flux_r"] >= min_psf_flux_r) & (
        df_db["psf_flux_r"] <= max_psf_flux_r
    )
    good_fluxstd_gaia = good_fstar_gaia & good_teff_gaia & good_flux_gaia

    if not good_fluxstd_ps1.any():
        logger.warning(
            "No good PS1-Gaia flux standard stars found after applying quality cuts."
        )

    if not good_fluxstd_gaia.any():
        logger.warning(
            "No good Gaia-only flux standard stars found after applying quality cuts."
        )

    # Exclusive priority selection: stars with PS1 data use ONLY PS1 criteria
    # Stars without PS1 data use ONLY Gaia criteria
    # Stars with PS1 data that fail PS1 criteria are excluded entirely
    good_fluxstd_ps1_exclusive = has_ps1_data & good_fluxstd_ps1
    good_fluxstd_gaia_exclusive = ~has_ps1_data & good_fluxstd_gaia
    good_fluxstd = good_fluxstd_ps1_exclusive | good_fluxstd_gaia_exclusive

    if not good_fluxstd.any():
        raise ValueError(
            "No good flux standard stars found after applying quality cuts"
        )

    # Create catalog source labels
    catalog_source = pd.Series([""] * len(df_db), index=df_db.index, dtype=str)
    catalog_source[good_fluxstd_ps1_exclusive] = "PS1+Gaia"
    catalog_source[good_fluxstd_gaia_exclusive] = "Gaia-only"

    logger.info(f"Number of good fluxstd objects: {np.sum(good_fluxstd)}/{len(df_db)}")
    logger.info(f"  - PS1+Gaia: {np.sum(good_fluxstd_ps1_exclusive)}")
    logger.info(f"  - Gaia-only: {np.sum(good_fluxstd_gaia_exclusive)}")

    required_cols = [
        "ra",
        "dec",
        "epoch",
        "pmra",
        "pmra_error",
        "pmdec",
        "pmdec_error",
        "parallax",
        "parallax_error",
    ]
    missing_cols = [col for col in required_cols if col not in df_db.columns]
    if missing_cols:
        raise ValueError(
            f"Required columns missing from df_db: {missing_cols}. "
            f"Available columns: {list(df_db.columns)}"
        )

    # from "J2016.0" to 2016.0
    epoch_fluxstd = df_db["epoch"].str.replace("J", "", regex=False).astype(float)

    # Set flag for good proper motion and parallax measurements
    logger.info("Setting flags for good proper motion and parallax measurements")
    good_motion_fluxstd = (
        (np.abs(df_db["pmra"] / df_db["pmra_error"]) > DEFAULT_PM_SNR_THRESHOLD)
        & (np.abs(df_db["pmdec"] / df_db["pmdec_error"]) > DEFAULT_PM_SNR_THRESHOLD)
        & (df_db["parallax"] > 0.0)
        & (df_db["parallax"] / df_db["parallax_error"] > DEFAULT_PM_SNR_THRESHOLD)
    )

    # Set bad proper motion and parallax measurements to zero
    logger.info(
        "Setting bad proper motion and parallax measurements to zero and very small values, respectively"
    )
    parallax_fluxstd = df_db["parallax"].to_numpy()
    parallax_fluxstd[~good_motion_fluxstd] = DEFAULT_PARALLAX_FALLBACK
    pmra_fluxstd = df_db["pmra"].to_numpy()
    pmra_fluxstd[~good_motion_fluxstd] = 0.0
    pmdec_fluxstd = df_db["pmdec"].to_numpy()
    pmdec_fluxstd[~good_motion_fluxstd] = 0.0

    coords_db = SkyCoord(
        ra=df_db["ra"].to_numpy() * u.deg,
        dec=df_db["dec"].to_numpy() * u.deg,
        frame="icrs",
        distance=Distance(parallax=parallax_fluxstd * u.mas),
        pm_ra_cosdec=pmra_fluxstd * (u.mas / u.yr),
        pm_dec=pmdec_fluxstd * (u.mas / u.yr),
        obstime=Time(epoch_fluxstd, format="jyear", scale="tcb"),
    )

    logger.info("Applying space motion to fluxstd objects to epoch 2000.0")
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", category=ErfaWarning, message=".*distance overridden.*"
        )
        coords_db_epoch2000 = coords_db.apply_space_motion(
            Time(2000.0, format="jyear", scale="tcb")
        )

    return coords_db_epoch2000, good_fluxstd, catalog_source


def match_catalog(
    df: pd.DataFrame,
    df_fluxstd: pd.DataFrame,
    sep: u.Quantity = DEFAULT_MATCH_SEPARATION,
) -> pd.DataFrame:
    """
    Match input catalog with flux standard catalog using sky coordinates.

    Performs spatial cross-matching using astropy's match_coordinates_sky,
    applying proper motion and parallax corrections to the flux standard catalog.
    Prioritizes PS1+Gaia flux standards, using Gaia-only as fallback.

    Parameters
    ----------
    df : pd.DataFrame
        Input catalog DataFrame with 'ra', 'dec', 'obj_id' columns
    df_fluxstd : pd.DataFrame
        Flux standard catalog DataFrame
    sep : u.Quantity, default=1*u.arcsec
        Maximum separation for a match

    Returns
    -------
    pd.DataFrame
        DataFrame with matched objects containing columns:

        - obj_id: Input object ID
        - ra, dec: Input object coordinates
        - fluxstd_obj_id: Matched flux standard object ID
        - fluxstd_ra, fluxstd_dec: Matched flux standard coordinates
        - sep_arcsec: Separation in arcseconds
        - catalog_source: Catalog source ("PS1+Gaia" or "Gaia-only")

    Raises
    ------
    ValueError
        If required columns are missing or DataFrames are empty
    """
    # Validate inputs
    if df.empty:
        raise ValueError("Input DataFrame df is empty")

    if df_fluxstd.empty:
        raise ValueError("Flux standard DataFrame df_fluxstd is empty")

    required_input_cols = ["ra", "dec", "obj_id"]
    missing_cols = [col for col in required_input_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(
            f"Required columns missing from input df: {missing_cols}. "
            f"Available columns: {list(df.columns)}"
        )

    coords_input = SkyCoord(
        ra=df["ra"].to_numpy() * u.deg,
        dec=df["dec"].to_numpy() * u.deg,
        frame="icrs",
    )

    coords_db_epoch2000, good_fluxstd, catalog_source = generate_fluxstd_coords(
        df_fluxstd
    )

    good_coords = coords_db_epoch2000[good_fluxstd]
    good_df_fluxstd = df_fluxstd.loc[good_fluxstd].reset_index(drop=True)
    good_catalog_source = catalog_source[good_fluxstd].reset_index(drop=True)

    # The result object has indices and angular_separation attributes
    # The length of these attributes is the same as the length of coords_input
    res_coord_match = match_coordinates_sky(coords_input, good_coords, nthneighbor=1)

    is_matched = res_coord_match.angular_separation <= sep

    matched_fluxstd_indices = res_coord_match.indices_to_catalog[is_matched].astype(int)

    df_fluxstd_matched = good_df_fluxstd.iloc[matched_fluxstd_indices].copy()

    df_out = df.loc[is_matched, ["obj_id", "ra", "dec"]].copy().reset_index(drop=True)
    df_out["fluxstd_obj_id"] = df_fluxstd_matched["obj_id"].to_numpy()
    # df_out["fluxstd_ra"] = df_fluxstd_matched["ra"].to_numpy()
    # df_out["fluxstd_dec"] = df_fluxstd_matched["dec"].to_numpy()
    df_out["fluxstd_ra"] = good_coords.ra[matched_fluxstd_indices].deg
    df_out["fluxstd_dec"] = good_coords.dec[matched_fluxstd_indices].deg
    df_out["sep_arcsec"] = (
        res_coord_match.angular_separation[is_matched].to(u.arcsec).value
    )
    df_out["catalog_source"] = good_catalog_source.to_numpy()[matched_fluxstd_indices]

    df_out["obj_id_str"] = df_out["obj_id"].astype(str)
    df_out["fluxstd_obj_id_str"] = df_out["fluxstd_obj_id"].astype(str)

    logger.info(f"Matched {len(df_out)} objects within {sep.to(u.arcsec).value} arcsec")

    logger.info(f"Showing first 5 matches:\n{df_out.head(5)}")

    return df_out


def load_input_list(input_file: str) -> pd.DataFrame:
    """
    Load input catalog from CSV file.

    Parameters
    ----------
    input_file : str
        Path to input CSV file

    Returns
    -------
    pd.DataFrame
        DataFrame containing input catalog

    Raises
    ------
    FileNotFoundError
        If input_file does not exist
    """

    def check_bigint(value: Any) -> int:
        try:
            int_value = int(value)
            if math.isclose(int_value, float(value)):
                if (
                    int_value >= -9223372036854775808
                    and int_value <= 9223372036854775807
                ):
                    return int_value
                else:
                    raise ValueError(
                        f"Out of range value detected (-9223372036854775808, 9223372036854775807): {value}"
                    )
            else:
                raise ValueError(f"Non integer value detected: {value}")
        except ValueError as e:
            raise ValueError(f"{e}") from e
        except TypeError as e:
            raise TypeError(f"{e}") from e

    validated_input_file = _validate_path(input_file, "input_file")

    logger.info(f"Loading input catalog from {validated_input_file}")
    df = pd.read_csv(
        validated_input_file,
        dtype={"ra": float, "dec": float},
        converters={"obj_id": check_bigint},
    )

    return df


def _get_neighbor_pixels(hpx_index: int, nside: int, order: str = "ring") -> list[int]:
    """
    Get neighboring HEALPix pixels.

    Uses astropy_healpix.neighbours() function to find the 8 neighboring pixels.

    Parameters
    ----------
    hpx_index : int
        HEALPix index of the central pixel
    nside : int
        HEALPix NSIDE parameter
    order : str, default="ring"
        HEALPix ordering scheme

    Returns
    -------
    list[int]
        List of unique neighbor pixel indices (always 8 neighbors in ring scheme)
    """
    neighbor_indices = neighbours(hpx_index, nside=nside, order=order)
    # Convert numpy array to list of integers
    return [int(idx) for idx in neighbor_indices]


def _expand_with_neighbors(
    hpx_indices: list[int], nside: int, order: str = DEFAULT_HEALPIX_ORDER
) -> list[int]:
    """
    Expand HEALPix pixel list to include neighboring pixels.

    Parameters
    ----------
    hpx_indices : list[int]
        List of HEALPix indices
    nside : int
        HEALPix NSIDE parameter
    order : str, default=DEFAULT_HEALPIX_ORDER
        HEALPix ordering scheme

    Returns
    -------
    list[int]
        Unique sorted list of pixels including originals and neighbors
    """
    expanded = set(hpx_indices)

    for hpx in hpx_indices:
        neighbor_list = _get_neighbor_pixels(hpx, nside, order)
        expanded.update(neighbor_list)

    return sorted(expanded)


def _process_healpix_batch(
    df_input_batch: pd.DataFrame,
    hpx_indices_batch: list[int],
    duckdb_data_root: str,
    nside: int,
    sep: u.Quantity = DEFAULT_MATCH_SEPARATION,
) -> pd.DataFrame:
    """
    Process a batch of HEALPix pixels for flux standard matching.

    Parameters
    ----------
    df_input_batch : pd.DataFrame
        Input objects belonging to this batch of pixels
    hpx_indices_batch : list[int]
        List of HEALPix indices to query (may include neighbors)
    duckdb_data_root : str
        Root directory containing Parquet files
    nside : int
        HEALPix NSIDE parameter
    sep : u.Quantity, default=DEFAULT_MATCH_SEPARATION
        Maximum separation for matching

    Returns
    -------
    pd.DataFrame
        Matched results for this batch (empty DataFrame if no matches or on error)
    """
    try:
        # Fetch flux standards for all pixels in batch
        df_fluxstd = fetch_fluxstd(hpx_indices_batch, duckdb_data_root, nside=nside)

        if df_fluxstd.empty:
            return pd.DataFrame()

        # Perform matching for this batch's data
        df_matched = match_catalog(df_input_batch, df_fluxstd, sep=sep)

        return df_matched

    except Exception as e:
        logger.error(f"Error processing HEALPix batch: {e}")
        return pd.DataFrame()


def process_fluxstd_lookup(df_input: pd.DataFrame) -> pd.DataFrame:
    """
    Process flux standard star lookup and matching.

    Parameters
    ----------
    df_input : pd.DataFrame
        Input catalog DataFrame with 'ra', 'dec', 'obj_id' columns

    Returns
    -------
    pd.DataFrame
        DataFrame with matched objects
    """
    # Add HEALPix columns
    df_input = add_healpix_columns(df_input, nside=DEFAULT_NSIDE)

    # Get unique HEALPix indices to query
    hpx_indices = df_input[f"hpx{DEFAULT_NSIDE}"].unique().tolist()

    logger.info(f"Processing {len(hpx_indices)} unique HEALPix pixels")

    duckdb_data_root = os.getenv("DUCKDB_DATA_ROOT")
    if not duckdb_data_root:
        logger.error("DUCKDB_DATA_ROOT is not set")
        raise ValueError(
            "DUCKDB_DATA_ROOT is not set. Create a .env file (see .env.example) or export DUCKDB_DATA_ROOT."
        )

    # Create batches of pixels
    batch_size = DEFAULT_BATCH_SIZE
    num_batches = (len(hpx_indices) + batch_size - 1) // batch_size
    batches = [
        hpx_indices[i * batch_size : (i + 1) * batch_size] for i in range(num_batches)
    ]

    logger.info(
        f"Processing in {num_batches} batches (batch size: {batch_size}, "
        f"include neighbors: {DEFAULT_INCLUDE_NEIGHBORS})"
    )

    # Start overall timing
    overall_start_time = time.time()

    # Process each batch
    results = []
    for batch_idx, batch_pixels in enumerate(batches, start=1):
        # Start batch timing
        batch_start_time = time.time()

        # Optionally expand with neighbors
        if DEFAULT_INCLUDE_NEIGHBORS:
            batch_pixels_expanded = _expand_with_neighbors(
                batch_pixels, DEFAULT_NSIDE, DEFAULT_HEALPIX_ORDER
            )
            logger.info(
                f"Batch {batch_idx}/{num_batches}: {len(batch_pixels)} pixels "
                f"expanded to {len(batch_pixels_expanded)} (with neighbors)"
            )
        else:
            batch_pixels_expanded = batch_pixels
            logger.info(
                f"Batch {batch_idx}/{num_batches}: {len(batch_pixels)} pixels "
                f"(neighbors disabled)"
            )

        # Get input objects for original pixels in this batch
        mask = df_input[f"hpx{DEFAULT_NSIDE}"].isin(batch_pixels)
        df_input_batch = df_input.loc[mask]

        logger.info(f"  Processing {len(df_input_batch)} input objects")

        # Process batch
        df_batch_result = _process_healpix_batch(
            df_input_batch, batch_pixels_expanded, duckdb_data_root, DEFAULT_NSIDE
        )

        # Calculate batch elapsed time
        batch_elapsed_time = time.time() - batch_start_time

        if not df_batch_result.empty:
            results.append(df_batch_result)
            logger.info(
                f"  Found {len(df_batch_result)} matches in this batch "
                f"(elapsed time: {batch_elapsed_time:.2f} sec)"
            )
        else:
            logger.info(
                f"  No matches in this batch (elapsed time: {batch_elapsed_time:.2f} sec)"
            )

    # Combine results from all pixels
    if not results:
        logger.warning("No matches found across any HEALPix pixels")
        return pd.DataFrame(
            columns=[
                "obj_id",
                "ra",
                "dec",
                "fluxstd_obj_id",
                "fluxstd_ra",
                "fluxstd_dec",
                "sep_arcsec",
                "catalog_source",
                "obj_id_str",
                "fluxstd_obj_id_str",
            ]
        )

    df_out = pd.concat(results, ignore_index=True)

    # Calculate overall elapsed time
    overall_elapsed_time = time.time() - overall_start_time

    logger.info(
        f"Total matches: {len(df_out)} objects across {len(results)}/{num_batches} batches "
        f"({len(hpx_indices)} unique input pixels)"
    )
    logger.info(f"Total matching time: {overall_elapsed_time:.2f} sec")

    return df_out


if __name__ == "__main__":

    duckdb_data_root = "./duckdb_root/fluxstd_v3.4_nside32/release/"
    input_file = "./tmp/gaia_dr3_cosmos_r30arcmin.csv"

    df = load_input_list(input_file)
    df = add_healpix_columns(df, nside=DEFAULT_NSIDE)

    # Get unique HEALPix indices to query
    hpx_indices = df[f"hpx{DEFAULT_NSIDE}"].unique().tolist()
    df_fluxstd = fetch_fluxstd(hpx_indices, duckdb_data_root, nside=DEFAULT_NSIDE)
    df_out = match_catalog(df, df_fluxstd)
