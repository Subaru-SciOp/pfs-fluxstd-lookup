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
import warnings
from pathlib import Path

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
from astropy_healpix import lonlat_to_healpix
from erfa import ErfaWarning
from loguru import logger

# Configuration constants
DEFAULT_DUCKDB_THREADS = 16
"""Default number of threads for DuckDB operations"""

DEFAULT_PARALLAX_FALLBACK = 1e-7
"""Default parallax value (mas) for objects with poor parallax measurements.
This very small but non-zero value prevents division by zero while maintaining
the correct sign for distance calculations."""

DEFAULT_PM_SNR_THRESHOLD = 3.0
"""Signal-to-noise ratio threshold for accepting proper motion and parallax measurements"""

DEFAULT_MATCH_SEPARATION = 1 * u.arcsec
"""Default maximum separation for catalog matching"""

DEFAULT_NSIDE = 32
"""Default HEALPix NSIDE parameter"""

DEFAULT_HEALPIX_ORDER = "ring"
"""Default HEALPix ordering scheme"""

# Flux standard quality selection criteria
DEFAULT_MIN_PROB_F_STAR = 0.5
"""Minimum stellar probability for flux standard selection"""

DEFAULT_MAX_PROB_F_STAR = 1.0
"""Maximum stellar probability for flux standard selection"""

DEFAULT_MIN_PSF_MAG_G = 17.0
"""Minimum PS1 g-band PSF magnitude for flux standard selection"""

DEFAULT_MAX_PSF_MAG_G = 19.0
"""Maximum PS1 g-band PSF magnitude for flux standard selection"""

DEFAULT_MIN_PSF_MAG_R = 17.0
"""Minimum Gaia G-filter PSF magnitude for flux standard selection"""

DEFAULT_MAX_PSF_MAG_R = 19.0
"""Maximum Gaia G-filter PSF magnitude for flux standard selection"""

DEFAULT_MIN_TEFF_BRUTUS = 5000.0
"""Minimum effective temperature (K) from Brutus for flux standard selection"""

DEFAULT_MAX_TEFF_BRUTUS = 8000.0
"""Maximum effective temperature (K) from Brutus for flux standard selection"""

DEFAULT_MIN_TEFF_GSPPHOT = 5000.0
"""Minimum effective temperature (K) from Gaia GSP-Phot for flux standard selection"""

DEFAULT_MAX_TEFF_GSPPHOT = 8000.0
"""Maximum effective temperature (K) from Gaia GSP-Phot for flux standard selection"""


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


def generate_fluxstd_coords(df_db: pd.DataFrame) -> SkyCoord:
    """
    Generate SkyCoord objects for flux standard catalog with proper motion and parallax.

    Applies proper motion and parallax corrections, propagating coordinates to epoch 2000.0.
    Bad proper motion and parallax measurements (SNR < 3.0) are set to zero/default values.

    Parameters
    ----------
    df_db : pd.DataFrame
        DataFrame with flux standard catalog data

    Returns
    -------
    SkyCoord
        SkyCoord object with coordinates propagated to epoch 2000.0

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
    good_fluxstd = good_fstar_brutus & good_flux_ps1 & good_teff_brutus

    if not good_fluxstd.any():
        logger.warning(
            "No good flux PS1-Gaia stars found after applying quality cuts. Move to Gaia-only selection"
        )
        good_fstar_gaia = df_db["is_fstar_gaia"]
        good_teff_gaia = (df_db["teff_gspphot"] >= min_teff_gspphot) & (
            df_db["teff_gspphot"] <= max_teff_gspphot
        )
        good_flux_gaia = (df_db["psf_flux_r"] >= min_psf_flux_r) & (
            df_db["psf_flux_r"] <= max_psf_flux_r
        )
        good_fluxstd = good_fstar_gaia & good_teff_gaia & good_flux_gaia

        if not good_fluxstd.any():
            raise ValueError(
                "No good flux standard stars found after applying Gaia-only quality cuts"
            )

    logger.info(f"Number of good fluxstd objects: {np.sum(good_fluxstd)}/{len(df_db)}")

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

    return coords_db_epoch2000[good_fluxstd]


def match_catalog(
    df: pd.DataFrame,
    df_fluxstd: pd.DataFrame,
    sep: u.Quantity = DEFAULT_MATCH_SEPARATION,
) -> pd.DataFrame:
    """
    Match input catalog with flux standard catalog using sky coordinates.

    Performs spatial cross-matching using astropy's match_coordinates_sky,
    applying proper motion and parallax corrections to the flux standard catalog.

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

    coords_db_epoch2000 = generate_fluxstd_coords(df_fluxstd)

    # The result object has indices and angular_separation attributes
    # The length of these attributes is the same as the length of coords_input
    res_coord_match = match_coordinates_sky(
        coords_input, coords_db_epoch2000, nthneighbor=1
    )

    is_matched = res_coord_match.angular_separation <= sep

    matched_fluxstd_indices = res_coord_match.indices_to_catalog[is_matched]

    df_fluxstd_matched = df_fluxstd.iloc[matched_fluxstd_indices].copy()

    df_out = df.loc[is_matched, ["obj_id", "ra", "dec"]].copy().reset_index(drop=True)
    df_out["fluxstd_obj_id"] = df_fluxstd_matched["obj_id"].to_numpy()
    df_out["fluxstd_ra"] = df_fluxstd_matched["ra"].to_numpy()
    df_out["fluxstd_dec"] = df_fluxstd_matched["dec"].to_numpy()
    df_out["sep_arcsec"] = (
        res_coord_match.angular_separation[is_matched].to(u.arcsec).value
    )

    logger.info(f"Matched {len(df_out)} objects within {sep.to(u.arcsec).value} arcsec")

    logger.info(f"\n{df_out.head()}")

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

    def check_bigint(value):
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


if __name__ == "__main__":

    duckdb_data_root = "./duckdb_root/fluxstd_v3.4_nside32/release/"
    input_file = "./tmp/gaia_dr3_cosmos_r30arcmin.csv"

    df = load_input_list(input_file)
    df = add_healpix_columns(df, nside=DEFAULT_NSIDE)

    # Get unique HEALPix indices to query
    hpx_indices = df[f"hpx{DEFAULT_NSIDE}"].unique().tolist()
    df_fluxstd = fetch_fluxstd(hpx_indices, duckdb_data_root, nside=DEFAULT_NSIDE)
    df_out = match_catalog(df, df_fluxstd)
