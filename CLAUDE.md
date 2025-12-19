# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Python-based web application for matching astronomical catalogs with flux standard star catalogs from Gaia DR3 and Pan-STARRS DR2. The application uses HEALPix spatial indexing for efficient sky coordinate matching and provides a Panel-based web interface for users to upload target lists and retrieve matched flux standards.

The application supports dual catalog selection: PS1+Gaia flux standards are preferred, with Gaia-only standards used as fallback for regions without Pan-STARRS coverage (Dec < -30°). Each matched star includes a `catalog_source` field indicating its origin.

## Architecture

### Core Components

The codebase consists of two main files:

- **[utils.py](utils.py)**: Core catalog matching utilities using HEALPix spatial indexing
  - Implements spatial matching with proper motion and parallax corrections
  - Uses DuckDB to query Parquet-based flux standard catalog partitioned by HEALPix pixels
  - Dual catalog selection: PS1+Gaia (preferred) and Gaia-only (fallback)
  - Applies quality cuts based on stellar probability, magnitude ranges, and effective temperature
  - Tracks catalog source for each matched star (PS1+Gaia or Gaia-only)
  - Processes catalogs in batches of HEALPix pixels with optional neighbor inclusion

- **[app.py](app.py)**: Panel web application interface
  - File upload widget for CSV input (requires `obj_id`, `ra`, `dec` columns)
  - Interactive result table with Tabulator widget
  - CSV download functionality for matched results
  - Uses Material template with custom Geist font styling

### Data Flow

1. User uploads CSV file with target coordinates (RA, Dec, obj_id)
2. Input coordinates are assigned HEALPix indices at configured NSIDE (default: 32)
3. Flux standard catalog is queried from DuckDB/Parquet files by HEALPix pixel
4. Proper motion and parallax corrections are applied to flux standards (epoch 2000.0)
5. Quality cuts are applied (stellar probability, magnitude, temperature)
6. Sky coordinate matching performed using astropy's `match_coordinates_sky`
7. Results displayed in interactive table and available for CSV download

### Configuration

All runtime parameters are configurable via environment variables (see [.env.example](.env.example)):

- **DuckDB settings**: `DUCKDB_THREADS`, `DUCKDB_DATA_ROOT` (required)
- **HEALPix parameters**: `FLUXSTD_NSIDE`, `FLUXSTD_HEALPIX_ORDER`
- **Matching parameters**: `FLUXSTD_MATCH_SEPARATION_ARCSEC`, `FLUXSTD_PM_SNR_THRESHOLD`
- **Quality cuts**: Magnitude ranges, temperature ranges, stellar probability thresholds
- **Batch processing**: `FLUXSTD_BATCH_SIZE`, `FLUXSTD_INCLUDE_NEIGHBORS`

Configuration is loaded via `python-dotenv` and parsed with type-safe helper functions (`_env_int`, `_env_float`, `_env_str`, `_env_bool`) in [utils.py](utils.py:52-114).

### Key Algorithms

**HEALPix Spatial Indexing** ([utils.py](utils.py:217-299)):
- Assigns HEALPix indices to input coordinates using `astropy_healpix.lonlat_to_healpix`
- Supports both ring and nested ordering schemes
- Optionally includes neighboring pixels to handle boundary cases

**Proper Motion Correction** ([utils.py](utils.py:382-524)):
- Reads proper motion (pmra, pmdec), parallax, and epoch from catalog
- Filters out low-SNR measurements (default: SNR < 3.0)
- Falls back to zero proper motion and very small parallax for poor measurements
- Propagates coordinates to epoch 2000.0 using `SkyCoord.apply_space_motion`

**Quality Selection** ([utils.py](utils.py:433-475)):
- **Primary selection (PS1+Gaia)**: Pan-STARRS + Gaia with Brutus stellar probability and temperature
  - Uses `prob_f_star`, `psf_flux_g`, and `teff_brutus`
  - Preferred for Dec >= -30° where PS1 coverage exists
- **Fallback selection (Gaia-only)**: Gaia with GSP-Phot temperature and stellar flags
  - Uses `is_fstar_gaia`, `psf_flux_r`, and `teff_gspphot`
  - Used for Dec < -30° or when PS1 data unavailable
- Combined with OR operation: both selections contribute to final catalog
- Each star labeled with `catalog_source` ("PS1+Gaia" or "Gaia-only")

**Batch Processing** ([utils.py](utils.py:769-871)):
- Divides HEALPix pixels into batches (default: 20 pixels/batch)
- Processes each batch independently to manage memory
- Concatenates results from all batches

## Development Commands

### Environment Setup

This project uses `uv` for dependency management:

```bash
# Install dependencies (automatically creates/syncs virtual environment)
uv sync

# Add a new dependency
uv add <package-name>

# Add a development dependency
uv add --dev <package-name>
```

### Code Quality

```bash
# Format code with Black
uv run black .

# Lint with Ruff (includes auto-fixes)
uv run ruff check --fix .

# Type checking (if needed - not currently configured)
# uv run mypy .

# Run Python scripts/modules
uv run python script.py
uv run python -m module_name
```

### Running the Application

The application is launched via [launch_app.bash](launch_app.bash):

```bash
# Production mode (multi-threaded, port 5107)
./launch_app.bash

# Development mode (auto-reload, port 5207)
./launch_app.bash dev
```

**Requirements**:
- Create a `.env` file from `.env.example` with `DUCKDB_DATA_ROOT` and `PFS_APP_HOSTNAME` set
- The script validates that `PFS_APP_HOSTNAME` matches `hostname -f` before launching
- Static assets (fonts) are served from `./assets` directory

**Configuration Notes**:
- Production: 8 threads, port 5107, WebSocket max 150 MB
- Development: auto-reload enabled, port 5207, WebSocket max 150 MB
- URL prefix: `/fluxstd-lookup` (configurable via `PREFIX` in script)

### Testing

There are no automated tests currently. Manual testing workflow:

1. Prepare a CSV file with `obj_id`, `ra`, `dec` columns (see examples in `tmp/`)
2. Launch the app in dev mode: `./launch_app.bash dev`
3. Navigate to `http://localhost:5207/fluxstd-lookup`
4. Upload the CSV and verify matching results
5. Check logs for errors or warnings

## Important Implementation Notes

### Security

- **Path validation** ([utils.py](utils.py:180-214)): All file paths are validated with `_validate_path()` to prevent path traversal attacks
- **SQL injection prevention** ([utils.py](utils.py:339-365)): HEALPix indices are explicitly cast to integers before SQL query construction
- **No credentials in code**: Database connection strings and sensitive config are environment-only

### Error Handling

- Input validation with descriptive error messages for missing columns, empty DataFrames, invalid parameters
- Dual catalog support: PS1+Gaia and Gaia-only selections run in parallel, combined with OR operation
- If neither PS1+Gaia nor Gaia-only selections yield results, raises ValueError
- Batch processing errors are logged but don't halt entire job
- Panel UI shows error notifications with `pn.state.notifications.error()`

### Performance Considerations

- **Batch processing**: Processes multiple HEALPix pixels per batch to reduce overhead
- **Neighbor pixel inclusion**: Can be disabled via `FLUXSTD_INCLUDE_NEIGHBORS=false` if boundary matching not needed
- **DuckDB threading**: Set `DUCKDB_THREADS` based on available cores
- **HEALPix NSIDE**: Lower NSIDE = larger pixels = fewer queries but more data per query. Default 32 is balanced for sky coverage.

### Coordinate Systems

- Input coordinates: ICRS frame in degrees (J2000.0)
- Flux standard catalog: Gaia DR3 coordinates with proper motion at catalog epoch (typically J2016.0)
- Matching performed at: Epoch 2000.0 after proper motion correction
- Separation: Angular separation in arcseconds (configurable via `FLUXSTD_MATCH_SEPARATION_ARCSEC`)

## Code Style

- **Formatting**: Black (line length 88, target Python 3.12)
- **Linting**: Ruff with rules E4/E7/E9/F/B/I/UP/D (NumPy-style docstrings)
- **Docstrings**: All public functions have NumPy-style docstrings with Parameters, Returns, Raises sections
- **Type hints**: Used throughout (typing module for generics)
- **Logging**: loguru logger for all informational and error messages
- **Constants**: Module-level constants with environment variable overrides, documented with inline comments

## File References

When working with this codebase, be aware of these key locations:

- Flux standard catalog data: `$DUCKDB_DATA_ROOT/**/*.parquet` (HEALPix-partitioned)
- Configuration: `.env` (not in git, see [.env.example](.env.example))
- Launch script: [launch_app.bash](launch_app.bash)
- Static assets: `assets/fonts/` (Geist and Geist Mono font files)
- Example data: `tmp/*.csv` (not in git, for local testing)
- Migration script: [scripts/postgres2duckdb.py](scripts/postgres2duckdb.py) (converts PostgreSQL to Parquet)
