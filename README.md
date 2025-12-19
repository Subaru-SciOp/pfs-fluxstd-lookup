# PFS Flux Standard Star Lookup

A web application for matching astronomical target catalogs with flux standard star catalogs from Gaia DR3 and Pan-STARRS DR2.

## Overview

This tool helps astronomers identify suitable flux standard stars for spectroscopic observations with the Prime Focus Spectrograph (PFS) at Subaru Telescope. Given a list of target coordinates, it finds nearby flux standard stars with well-characterized photometry and astrometry.

The application uses HEALPix spatial indexing for efficient catalog queries and applies proper motion corrections to ensure accurate coordinate matching.

## Features

- Upload target lists in CSV format
- Match targets with flux standard stars from Gaia DR3 × Pan-STARRS DR2
- Dual catalog support: prioritize PS1+Gaia, fallback to Gaia-only for southern sky (Dec < -30°)
- Apply proper motion and parallax corrections to epoch 2000.0
- Quality filtering based on stellar probability, magnitude, and temperature
- Track catalog source for each matched star (PS1+Gaia or Gaia-only)
- Interactive result visualization
- Download matched results as CSV

## Requirements

- Python 3.12 or later
- uv (Python package manager)
- Flux standard catalog data in Parquet format (partitioned by HEALPix)
- PostgreSQL database with flux standard catalog (for initial data preparation)

## Installation

1. Clone the repository:

```bash
git clone https://github.com/Subaru-PFS/pfs-fluxstd-lookup.git
cd pfs-fluxstd-lookup
```

2. Install dependencies using uv:

```bash
uv sync
```

3. Create a `.env` file from the example:

```bash
cp .env.example .env
```

4. Edit `.env` and configure the required variables:

```bash
# Required: Path to flux standard catalog data
DUCKDB_DATA_ROOT="/path/to/fluxstd_catalog/parquet/"

# Required: Hostname where the app will run
PFS_APP_HOSTNAME="your-hostname.example.org"
```

## Data Preparation

Before running the application, you need to prepare the flux standard catalog data in HEALPix-partitioned Parquet format.

### Exporting from PostgreSQL

If your flux standard catalog is stored in PostgreSQL, use the provided migration script:

```bash
# Create a TOML configuration file with database credentials
cat > config.toml << 'EOF'
[targetdb.db]
host = "your-db-host"
port = 5432
dbname = "your-database"
user = "your-username"
password = "your-password"

[targetdb.fluxstd]
# Optional: Filter by input_catalog_id or version
# input_catalog_id = [1, 2, 3]
# version = ["v3.4"]
EOF

# Export to HEALPix-partitioned Parquet
uv run python scripts/postgres2duckdb.py \
  --config config.toml \
  --outdir ./duckdb_root/fluxstd_v3.4_nside32/ \
  --nside 32 \
  --threads 16 \
  --clean
```

This will create a directory structure like:

```
./duckdb_root/fluxstd_v3.4_nside32/
└── release/
    ├── hpx32_group=0/
    │   └── data.parquet
    ├── hpx32_group=1/
    │   └── data.parquet
    ...
```

### Export Options

Key parameters for `postgres2duckdb.py`:

- `--config`: TOML file with database credentials (required)
- `--outdir`: Output directory for Parquet files (required)
- `--nside`: HEALPix NSIDE parameter (default: 32)
- `--group-size`: HEALPix grouping factor (default: 64, groups 64 pixels per partition)
- `--threads`: DuckDB threads for parallel processing (default: 16)
- `--chunk-rows`: Rows per chunk from PostgreSQL (default: 2,000,000)
- `--clean`: Remove temporary directory after compaction
- `--compact-only`: Only perform compaction on existing tmp data (skip PostgreSQL fetch)

The script performs two phases:
1. **Export**: Fetches data from PostgreSQL in chunks and writes to `tmp/` directory
2. **Compaction**: Consolidates multiple files per partition into single files in `release/` directory

### Setting DUCKDB_DATA_ROOT

After export, update your `.env` file to point to the release directory:

```bash
DUCKDB_DATA_ROOT="/absolute/path/to/duckdb_root/fluxstd_v3.4_nside32/release/"
```

## Usage

### Running the Application

Start the application in production mode:

```bash
./launch_app.bash
```

Or in development mode with auto-reload:

```bash
./launch_app.bash dev
```

The application will be available at:
- Production: `http://your-hostname:5107/fluxstd-lookup`
- Development: `http://localhost:5207/fluxstd-lookup`

### Input Format

Upload a CSV file containing your target list with the following columns:

- `obj_id`: Unique identifier for each object (integer)
- `ra`: Right Ascension in degrees (float)
- `dec`: Declination in degrees (float)

Example:

```csv
obj_id,ra,dec
1,150.1234,2.5678
2,150.2345,2.6789
```

The format follows the [PFS Target Uploader specification](https://pfs-etc.naoj.hawaii.edu/uploader/doc/inputs.html).

### Output Format

The matched results include:

- `obj_id`: Original target ID
- `ra`, `dec`: Original target coordinates (degrees)
- `fluxstd_obj_id`: Matched flux standard star ID (Gaia source_id)
- `fluxstd_ra`, `fluxstd_dec`: Matched flux standard coordinates at epoch 2000.0 (degrees)
- `sep_arcsec`: Angular separation between target and flux standard (arcseconds)
- `catalog_source`: Catalog source ("PS1+Gaia" or "Gaia-only")

## Configuration

All parameters can be customized via environment variables in the `.env` file:

### Matching Parameters

- `FLUXSTD_MATCH_SEPARATION_ARCSEC`: Maximum separation for a match (default: 1.0)
- `FLUXSTD_NSIDE`: HEALPix NSIDE parameter for spatial indexing (default: 32)
- `FLUXSTD_PM_SNR_THRESHOLD`: SNR threshold for proper motion acceptance (default: 3.0)

### Quality Cuts

**PS1+Gaia selection (preferred)**:
- `FLUXSTD_MIN_PSF_MAG_G`, `FLUXSTD_MAX_PSF_MAG_G`: Pan-STARRS g-band magnitude range (default: 17.0-19.0)
- `FLUXSTD_MIN_TEFF_BRUTUS`, `FLUXSTD_MAX_TEFF_BRUTUS`: Effective temperature range from Brutus (default: 5000-8000 K)
- `FLUXSTD_MIN_PROB_F_STAR`: Minimum stellar probability (default: 0.5)

**Gaia-only selection (fallback)**:
- `FLUXSTD_MIN_PSF_MAG_R`, `FLUXSTD_MAX_PSF_MAG_R`: Gaia G-band magnitude range (default: 17.0-19.0)
- `FLUXSTD_MIN_TEFF_GSPPHOT`, `FLUXSTD_MAX_TEFF_GSPPHOT`: Effective temperature range from GSP-Phot (default: 5000-8000 K)

### Performance

- `DUCKDB_THREADS`: Number of threads for DuckDB queries (default: 16)
- `FLUXSTD_BATCH_SIZE`: Number of HEALPix pixels per batch (default: 20)
- `FLUXSTD_INCLUDE_NEIGHBORS`: Include neighboring pixels in search (default: true)

See [.env.example](.env.example) for complete configuration options.

## Development

### Code Quality

Format code with Black:

```bash
uv run black .
```

Lint with Ruff:

```bash
uv run ruff check --fix .
```

### Project Structure

```
.
├── app.py                  # Panel web application
├── utils.py                # Core catalog matching utilities
├── launch_app.bash         # Application launcher script
├── .env.example            # Configuration template
├── pyproject.toml          # Project metadata and dependencies
├── assets/fonts/           # Web fonts (Geist)
├── scripts/                # Utility scripts
│   └── postgres2duckdb.py  # Catalog migration tool
└── tmp/                    # Example data (not in git)
```

## Technical Details

### Catalog Matching Process

1. Input coordinates are assigned HEALPix pixel indices
2. Flux standard catalog is queried from Parquet files by HEALPix pixel
3. Proper motion and parallax corrections are applied to flux standards
4. Dual quality selection applied:
   - **PS1+Gaia**: Uses Brutus stellar probability and temperature (preferred)
   - **Gaia-only**: Uses GSP-Phot temperature and Gaia stellar flags (fallback for Dec < -30°)
5. Catalog source tracked for each star (PS1+Gaia or Gaia-only)
6. Angular separation matching using astropy's `match_coordinates_sky`
7. Results include only matches within the specified separation threshold

### Coordinate Systems

- Input coordinates: ICRS frame (J2000.0) in degrees
- Flux standard catalog: Gaia DR3 coordinates at catalog epoch (~J2016.0)
- Matching performed: Epoch 2000.0 after proper motion correction
- Proper motion corrections applied when SNR > 3.0 for pmra, pmdec, and parallax

### Data Requirements

The flux standard catalog must be stored as Parquet files partitioned by HEALPix pixel:

```
$DUCKDB_DATA_ROOT/
├── hpx32=0/
│   └── data.parquet
├── hpx32=1/
│   └── data.parquet
...
```

Each Parquet file should contain columns for:
- Gaia DR3 astrometry: ra, dec, pmra, pmdec, parallax, epoch
- Photometry: psf_flux_g (PS1 g-band), psf_flux_r (Gaia G-band)
- Stellar parameters: teff_brutus (Brutus), teff_gspphot (GSP-Phot), prob_f_star (Brutus stellar probability), is_fstar_gaia (Gaia stellar flag)

## License

MIT License. See [LICENSE](LICENSE) for details.

## Authors

- Masato Onodera (monodera@naoj.org)
- Subaru/PFS obsproc team (pfs-obs-help@naoj.org)

## Acknowledgments

This tool uses data from:
- Gaia DR3 (ESA/Gaia/DPAC)
- Pan-STARRS DR2 (PS1 Science Consortium)
