#!/usr/bin/env python3

import io
import os
import tempfile

import numpy as np
import pandas as pd
import panel as pn
from dotenv import load_dotenv
from loguru import logger

from utils import (
    add_healpix_columns,
    fetch_fluxstd,
    generate_fluxstd_coords,
    load_input_list,
    match_catalog,
)

pn.extension("tabulator", notifications=True)

logger.info("PFS Flux Standard Star Lookup Tool")


load_dotenv()

# Widgets
file_input = pn.widgets.FileInput(
    accept=".csv",
    width=500,
    height=50,
    multiple=False,
)

process_button = pn.widgets.Button(
    name="Match Flux Standard Stars",
    button_type="primary",
    # button_type="default",
    # button_style="outline",
    icon="search",
    icon_size="1.5em",
    width=300,
    height=50,
)

# Download widget
download_button = pn.widgets.FileDownload(
    callback=lambda: None,  # Will be updated dynamically
    filename="matched_fluxstd.csv",
    button_type="primary",
    # button_style="outline",
    icon="download",
    icon_size="1.5em",
    width=300,
    height=50,
    label="Download CSV",
)
download_button.disabled = True

result_table = pn.widgets.Tabulator(
    name="Matched Flux Standard Stars",
    width=1280,
    height=640,
    pagination="remote",
    page_size=250,
    show_index=False,
    disabled=True,
    selectable=False,
    layout="fit_data_fill",
    titles={
        "obj_id_str": "obj_id",
        "ra": "RA (deg)",
        "dec": "Dec (deg)",
        "fluxstd_obj_id_str": "fluxstd_obj_id",
        "fluxstd_ra": "fluxstd_RA (deg)",
        "fluxstd_dec": "fluxstd_Dec (deg)",
        "sep_arcsec": "separation (arcsec)",
    },
    stylesheets=[
        """
        .tabulator-cell {
            font-family: 'Geist Mono', 'Courier New', monospace !important;
            /* Enables slashed zero */
            font-feature-settings: "zero" !important;
        }
        """
    ],
)
result_table.visible = False
process_button.disabled = True


# Callback
def match_fluxstd_callback(event):
    """Process the uploaded file and match flux standard stars."""
    process_button.disabled = True

    process_button.disabled = True

    result_table.value = pd.DataFrame()  # Clear previous results

    if not file_input.value:
        pn.state.notifications.error("Please upload a CSV file before processing.")
        return

    # Read the uploaded file
    # file_input.value is bytes
    filename = file_input.filename

    logger.info(f"Processing uploaded file: {filename}")

    try:
        # FileInput provides bytes content, so we write to a temporary file
        # and pass the path to load_input_list()
        with tempfile.NamedTemporaryFile(
            mode="wb", suffix=".csv", delete=True
        ) as tmp_file:
            tmp_file.write(file_input.value)
            tmp_file.flush()  # Ensure data is written to disk
            # Load using the existing utility function while file still exists
            df_input = load_input_list(tmp_file.name)

        # Temporary file is automatically deleted here
        logger.info(f"Successfully loaded {len(df_input)} rows from uploaded CSV")

        # Perform matching
        df_input = add_healpix_columns(df_input)
        hpx_indices = df_input["hpx32"].unique().tolist()

        duckdb_data_root = os.getenv("DUCKDB_DATA_ROOT")
        if not duckdb_data_root:
            pn.state.notifications.error(
                "DUCKDB_DATA_ROOT is not set. Create a .env file (see .env.example) or export DUCKDB_DATA_ROOT.",
                duration=0,
            )
            logger.error("DUCKDB_DATA_ROOT is not set")
            return

        df_fluxstd = fetch_fluxstd(hpx_indices, duckdb_data_root)
        df_out = match_catalog(df_input, df_fluxstd)

        # Update the result table
        result_table.value = df_out.loc[
            :,
            [
                "obj_id_str",
                "ra",
                "dec",
                "fluxstd_obj_id_str",
                "fluxstd_ra",
                "fluxstd_dec",
                "sep_arcsec",
            ],
        ]
        result_table.visible = True

        # Update download button with the results
        # NOTE: FileDownload treats a returned `str` as a file path.
        # Return bytes (or a file-like object) to download generated content.
        def get_csv_data():
            csv_bytes = (
                df_out.loc[
                    :,
                    [
                        "obj_id",
                        "ra",
                        "dec",
                        "fluxstd_obj_id",
                        "fluxstd_ra",
                        "fluxstd_dec",
                        "sep_arcsec",
                    ],
                ]
                .to_csv(index=False)
                .encode("utf-8")
            )
            return io.BytesIO(csv_bytes)

        download_button.callback = get_csv_data
        download_button.disabled = False

        pn.state.notifications.success("Flux standard star matching completed.")
    except Exception as e:
        pn.state.notifications.error(
            f"Error reading the uploaded file: {e}", duration=0
        )
        logger.error(f"Error reading uploaded file: {e}")
        return
    finally:
        process_button.disabled = file_input.value is None


def _toggle_button_on_file(event):
    process_button.disabled = event.new is None or len(event.new) == 0


# watch value
file_input.param.watch(_toggle_button_on_file, "value")
# connect callbacks
process_button.on_click(match_fluxstd_callback)


# Layout
main = pn.Column(
    pn.pane.Markdown(
        """
### Usage
- Prepare a CSV file containing your targets following [the format required for PFS Target Uploader](https://pfs-etc.naoj.hawaii.edu/uploader/doc/inputs.html).
- The minimum requirement is to have **`obj_id`**, **`ra`**, and **`dec`** columns.
- Press the **Match Flux Standard Stars** button to start the matching process to see the result.
- You can download the matched results as a CSV file by clicking the **Download CSV** button

### Note
- Flux standard stars are selected from Gaia DR3 catalog cross-matched with Pan-STARRS (PS1) DR2 when available.
- Coordinates are matched after applying proper motion correction to the epoch of 2000.0.
""",
        disable_anchors=True,
        styles={"font-size": "medium"},
    ),
    pn.Row(file_input, process_button, download_button),
    result_table,
)

# pn.template.FastListTemplate(
pn.template.MaterialTemplate(
    title="PFS Flux Standard Star Lookup",
    main=main,
    raw_css=[
        """
        /* Local Variable Fonts */
        @font-face {
            font-family: 'Geist';
            src: url('assets/fonts/Geist-Variable.woff2') format('woff2');
            font-weight: 100 900;
            font-style: normal;
            font-display: swap;
        }

        @font-face {
            font-family: 'Geist Mono';
            src: url('assets/fonts/GeistMono-Variable.woff2') format('woff2');
            font-weight: 100 900;
            font-style: normal;
            font-display: swap;
        }

        /* Global Reset & Override */
        :root, body, .bk, .bk-root, :host {
            font-family: 'Geist', sans-serif !important;
            --body-font: 'Geist', sans-serif !important;
            --bk-font-family: 'Geist', sans-serif !important;

            /* monospace numbers and better visibility (tnum, cv05) */
            /* font-feature-settings: "tnum", "cv05" !important; */
        }
        """
    ],
).servable()
