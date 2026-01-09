#!/usr/bin/env python3

import io
import os
import tempfile
from typing import Any

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
    process_fluxstd_lookup,
)

load_dotenv()

PRIMARY_COLOR = os.environ.get("PRIMARY_COLOR", "#5a79ba")
pn.extension(
    "tabulator",
    design="material",
    notifications=True,
    global_css=[f":root {{ --design-primary-color: {PRIMARY_COLOR}; }}"],
)

logger.info("PFS Flux Standard Star Lookup Tool")


# Widgets
file_input = pn.widgets.FileInput(
    accept=".csv",
    width=300,
    height=50,
    multiple=False,
)

process_button = pn.widgets.Button(
    name="Match Flux Standard Stars",
    button_type="primary",
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
    height=720,
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
        "catalog_source": "catalog",
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

# Object count display
object_count_text = pn.pane.Markdown(
    "",
    styles={"font-size": "large", "font-weight": "bold"},
    width=1280,
    # margin=(10, 0, 10, 0),
)
object_count_text.visible = False

# Loading spinner for processing state
loading_spinner_size = 120
loading_spinner = pn.indicators.LoadingSpinner(
    value=True,
    name="Processing...",
    size=loading_spinner_size,
    color="secondary",
    bgcolor="light",
    margin=(
        (720 - loading_spinner_size) // 8,
        0,
        0,
        (1280 - loading_spinner_size) // 10,
    ),
)

# Container to switch between spinner and result table
result_area = pn.Column(width=1280, height=720, visible=False, margin=(0, 0, 0, 0))


# Callback
def match_fluxstd_callback(event: Any) -> None:
    """Process the uploaded file and match flux standard stars."""
    process_button.disabled = True

    result_table.value = pd.DataFrame()  # Clear previous results

    if not file_input.value:
        pn.state.notifications.error("Please upload a CSV file before processing.")
        return

    # Show loading spinner
    result_area.objects = [loading_spinner]
    result_area.visible = True

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
        df_out = process_fluxstd_lookup(df_input)

        # Count matched objects
        num_objects = len(df_out)
        object_count_text.object = f"**Detected objects: {num_objects}**"
        object_count_text.visible = True

        # Update the result table and decide whether to show it
        if num_objects > 0:
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
                    "catalog_source",
                ],
            ]
            result_table.visible = True
            # Show both count and table
            result_area.objects = [object_count_text, result_table]
        else:
            # No matches - show only the count, hide table
            result_table.visible = False
            result_area.objects = [object_count_text]

        # Update download button with the results
        # NOTE: FileDownload treats a returned `str` as a file path.
        # Return bytes (or a file-like object) to download generated content.
        def get_csv_data() -> io.BytesIO:
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
                        "catalog_source",
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
            f"Error processing flux standard matching: {e}", duration=0
        )
        logger.error(f"Error during flux standard matching: {e}")
        # Hide result area on error
        result_area.visible = False
        return
    finally:
        process_button.disabled = file_input.value is None


def _toggle_button_on_file(event: Any) -> None:
    process_button.disabled = event.new is None or len(event.new) == 0
    # Clear previous results when new file is uploaded
    if event.new is not None and len(event.new) > 0:
        result_table.value = pd.DataFrame()
        result_table.visible = False
        object_count_text.visible = False
        result_area.visible = False
        download_button.disabled = True


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
- Cross-matching takes long time for objects distributing over large sky area. It is strongly recommended to limit the input list to objects within a small sky area (e.g., several sq. degrees) for efficient processing.
""",
        disable_anchors=True,
        styles={"font-size": "medium"},
        width=960,
    ),
    pn.Row(file_input, process_button, download_button),
    result_area,
)

# pn.template.FastListTemplate(
pn.template.MaterialTemplate(
    title="PFS Flux Standard Star Lookup",
    main=[main],
    # header_background="#5a79ba",
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
        * {
            font-family: 'Geist', sans-serif !important;
        }

        :root {
            --body-font: 'Geist', sans-serif !important;
            --bk-font-family: 'Geist', sans-serif !important;
        }

        body, .bk, .bk-root, :host,
        .mdc-typography, .mdc-typography--body1, .mdc-typography--body2,
        .markdown, .bk-input, button, input, select, textarea {
            font-family: 'Geist', sans-serif !important;
        }

        /* monospace numbers and better visibility (tnum, cv05) */
        /* font-feature-settings: "tnum", "cv05" !important; */
        """
    ],
).servable()
