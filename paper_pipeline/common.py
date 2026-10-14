"""Shared paths, constants and plotting style for the manuscript pipeline.

Run the stages from the repository root, in order:
    python paper_pipeline/01_build_dataset.py
    python paper_pipeline/02_fdi_ifs.py
    python paper_pipeline/03_complexity.py
    python paper_pipeline/04_link_and_maps.py
"""

from pathlib import Path

import matplotlib as mpl

ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "table_csv" / "raw"
DATA_DIR = ROOT / "table_csv"
OUT_DIR = ROOT / "outputs_paper"
TAB_DIR = OUT_DIR / "tables"
FIG_DIR = OUT_DIR / "figures"
PAPER_FIG_DIR = ROOT / "docs" / "Clasiffication_Pollution_Sources_EPA" / "figs"
SHAPEFILE = ROOT / "paper_pipeline" / "shapes" / "cb_2018_us_state_5m.shp"

for folder in (RAW_DIR, OUT_DIR, TAB_DIR, FIG_DIR, PAPER_FIG_DIR):
    folder.mkdir(parents=True, exist_ok=True)

YEARS = [2023, 2024, 2025]
TOTAL_DAYS = 365 + 366 + 365  # 2024 is a leap year
COMPLETENESS = 0.75  # fraction of calendar days required per station-parameter
MIN_VALID_DAYS = int(round(COMPLETENESS * TOTAL_DAYS))  # 822 days
MIN_OBS_PERCENT = 75.0  # AQS daily records built from <75 % of sub-daily values are discarded

PARAMETERS = {
    "88101": "PM2.5 FRM/FEM",
    "88502": "PM2.5 non-FRM/FEM",
    "81102": "PM10",
    "86101": "PMc",
}
POLLUTANT_ORDER = ["PM2.5 FRM/FEM", "PM2.5 non-FRM/FEM", "PM10", "PMc"]
POLLUTANT_LABEL = {
    "PM2.5 FRM/FEM": r"PM$_{2.5}$ FRM/FEM",
    "PM2.5 non-FRM/FEM": r"PM$_{2.5}$ non-FRM/FEM",
    "PM10": r"PM$_{10}$",
    "PMc": r"PM$_{c}$",
}
# Within each site/POC/day we keep a single record, preferring the daily block average.
DURATION_PRIORITY = ["24-HR BLK AVG", "24 HOUR", "1 HOUR"]

LAND_USE_ORDER = [
    "Commercial", "Residential", "Industrial", "Agricultural",
    "Desert", "Military reservation", "Forest",
]
LOCATION_ORDER = ["Rural", "Suburban", "Urban and center city"]
SEASON_ORDER = ["Winter", "Spring", "Summer", "Autumn"]
SEASON_OF_MONTH = {12: "Winter", 1: "Winter", 2: "Winter", 3: "Spring", 4: "Spring", 5: "Spring",
                   6: "Summer", 7: "Summer", 8: "Summer", 9: "Autumn", 10: "Autumn", 11: "Autumn"}

# Colours kept from the original manuscript figures (one per land-use class).
LAND_USE_COLORS = {
    "Commercial": "rebeccapurple",
    "Residential": "red",
    "Industrial": "blue",
    "Agricultural": "darkgreen",
    "Desert": "darkgoldenrod",
    "Military reservation": "black",
    "Forest": "forestgreen",
}
LAND_USE_MARKERS = {
    "Commercial": ">", "Residential": "s", "Industrial": "*", "Agricultural": "^",
    "Desert": "D", "Military reservation": "o", "Forest": "o",
}
LOCATION_COLORS = {"Rural": "#1b9e77", "Suburban": "#d95f02", "Urban and center city": "#7570b3"}
POLLUTANT_COLORS = {
    "PM2.5 FRM/FEM": "#c0392b",
    "PM2.5 non-FRM/FEM": "#e67e22",
    "PM10": "#2c3e50",
    "PMc": "#b7950b",
}

SEED = 20260513


def set_style():
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 9.5,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "figure.dpi": 110,
    })


def save_figure(fig, name):
    """Save a figure as PNG in outputs_paper/figures and as PDF next to the manuscript."""
    fig.savefig(FIG_DIR / f"{name}.png")
    fig.savefig(PAPER_FIG_DIR / f"{name}.pdf")
