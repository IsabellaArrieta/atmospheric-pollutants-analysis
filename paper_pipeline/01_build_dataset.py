"""Stage 1: download EPA AirData daily files (2023-2025) and build clean daily site series.

Outputs (table_csv/):
    sites_meta.parquet        one row per AQS site with land use / location setting
    daily_long.parquet        one value per site x parameter x day (POCs averaged)
    screening_summary.csv     completeness counts per site and parameter
"""

import io
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (DATA_DIR, DURATION_PRIORITY, LAND_USE_ORDER, LOCATION_ORDER, MIN_OBS_PERCENT,
                    PARAMETERS, RAW_DIR, TAB_DIR, YEARS)

BASE_URL = "https://aqs.epa.gov/aqsweb/airdata"
# 50 states + DC (FIPS); territories, Mexico (80) and Canada (CC) are excluded.
US_STATE_FIPS = {f"{c:02d}" for c in [1, 2, 4, 5, 6, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21, 22,
                                     23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
                                     40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 53, 54, 55, 56]}


def download(name):
    target = RAW_DIR / name
    if not target.exists():
        print(f"downloading {name}")
        response = requests.get(f"{BASE_URL}/{name}", timeout=300)
        response.raise_for_status()
        target.write_bytes(response.content)
    return target


def read_zip_csv(path, usecols=None, dtype=None):
    with zipfile.ZipFile(path) as archive:
        member = next(n for n in archive.namelist() if n.lower().endswith(".csv"))
        with archive.open(member) as handle:
            return pd.read_csv(handle, usecols=usecols, dtype=dtype, low_memory=False)


def site_id(state, county, site):
    return (state.astype(str).str.zfill(2) + "-" + county.astype(str).str.zfill(3) + "-"
            + site.astype(str).str.zfill(4))


def build_sites():
    sites = read_zip_csv(download("aqs_sites.zip"), dtype={"State Code": str})
    sites["site_id"] = site_id(sites["State Code"], sites["County Code"], sites["Site Number"])
    sites["Land Use"] = sites["Land Use"].str.strip().str.capitalize()
    sites["Location Setting"] = sites["Location Setting"].str.strip().str.capitalize()
    sites["is_us"] = sites["State Code"].str.zfill(2).isin(US_STATE_FIPS)
    sites["is_open"] = sites["Site Closed Date"].isna()
    keep = ["site_id", "State Code", "State Name", "Latitude", "Longitude", "Land Use",
            "Location Setting", "is_us", "is_open"]
    return sites[keep].drop_duplicates("site_id")


def build_parameter(code):
    usecols = ["State Code", "County Code", "Site Num", "POC", "Sample Duration", "Date Local",
               "Event Type", "Observation Percent", "Arithmetic Mean"]
    frames = []
    for year in YEARS:
        frame = read_zip_csv(download(f"daily_{code}_{year}.zip"), usecols=usecols,
                             dtype={"State Code": str})
        frames.append(frame)
    raw = pd.concat(frames, ignore_index=True)
    stats = {"parameter": PARAMETERS[code], "raw_rows": len(raw)}
    raw["site_id"] = site_id(raw["State Code"], raw["County Code"], raw["Site Num"])
    raw["date"] = pd.to_datetime(raw["Date Local"])
    raw = raw[raw["date"].dt.year.isin(YEARS)]
    # Exceptional-event days appear twice ("Included"/"Excluded"); keep the measured concentration.
    raw = raw[raw["Event Type"].fillna("None").ne("Excluded")]
    raw = raw[raw["Observation Percent"].fillna(100) >= MIN_OBS_PERCENT]
    raw = raw[np.isfinite(raw["Arithmetic Mean"])]
    # Identical daily values are repeated once per pollutant standard: collapse them first.
    per_duration = (raw.groupby(["site_id", "POC", "Sample Duration", "date"], observed=True)
                    ["Arithmetic Mean"].mean().reset_index())
    rank = {name: i for i, name in enumerate(DURATION_PRIORITY)}
    per_duration["priority"] = per_duration["Sample Duration"].map(rank).fillna(len(rank))
    per_poc = (per_duration.sort_values("priority")
               .drop_duplicates(["site_id", "POC", "date"], keep="first"))
    # Averaging collocated instruments with different sampling schedules (e.g. 1-in-3-day filters
    # and continuous monitors) injects artificial day-to-day jumps, so each site keeps the single
    # most complete instrument record (POC); ties go to the lowest POC.
    poc_days = per_poc.groupby(["site_id", "POC"])["date"].nunique().rename("n").reset_index()
    primary = (poc_days.sort_values(["site_id", "n", "POC"], ascending=[True, False, True])
               .drop_duplicates("site_id"))
    daily = per_poc.merge(primary[["site_id", "POC"]], on=["site_id", "POC"], how="inner")
    daily = daily[["site_id", "POC", "date", "Arithmetic Mean"]].rename(columns={"Arithmetic Mean": "value"})
    daily["parameter"] = PARAMETERS[code]
    stats.update({
        "rows_after_qc": len(raw),
        "site_poc_days": len(per_poc),
        "site_days": len(daily),
        "sites": daily["site_id"].nunique(),
        "sites_with_multiple_poc": int((poc_days.groupby("site_id").size() > 1).sum()),
        "max_records_per_site_day_raw": int(raw.groupby(["site_id", "date"]).size().max()),
    })
    print(stats)
    return daily, stats


def main():
    sites = build_sites()
    sites.to_parquet(DATA_DIR / "sites_meta.parquet", index=False)

    daily_frames, stats = [], []
    for code in PARAMETERS:
        daily, info = build_parameter(code)
        daily_frames.append(daily)
        stats.append(info)
    daily_long = pd.concat(daily_frames, ignore_index=True)
    daily_long = daily_long.merge(sites, on="site_id", how="left")
    daily_long.to_parquet(DATA_DIR / "daily_long.parquet", index=False)
    pd.DataFrame(stats).to_csv(TAB_DIR / "build_stats.csv", index=False)

    valid_days = (daily_long.groupby(["site_id", "parameter"])["date"].nunique()
                  .rename("valid_days").reset_index())
    valid_days = valid_days.merge(sites, on="site_id", how="left")
    valid_days["in_scope"] = (valid_days["is_us"].fillna(False)
                              & valid_days["Land Use"].isin(LAND_USE_ORDER)
                              & valid_days["Location Setting"].isin(LOCATION_ORDER))
    valid_days.to_csv(TAB_DIR / "screening_summary.csv", index=False)
    print(valid_days.groupby("parameter")["valid_days"].describe())


if __name__ == "__main__":
    main()
