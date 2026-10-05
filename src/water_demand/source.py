"""NASA POWER daily gridded products with pinned scientific content and units."""

import hashlib
import json
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

UNITS = {
    "T2M": "C",
    "T2M_MAX": "C",
    "T2M_MIN": "C",
    "RH2M": "%",
    "WS2M": "m/s",
    "ALLSKY_SFC_SW_DWN": "MJ/m^2/day",
    "PRECTOTCORR": "mm/day",
    "PS": "kPa",
}


def scientific_hash(payload):
    # Response timing is volatile metadata, not part of the scientific data.
    stable = {k: v for k, v in payload.items() if k != "times"}
    return hashlib.sha256(
        json.dumps(stable, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def read_sources(root: Path):
    raw = root / "data/raw"
    raw.mkdir(parents=True, exist_ok=True)
    frames = []
    for item in json.loads((root / "data/source-manifest.json").read_text()):
        file = raw / item["file"]
        if not file.exists():
            with urllib.request.urlopen(item["url"], timeout=90) as source:
                payload = json.loads(source.read())
            if scientific_hash(payload) != item["scientific_sha256"]:
                raise ValueError(
                    "POWER scientific content changed; explicitly review source version"
                )
            temporary = file.with_suffix(".part")
            temporary.write_text(json.dumps(payload, separators=(",", ":")))
            temporary.replace(file)
        payload = json.loads(file.read_text())
        if scientific_hash(payload) != item["scientific_sha256"]:
            raise ValueError("POWER scientific content checksum mismatch")
        if payload["header"]["time_standard"] != "LST":
            raise ValueError("Local standard time required")
        for name, unit in UNITS.items():
            if payload["parameters"][name]["units"] != unit:
                raise ValueError("Unexpected parameter units: " + name)
        frame = pd.DataFrame(payload["properties"]["parameter"])
        frame.index = pd.to_datetime(frame.index, format="%Y%m%d")
        if frame.index.duplicated().any():
            raise ValueError("Duplicate daily dates")
        frame = frame.reindex(pd.date_range("2014-01-01", "2024-12-31", freq="D"))
        frame = frame.mask(frame == payload["header"]["fill_value"])
        frame = frame.where(np.isfinite(frame))
        frame["site"] = item["name"]
        frame["latitude"] = item["latitude"]
        frame["longitude"] = item["longitude"]
        frame["elevation_m"] = payload["geometry"]["coordinates"][2]
        frame.index.name = "date"
        frames.append(frame.reset_index())
    full = pd.concat(frames, ignore_index=True)
    if (
        not full.RH2M.dropna().between(0, 100).all()
        or (full.WS2M.dropna() < 0).any()
        or (full.PRECTOTCORR.dropna() < 0).any()
        or (full.ALLSKY_SFC_SW_DWN.dropna() < 0).any()
        or (full.PS.dropna() <= 0).any()
        or (full.T2M_MIN > full.T2M_MAX).any()
    ):
        raise ValueError("Invalid physical source range; review before analysis")
    return full
