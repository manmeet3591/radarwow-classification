"""Quasi-vertical profiles (Ryzhkov et al. 2016, JTECH, doi:10.1175/JTECH-D-15-0020.1).

Averaging a high-elevation PPI over all azimuths turns it into one profile versus height. Ryzhkov et al.
recommend elevations of 10-30 deg; Sanchez-Rivas & Rico-Ramirez (2021) used 9 deg when nothing higher
was available. Each profile point averages gates at one range, i.e. one height and one horizontal
distance (range x cos(elevation)), so storm structure inside that ring is mixed into the profile.
"""

from __future__ import annotations

import numpy as np
import xarray as xr

from radarwow.core.noise import signal_mask

KE_RE_M = 4.0 / 3.0 * 6371e3


def beam_height_m(r_m, elevation_deg):
    """Beam-centre height above the radar for the 4/3-Earth model."""
    el = np.deg2rad(elevation_deg)
    return np.sqrt(r_m ** 2 + KE_RE_M ** 2 + 2 * r_m * KE_RE_M * np.sin(el)) - KE_RE_M


def sweep_elevation(ds) -> float:
    if "sweep_fixed_angle" in ds:
        return float(np.asarray(ds["sweep_fixed_angle"]).ravel()[0])
    return float(np.nanmedian(ds["elevation"].values))


def quasi_vertical_profile(ds: xr.Dataset, moments=("DBZH", "ZDR", "RHOHV", "PHIDP"), min_fraction: float = 0.3) -> xr.Dataset:
    """Azimuthal mean of each moment over signal gates, versus beam height (km above the radar).

    Heights where fewer than ``min_fraction`` of the azimuths have signal are set to NaN.
    PHIDP is averaged as given (offset not removed); only its shape is used downstream.
    """
    el = sweep_elevation(ds)
    sig = signal_mask(ds)
    frac = sig.mean(axis=0)
    h_km = beam_height_m(ds["range"].values, el) / 1e3
    out = {}
    for m in moments:
        if m not in ds:
            continue
        v = np.where(sig, ds[m].values.astype(float), np.nan)
        prof = np.nanmean(v, axis=0)
        out[m] = (("height",), np.where(frac >= min_fraction, prof, np.nan).astype("float32"))
    out["FRACTION"] = (("height",), frac.astype("float32"))
    return xr.Dataset(out, coords={"height": h_km, "range": ("height", ds["range"].values)},
                      attrs={"elevation_deg": el, "reference": "Ryzhkov et al. 2016, doi:10.1175/JTECH-D-15-0020.1"})


def highest_sweep(dt, min_elevation: float = 9.0):
    """Name and dataset of the highest-elevation sweep with elevation >= min_elevation (or None)."""
    best = None
    for name, node in dt.children.items():
        if not name.startswith("sweep"):
            continue
        ds = node.to_dataset()
        if "DBZH" not in ds or "RHOHV" not in ds:
            continue
        el = sweep_elevation(ds)
        if el >= min_elevation and (best is None or el > best[2]):
            best = (name, ds, el)
    return (best[0], best[1]) if best else (None, None)
