"""Melting-layer detection from QVPs: Sanchez-Rivas & Rico-Ramirez 2021, AMT 14, 2873-2890,
doi:10.5194/amt-14-2873-2021 (QVP variant of their melting level algorithm, MLA).

Implemented as described in their Sect. 4 (read from the open-access article):
Part 1  Z* = Z normalised from [5, 60] dBZ to [0, 1]; rho* = rhoHV normalised from [0.85, 1] (values outside
        clipped). P_comb = Z* (1 - rho*). Its largest local maximum is P_peak; if P_peak < k (0.08 for QVPs)
        the profile has no melting-layer signature. Otherwise the upper limit UL = height(P_peak) + 750 m.
Part 2  Profiles are cut at UL; ZDR is min-max normalised within the profile. The combination they selected
        for QVPs is P14* = Z* ZDR* (1 - rho*). It is sharpened with P = P* - w P*'' (w = 0.75), the largest
        peak is checked against k again, and the valleys directly below and above it are the ML bottom and
        top (ML top = estimated melting level). Their reported error vs radiosondes is ~200 m.
They limit the search to heights below 5 km to reduce beam broadening (UK climate); we keep that as the
parameter ``max_height_km`` because the melting layer in the tropics and in Indian monsoon reaches ~5 km.
The second derivative is taken with respect to the sample index (the paper does not state the unit).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .qvp import highest_sweep, quasi_vertical_profile


def _norm(x, lo, hi):
    return np.clip((x - lo) / (hi - lo), 0.0, 1.0)


def _peaks(p):
    """Indices of local maxima (direct neighbours smaller), as in the paper's step 1b."""
    i = np.arange(1, p.size - 1)
    return i[(p[i] > p[i - 1]) & (p[i] > p[i + 1])]


def _valleys(p):
    i = np.arange(1, p.size - 1)
    return i[(p[i] < p[i - 1]) & (p[i] < p[i + 1])]


@dataclass
class MeltingLayer:
    detected: bool
    bottom_km: float = np.nan
    top_km: float = np.nan
    peak_km: float = np.nan
    p_peak: float = np.nan
    elevation_deg: float = np.nan
    reason: str = ""


def detect_melting_layer(qvp, k: float = 0.08, w: float = 0.75, ul_margin_km: float = 0.75,
                         max_height_km: float = 6.0, min_height_km: float = 0.2) -> MeltingLayer:
    """Apply the Sanchez-Rivas & Rico-Ramirez (2021) QVP algorithm to one profile (heights above the radar)."""
    h = qvp["height"].values
    keep = (h >= min_height_km) & (h <= max_height_km) & np.isfinite(qvp["DBZH"].values) & np.isfinite(qvp["RHOHV"].values)
    if keep.sum() < 10:
        return MeltingLayer(False, reason="too few valid heights in the profile")
    h = h[keep]
    z = _norm(qvp["DBZH"].values[keep], 5.0, 60.0)
    rho = _norm(qvp["RHOHV"].values[keep], 0.85, 1.0)
    pcomb = z * (1 - rho)
    pk = _peaks(pcomb)
    if pk.size == 0 or pcomb[pk].max() < k:
        return MeltingLayer(False, p_peak=float(pcomb[pk].max()) if pk.size else 0.0, elevation_deg=qvp.attrs.get("elevation_deg", np.nan),
                            reason=f"P_peak below k = {k}: no melting-layer signature")
    ipk = pk[np.argmax(pcomb[pk])]
    cut = h <= h[ipk] + ul_margin_km
    hh, zz, rr = h[cut], z[cut], rho[cut]
    if "ZDR" in qvp and np.isfinite(qvp["ZDR"].values[keep][cut]).sum() > 5:
        zdr = qvp["ZDR"].values[keep][cut]
        zdr = (zdr - np.nanmin(zdr)) / max(np.nanmax(zdr) - np.nanmin(zdr), 1e-6)
        pstar = zz * np.nan_to_num(zdr) * (1 - rr)
    else:
        pstar = zz * (1 - rr)  # Part 1 combination if ZDR is unavailable
    p = pstar - w * np.gradient(np.gradient(pstar))
    pk2 = _peaks(p)
    if pk2.size == 0 or p[pk2].max() < k:
        return MeltingLayer(False, p_peak=float(p[pk2].max()) if pk2.size else 0.0, elevation_deg=qvp.attrs.get("elevation_deg", np.nan),
                            reason="refined profile has no peak above k")
    i2 = pk2[np.argmax(p[pk2])]
    vl = _valleys(p)
    below, above = vl[vl < i2], vl[vl > i2]
    if below.size == 0 or above.size == 0:
        return MeltingLayer(False, peak_km=float(hh[i2]), p_peak=float(p[i2]), elevation_deg=qvp.attrs.get("elevation_deg", np.nan),
                            reason="no valley below or above the peak")
    return MeltingLayer(True, bottom_km=float(hh[below[-1]]), top_km=float(hh[above[0]]), peak_km=float(hh[i2]),
                        p_peak=float(p[i2]), elevation_deg=qvp.attrs.get("elevation_deg", np.nan))


def melting_layer_from_volume(dt, min_elevation: float = 9.0, **kw):
    """QVPs from sweeps at or above ``min_elevation``, highest first; the first one with a detected melting
    layer wins. Returns (MeltingLayer, qvp of the sweep used or of the highest sweep if none detected)."""
    from .qvp import sweep_elevation

    sweeps = []
    for name, node in dt.children.items():
        ds = node.to_dataset()
        if name.startswith("sweep") and "DBZH" in ds and "RHOHV" in ds:
            el = sweep_elevation(ds)
            if el >= min_elevation:
                sweeps.append((el, name, ds))
    if not sweeps:
        return MeltingLayer(False, reason=f"no sweep at or above {min_elevation} deg with DBZH and RHOHV"), None
    first = None
    seen = set()
    for el, name, ds in sorted(sweeps, key=lambda x: -x[0]):
        if round(el, 1) in seen:  # repeated cuts (e.g. NEXRAD SAILS) add nothing
            continue
        seen.add(round(el, 1))
        qvp = quasi_vertical_profile(ds)
        ml = detect_melting_layer(qvp, **kw)
        first = first or (ml, qvp)
        if ml.detected:
            return ml, qvp
    return first


class MeltingLayerTracker:
    """Keeps a melting-layer estimate through a sequence of volumes.

    Detections that jump by more than ``max_jump_km`` from the last accepted one are rejected (the melting
    layer rarely moves that fast; one KTLX volume gave a spurious layer 2 km low). If ``prior_top_km`` is
    given (e.g. a model 0 C height above the radar), detections whose top is more than ``max_prior_km`` away
    from it are also rejected. Without a new accepted detection the last one is kept for ``max_age_min``.
    These thresholds are radar-wow assumptions, not from the paper.
    """

    def __init__(self, max_jump_km: float = 0.5, max_age_min: float = 180.0, max_prior_km: float = 1.0):
        self.max_jump_km, self.max_age_min, self.max_prior_km = max_jump_km, max_age_min, max_prior_km
        self.current, self.time = None, None

    def update(self, time, ml: MeltingLayer, prior_top_km=None):
        import pandas as pd

        t = pd.Timestamp(time)
        status = "kept"
        if ml.detected:
            if prior_top_km is not None and abs(ml.top_km - prior_top_km) > self.max_prior_km:
                status = "rejected (far from prior)"
            elif self.current is not None and abs(ml.top_km - self.current.top_km) > self.max_jump_km:
                status = "rejected (jump)"
            else:
                self.current, self.time, status = ml, t, "accepted"
        if self.current is not None and (t - self.time).total_seconds() / 60 > self.max_age_min:
            self.current, status = None, "expired"
        return self.current, status
