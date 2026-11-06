import numpy as np
import xarray as xr

from radarwow.classification import MeltingLayer, MeltingLayerTracker, detect_melting_layer


def _qvp(ml_bottom=2.4, ml_top=3.1, with_ml=True):
    h = np.arange(0.1, 6.0, 0.05)
    z = np.where(h < ml_bottom, 32.0, np.where(h < ml_top, 32.0, 32.0 - 6 * (h - ml_top)))
    rho = np.full(h.shape, 0.99)
    zdr = np.full(h.shape, 0.6)
    if with_ml:
        core = (h > ml_bottom) & (h < ml_top)
        bump = np.sin(np.pi * (h - ml_bottom) / (ml_top - ml_bottom))
        z = np.where(core, z + 8 * bump, z)
        rho = np.where(core, 0.99 - 0.08 * bump, rho)
        zdr = np.where(core, 0.6 + 1.2 * bump, zdr)
    return xr.Dataset({"DBZH": ("height", z), "RHOHV": ("height", rho), "ZDR": ("height", zdr)},
                      coords={"height": h}, attrs={"elevation_deg": 19.5})


def test_detects_bright_band_boundaries():
    ml = detect_melting_layer(_qvp())
    assert ml.detected
    assert abs(ml.bottom_km - 2.4) < 0.25 and abs(ml.top_km - 3.1) < 0.25 and ml.bottom_km < ml.peak_km < ml.top_km


def test_no_signature_no_detection():
    assert not detect_melting_layer(_qvp(with_ml=False)).detected


def test_tracker_rejects_jumps_and_keeps_last():
    tr = MeltingLayerTracker(max_jump_km=0.5)
    good = MeltingLayer(True, bottom_km=2.5, top_km=3.2)
    bad = MeltingLayer(True, bottom_km=1.0, top_km=1.3)
    assert tr.update("2022-05-05 07:00", good)[1] == "accepted"
    cur, status = tr.update("2022-05-05 07:10", bad)
    assert status == "rejected (jump)" and cur.top_km == 3.2
    assert tr.update("2022-05-05 12:00", MeltingLayer(False))[1] == "expired"
