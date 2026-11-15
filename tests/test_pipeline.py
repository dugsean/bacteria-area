"""Smoke tests for the area measurement.

Run with:  python -m pytest tests -q
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import bacteria_area as ba  # noqa: E402


UM_PER_PX = 0.5


def synthetic_field(shape=(512, 512), bars=((60, 60, 300, 12), (60, 200, 300, 12)),
                    seed=0):
    """A dark field with a few bright bars of known geometry (px)."""
    img = np.zeros(shape, dtype=np.float32)
    for y0, x0, length, width in bars:
        img[y0:y0 + width, x0:x0 + length] = 200.0
    rng = np.random.default_rng(seed)
    img += rng.normal(0, 3, img.shape)
    return np.clip(img, 0, 255).astype(np.uint8)


def test_blank_image_gives_no_area():
    img = np.full((256, 256), 20, dtype=np.uint8)
    mask, t = ba.segment(ba.preprocess(img, ba.PARAMS), ba.PARAMS)
    assert mask.sum() == 0


def test_bars_are_segmented_and_area_is_close_to_truth():
    img = synthetic_field()
    params = dict(ba.PARAMS)
    mask, _ = ba.segment(ba.preprocess(img, params), params)
    measured_px = float(mask.sum())
    truth_px = 2 * 300 * 12
    # the half-maximum rule shrinks the bar slightly, the closing adds a little
    assert 0.6 * truth_px < measured_px < 1.4 * truth_px


def test_area_conversion_and_overlay(tmp_path):
    img = synthetic_field()
    params = dict(ba.PARAMS)
    x = ba.preprocess(img, params)
    mask, _ = ba.segment(x, params)
    area_px = float(mask.sum())
    area_um2 = area_px * UM_PER_PX ** 2
    assert area_um2 == pytest.approx(area_px * 0.25)

    png = tmp_path / "overlay.png"
    ba.draw_overlay(img, mask, UM_PER_PX, str(png), "test", area_um2, area_px)
    assert png.exists() and png.stat().st_size > 0


def test_cli_runs_on_a_directory(tmp_path):
    import tifffile

    d = tmp_path / "imgs"
    d.mkdir()
    for i in range(2):
        tifffile.imwrite(str(d / ("field%02d.tif" % i)), synthetic_field(seed=i),
                         resolution=(1e4 / UM_PER_PX, 1e4 / UM_PER_PX),
                         resolutionunit="CENTIMETER")
    out = tmp_path / "out"
    argv = sys.argv
    sys.argv = ["bacteria-area", "--input", str(d), "--outdir", str(out),
                 "--no-overlays"]
    try:
        ba.main()
    finally:
        sys.argv = argv
    csv = out / "areas.csv"
    assert csv.exists()
    text = csv.read_text()
    assert "field00" in text and "field01" in text
    assert "mask_area_um2" in text


def test_calibration_is_read_from_tiff_tags(tmp_path):
    import tifffile

    p = tmp_path / "cal.tif"
    tifffile.imwrite(str(p), synthetic_field(),
                     resolution=(1e4 / UM_PER_PX, 1e4 / UM_PER_PX),
                     resolutionunit="CENTIMETER")
    img, cal = ba.read_image(str(p))
    assert img.shape == (512, 512)
    assert cal is not None
    assert cal[0] == pytest.approx(UM_PER_PX, rel=1e-3)
