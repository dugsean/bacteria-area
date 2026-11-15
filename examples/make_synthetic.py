"""Write a synthetic example dataset so the pipeline can be tried out.

    python examples/make_synthetic.py

Creates examples/data/synthetic_field.tif (with a pixel size in the TIFF
tags), examples/data/synthetic_field.png and a second field in a
sub-directory, so that the per-group summary has something to show.
"""

import os

import numpy as np
import tifffile
from PIL import Image

UM_PER_PX = 0.5
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "data")


def field(seed=0, shape=(512, 512), bars=((80, 60, 320, 10),
                                          (200, 120, 260, 10),
                                          (320, 200, 200, 10))):
    """Dark field with bright bars (filament-like objects) of known geometry."""
    img = np.zeros(shape, dtype=np.float32)
    for y0, x0, length, width in bars:
        img[y0:y0 + width, x0:x0 + length] = 220.0
    rng = np.random.default_rng(seed)
    img += rng.normal(0, 3, img.shape)
    return np.clip(img, 0, 255).astype(np.uint8)


def main():
    os.makedirs(OUT, exist_ok=True)

    a = field(seed=0)
    tifffile.imwrite(os.path.join(OUT, "synthetic_field.tif"), a,
                     resolution=(1e4 / UM_PER_PX, 1e4 / UM_PER_PX),
                     resolutionunit="CENTIMETER")
    Image.fromarray(a).save(os.path.join(OUT, "synthetic_field.png"))

    b = field(seed=1, bars=((80, 60, 320, 10), (200, 120, 260, 10),
                             (320, 200, 200, 10), (150, 300, 180, 10)))
    sub = os.path.join(OUT, "second_group")
    os.makedirs(sub, exist_ok=True)
    tifffile.imwrite(os.path.join(sub, "synthetic_field_b.tif"), b,
                     resolution=(1e4 / UM_PER_PX, 1e4 / UM_PER_PX),
                     resolutionunit="CENTIMETER")

    truth_um2 = (320 * 10 + 260 * 10 + 200 * 10) * UM_PER_PX ** 2
    print("wrote", OUT)
    print("expected area of the first field: ~%.0f um^2 (bars only)" % truth_um2)


if __name__ == "__main__":
    main()
