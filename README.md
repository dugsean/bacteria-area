# bacteria-area

Measure the **area covered by bacteria** in fluorescence micrographs — one
number per image, plus an overlay that shows exactly which pixels were
counted.

```
A_mask = mask pixels × (µm/px)²      [µm²]
```

No filament tracing, no length estimate, no skeleton: the quantity reported
is the projected area of the fluorescent biomass.

The segmentation was tuned on filamentous cyanobacteria (Leica SP8 wide-field
fluorescence, 8-bit, 0.87 µm/px), but the script is generic: it segments
bright objects on a dark background in any TIFF/PNG/JPEG image.

---

## Install (Linux + macOS, venv — no conda)

```bash
git clone https://github.com/dugsean/bacteria-area.git
cd bacteria-area
./install.sh            # add --dev to also install pytest
```

The script creates `./.venv`, installs the dependencies and installs the
package in editable mode. Nothing is installed system-wide, no root needed.

```bash
source .venv/bin/activate      # optional
bacteria-area --help
```

Requirements: Python ≥ 3.9 and a C compiler is *not* needed (binary wheels
exist for numpy/scipy/scikit-image on macOS and Linux).

If the virtual environment cannot be created on Debian/Ubuntu:
`sudo apt install python3-venv` (or `python3.13-venv`).

---

## Use

```bash
# one image
bacteria-area --input field.tif --outdir results

# a directory of images (searched recursively, grouped by sub-directory)
bacteria-area --input images/ --outdir results

# explicit pixel size (overrides the file metadata)
bacteria-area --input images/ --um-per-px 0.871175 --outdir results
```

**Pixel size** is taken from (in this order): `--um-per-px`, the TIFF
resolution tags, or `manifest.tsv` of a Bio-Formats export. If none is
available the area is reported in **pixels only** and the table says so —
never guess a calibration.

**Input formats**: TIFF (`.tif/.tiff`, any bit depth, multi-channel — pick the
channel with `--channel`), PNG, JPEG. Multi-channel images are reduced to the
chosen channel (default 0).

### Try it without data

```bash
python examples/make_synthetic.py
bacteria-area --input examples/data --outdir results/example
```

### Output

```
results/
  areas.csv                 one row per image
  area_summary.csv/.md      per-group statistics
  figure_area.png           per-image area, one panel per group
  parameters.json           the parameters used
  overlays/<image>_area.png per-image overlay of the counted area
```

`areas.csv` columns: `group`, `image`, `path`, `width`, `height`, `um_per_px`,
`mask_area_px`, `mask_area_um2`, `coverage_pct` (share of the field covered),
`n_components` (connected structures in the mask), `threshold`, `sat_pct`
(% clipped pixels), `params` (the full parameter set, JSON).

The overlay shows the raw image (contrast-stretched), the counted area as a
translucent cyan fill with its border in white, the area in µm² and a 100 µm
scale bar.

---

## Method

1. **Pre-processing** — 3 × 3 median filter (hot-pixel removal); rolling-ball
   background subtraction (radius 40 px); removal of small saturated
   artefacts (compact clusters ≤ 80 px above 80 % of the image maximum on a
   dim surround, replaced by the local median) — otherwise a single artefact
   sets the image maximum and hence every threshold derived from it.
2. **Segmentation, two stages**
   * *detection*: threshold = max(log-Otsu, 0.15 × image maximum) — generous
     and image-relative, so dim objects are not missed;
   * *refinement*: a pixel is kept only if it is brighter than **half of the
     local maximum** (21 × 21 px window), which removes the out-of-focus halo
     around bright objects while keeping the dim end of a fading one;
   * objects < 20 px removed, mask closed (radius 2 px), holes filled.
3. **Area** — `A_mask` = mask pixels × (µm/px)²; `coverage_pct` =
   100 × mask pixels / image pixels.

All parameters are constants at the top of `bacteria_area.py`, can be
overridden from the command line (`--min-size`, `--ball-radius`,
`--cand-frac`, `--half-max`, `--local-win`, `--close-radius`, `--channel`),
and are stored in every output row.

### Caveats

* **Saturated images inflate the mask.** A clipped intensity profile is
  flattened, which widens the segmented region: for such fields the measured
  area is an **upper** bound. `sat_pct` flags them.
* **Projection.** The area is the area *in the plane*: filaments stacked along
  the optical axis occupy the same pixels and are not counted separately.
* **Detection limit.** Objects < 20 px (≈ 15 µm² at 0.87 µm/px) are discarded.
* **Parameters** were checked against alternatives: totals change by ≤ 5 %
  for background radii of 60–100 px and ≤ 2 % for a half-maximum factor of
  0.4 instead of 0.5.

---

## Leica `.lif` files (optional, needs Fiji)

Importing `.lif` needs [Fiji](https://fiji.sc/) with Bio-Formats. `lif_export.py`
runs inside Fiji and writes the planes of one `.lif` to raw files plus a
`manifest.tsv` (which carries the pixel size):

```bash
IMG_PATH=/path/to/file.lif OUT_DIR=analysis/file \
  fiji --headless --run lif_export.py
bacteria-area --input analysis/file --outdir results
```

`bacteria-area` recognises a directory that contains `manifest.tsv` and reads
the raw planes it lists.

---

## Tests

```bash
.venv/bin/python -m pytest tests -q
```

The tests check that a blank image yields no area, that synthetic bars are
segmented with an area close to the truth, that the µm² conversion and the
overlay writer work, that the CLI produces a table, and that the pixel size is
read from the TIFF tags.

---

## Repository layout

```
bacteria_area.py      the measurement (single file, no hidden state)
lif_export.py         optional Fiji/Bio-Formats step for .lif input
install.sh            venv installer (Linux + macOS)
pyproject.toml        package metadata, dependencies, console script
requirements.txt      the same dependencies for pip -r
examples/             synthetic example data
tests/                pytest smoke tests
docs/METHODS.md       manuscript methods paragraph
```

## License

MIT — see `LICENSE`.
