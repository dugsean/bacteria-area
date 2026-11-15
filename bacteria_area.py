"""Area of bacteria (fluorescent biomass) per image.

For every image the segmented mask is measured and reported as

    A_mask = number of mask pixels x (um/px)^2      [um^2]

and, if the pixel size is known, also as the covered fraction of the field.

Pipeline (identical for every image, all parameters are constants below):

  1. 3x3 median filter (hot-pixel removal)
  2. rolling-ball background subtraction (r = 40 px)
  3. small saturated artefacts removed (compact clusters <= 80 px above
     80 % of the image maximum on a dim surround -> local median)
  4. threshold = max(log-Otsu, 0.15 x image maximum)         [detection]
  5. keep pixels brighter than half of the LOCAL maximum (21 x 21 px) [refinement]
  6. discard objects < 20 px, close the mask (r = 2 px), fill holes

The script also writes an overlay PNG per image showing exactly which pixels
are counted.

Usage
    # one image
    python bacteria_area.py --input field.tif --outdir results/area

    # a directory of images (searched recursively)
    python bacteria_area.py --input images/ --outdir results/area

    # a Bio-Formats export (manifest.tsv + raw/ planes) - calibration is
    # taken from the manifest
    python bacteria_area.py --input ../analysis/5A1_DPI1 --outdir results/area

    # explicit calibration
    python bacteria_area.py --input images/ --um-per-px 0.871175

    # installed as a console script
    bacteria-area --input images/ --outdir results/area
"""

import argparse
import inspect
import json
import os

import numpy as np
import pandas as pd
import scipy.ndimage as ndi
from skimage.filters import threshold_otsu
from skimage.measure import label, regionprops
from skimage.morphology import closing, disk, remove_small_objects
from skimage.restoration import rolling_ball

__version__ = "0.1.0"

# ---------------------------------------------------------------- parameters
PARAMS = dict(
    channel=0,          # channel of a multi-channel image that carries the signal
    denoise=3,          # median filter size (px)
    ball_radius=40,     # rolling-ball background radius (px); 0 disables
    cand_frac=0.15,     # detection threshold as a fraction of the image maximum
    half_max=0.50,      # refinement threshold as a fraction of the LOCAL maximum
    local_win=21,       # local-maximum window (px)
    min_size=20,        # discard connected objects smaller than this (px)
    close_radius=2,     # morphological closing radius (px)
    fill_holes=True,
    remove_hot=True,    # remove small saturated hot pixels before thresholding
    hot_frac=0.80,      # pixels above this fraction of the max count as hot
    hot_max_area=80,    # hot clusters smaller than this (px) are artefacts
    hot_max_dim=12,     # ... and with a bounding box no larger than this (px)
)

IMAGE_EXT = (".tif", ".tiff", ".png", ".jpg", ".jpeg", ".raw")


# scikit-image renamed the argument of remove_small_objects in 0.26
# (min_size -> max_size, and the comparison became "<="); support both.
_ROSO_HAS_MAX_SIZE = "max_size" in inspect.signature(remove_small_objects).parameters


def remove_small_objects_compat(mask, min_size):
    """remove_small_objects with the same meaning on old and new scikit-image."""
    if _ROSO_HAS_MAX_SIZE:
        return remove_small_objects(mask, max_size=min_size)
    return remove_small_objects(mask, min_size=min_size)


# ------------------------------------------------------------------- helpers
def remove_hot_pixels(x, hot_frac=0.80, hot_max_area=80, hot_max_dim=12,
                      inpaint_radius=4):
    """Replace small saturated specks (camera hot pixels) by the local median.

    A hot artefact is a *small compact* cluster of near-saturated pixels on a
    much dimmer surround; saturated filament cores (long strips) and compact
    bright objects (equally bright surround) are kept.  Without this step the
    image maximum - and with it every threshold derived from the image
    maximum - would be set by the artefact rather than by the sample.
    """
    if x.max() <= 0:
        return x
    level = hot_frac * float(x.max())
    if level < 1:
        return x
    hot = x >= level
    lab = label(hot, connectivity=2)
    out = x.copy()
    hotmask = np.zeros_like(out, dtype=bool)
    for rp in regionprops(lab):
        h_, w_ = rp.bbox[2] - rp.bbox[0], rp.bbox[3] - rp.bbox[1]
        if rp.area > hot_max_area or max(h_, w_) > hot_max_dim:
            continue
        ys, xs = rp.coords[:, 0], rp.coords[:, 1]
        ring = out[max(0, rp.bbox[0] - 4):rp.bbox[2] + 4,
                   max(0, rp.bbox[1] - 4):rp.bbox[3] + 4]
        if np.median(ring) > 0.5 * level:
            continue            # surround is bright -> a genuine small object
        hotmask[ys, xs] = True
    if not hotmask.any():
        return out
    h, w = out.shape
    yy, xx = np.mgrid[0:h, 0:w]
    for y, x_ in zip(*np.nonzero(hotmask)):
        sel = ((yy - y) ** 2 + (xx - x_) ** 2 <= inpaint_radius ** 2) & (~hotmask)
        vals = out[sel]
        if vals.size:
            out[y, x_] = np.median(vals)
    return out


def preprocess(img, params):
    """Median filter -> rolling-ball background subtraction -> hot pixels."""
    x = img.astype(np.float32)
    if params["denoise"] and params["denoise"] > 1:
        x = ndi.median_filter(x, size=params["denoise"])
    if params["ball_radius"] and params["ball_radius"] > 0:
        bg = rolling_ball(x, radius=params["ball_radius"])
        x = np.clip(x - bg, 0, None)
    if params["remove_hot"]:
        x = remove_hot_pixels(x, params["hot_frac"], params["hot_max_area"],
                              params["hot_max_dim"])
    return x


def segment(x, params):
    """Two-stage segmentation -> (mask, detection threshold).

    Stage 1 (detection): a generous, image-relative threshold (the larger of
    a log-transformed Otsu value and a fixed fraction of the image maximum)
    finds every candidate filament, including dim ones that a global Otsu
    threshold misses.

    Stage 2 (local half-maximum refinement): a pixel is kept only if it is
    brighter than half of the *local* maximum intensity.  The reference is
    local, so saturated filaments lose their out-of-focus halo while
    filaments that fade along their length keep their dim end; the result is
    independent of the overall brightness of the sample.
    """
    t_log = float(threshold_otsu(np.log1p(x))) if x.max() > x.min() else x.max() + 1
    t = max(np.expm1(t_log), params["cand_frac"] * float(x.max()))
    t = max(t, 1.0)

    local_max = ndi.maximum_filter(x, size=params["local_win"])
    mask = (x > t) & (x > params["half_max"] * local_max)
    if params["min_size"] > 0:
        mask = remove_small_objects_compat(mask, params["min_size"])
    if params["close_radius"] > 0:
        mask = closing(mask, disk(params["close_radius"]))
    if params["fill_holes"]:
        mask = ndi.binary_fill_holes(mask)
    return mask, t


# ------------------------------------------------------------------- reading
def read_manifest(path):
    """Read the manifest.tsv written by export_lif.py."""
    lines = open(path).read().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, l.split("\t"))) for l in lines[1:] if l.strip()]


def _um_per_px_from_tif(page):
    """um/px from the TIFF resolution tags, or None."""
    try:
        xres = page.tags["XResolution"].value      # (numerator, denominator)
        yres = page.tags["YResolution"].value
        unit = int(page.tags["ResolutionUnit"].value)
    except Exception:
        return None
    px_per_cm = float(xres[0]) / float(xres[1] or 1)
    py_per_cm = float(yres[0]) / float(yres[1] or 1)
    if unit == 3:      # centimetre
        return 1e4 / px_per_cm, 1e4 / py_per_cm
    if unit == 2:      # inch
        return 25400.0 / px_per_cm, 25400.0 / py_per_cm
    return None


def read_image(path, channel=0):
    """Read one image -> (2-D array of the chosen channel, um/px or None).

    Supports TIFF (tifffile), PNG/JPEG (Pillow) and the header-less raw
    planes written by export_lif.py (only with a manifest).
    """
    ext = os.path.splitext(path)[1].lower()
    if ext in (".tif", ".tiff"):
        import tifffile
        with tifffile.TiffFile(path) as tf:
            arr = tf.asarray()
            cal = None
            try:
                cal = _um_per_px_from_tif(tf.pages[0])
            except Exception:
                cal = None
    elif ext in (".png", ".jpg", ".jpeg"):
        from PIL import Image
        arr = np.asarray(Image.open(path))
        cal = None
    elif ext == ".raw":
        raise ValueError("raw planes need a manifest.tsv (width/height)")
    else:
        raise ValueError("unsupported image type: %s" % ext)

    arr = np.squeeze(arr)
    if arr.ndim == 3:
        # (h, w, c) or (c, h, w) or (h, w, z)
        if arr.shape[0] <= 4 and arr.shape[2] > 4:
            arr = np.moveaxis(arr, 0, -1)      # (c, h, w) -> (h, w, c)
        c = int(channel) % arr.shape[2]
        arr = arr[..., c]
    elif arr.ndim != 2:
        raise ValueError("unsupported image shape %s in %s" % (arr.shape, path))
    return arr, cal


# ------------------------------------------------------------------ overlays
def mask_border(mask):
    """1-px border of a binary mask (the pixels at the edge of the area)."""
    if not mask.any():
        return np.zeros_like(mask)
    return mask & ~ndi.binary_erosion(mask, structure=np.ones((3, 3), bool),
                                      border_value=0)


def draw_overlay(img, mask, um, path, title, area_um2, area_px, dpi=150):
    """Render one overlay: image + counted area (fill + border) + scale bar."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    base = np.stack([img] * 3, axis=-1).astype(np.float32)
    hi = float(np.percentile(base, 99.9))
    if hi > 0:
        base = np.clip(base / hi, 0, 1)
    h, w = img.shape

    fig, ax = plt.subplots(figsize=(10, 6.7))
    ax.imshow(base, cmap="gray", vmin=0, vmax=1)
    ax.imshow(np.where(mask, 1.0, np.nan), cmap=ListedColormap(["#00ccff"]),
              vmin=0, vmax=1, alpha=0.35, interpolation="nearest")
    ax.imshow(np.where(mask_border(mask), 1.0, np.nan),
              cmap=ListedColormap(["#ffffff"]),
              vmin=0, vmax=1, alpha=1.0, interpolation="nearest")

    if um:
        bar_um = 100
        bar_px = bar_um / um
        x0, y0 = w - bar_px - 40, h - 40
        ax.plot([x0, x0 + bar_px], [y0, y0], color="white", lw=4)
        ax.text(x0 + bar_px / 2, y0 - 12, "%d µm" % bar_um, color="white",
                ha="center", fontsize=10)

    if um:
        txt = ("%s\narea (mask) = %.0f µm²  = %.0f px  = %.2f %% of the field"
               % (title, area_um2, area_px, 100.0 * area_px / mask.size))
    else:
        txt = ("%s\narea (mask) = %.0f px = %.2f %% of the field\n"
               "(pixel size unknown - area in px only)"
               % (title, area_px, 100.0 * area_px / mask.size))
    ax.text(0.012, 0.985, txt, transform=ax.transAxes, color="white",
            fontsize=9, va="top", ha="left",
            bbox=dict(facecolor="black", alpha=0.45, edgecolor="none", pad=4))
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


# -------------------------------------------------------------- image lookup
def find_images(root, recursive=True):
    """All image files below root, sorted."""
    out = []
    if recursive:
        for dirpath, _, names in os.walk(root):
            for n in sorted(names):
                if n.lower().endswith(IMAGE_EXT):
                    out.append(os.path.join(dirpath, n))
    else:
        for n in sorted(os.listdir(root)):
            p = os.path.join(root, n)
            if os.path.isfile(p) and n.lower().endswith(IMAGE_EXT):
                out.append(p)
    return sorted(out)


def group_of(path, root):
    """Group label: the directory of the file relative to the input root."""
    if os.path.isfile(root):
        return os.path.splitext(os.path.basename(path))[0]
    rel = os.path.relpath(os.path.dirname(path), root)
    if rel == ".":
        return os.path.basename(os.path.normpath(root))
    return rel.split(os.sep)[0]


def make_figure(df, outdir):
    """Per-image area, one panel per group."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    mpl.rcParams.update({
        "figure.dpi": 300, "savefig.dpi": 300, "font.size": 8,
        "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
    })
    groups = list(df["group"].unique())
    if not groups:
        return
    ncol = min(3, len(groups))
    nrow = int(np.ceil(len(groups) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 2.8 * nrow),
                             squeeze=False)
    for ax, g in zip(axes.ravel(), groups):
        sub = df[df["group"] == g].reset_index(drop=True)
        cal = sub[sub.um_per_px.notna()]
        y = cal.mask_area_um2 if len(cal) else sub.mask_area_px
        unit = "µm²" if len(cal) else "px"
        ax.bar(np.arange(len(sub)), y, color="#0072B2", width=0.7)
        ax.set_ylabel("mask area (%s)" % unit)
        ax.set_title(g, fontsize=8)
        ax.set_xticks([])
    for ax in axes.ravel()[len(groups):]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "figure_area.png"))
    plt.close(fig)


# ---------------------------------------------------------------------- main
def process_image(path, params, um=None):
    """Segment one image -> (mask, measurement dict)."""
    img, cal = read_image(path, params["channel"])
    if um is None and cal is not None:
        um = (cal[0] + cal[1]) / 2.0
    x = preprocess(img, params)
    mask, t = segment(x, params)
    area_px = float(mask.sum())
    full = float(np.prod(img.shape))
    lab = label(mask, connectivity=2)
    sat = 100.0 * float((img >= np.iinfo(img.dtype).max).sum()) / img.size
    m = dict(
        width=int(img.shape[1]), height=int(img.shape[0]),
        um_per_px=um,
        mask_area_px=area_px,
        mask_area_um2=area_px * um ** 2 if um else np.nan,
        coverage_pct=100.0 * area_px / full,
        n_components=int(lab.max()) if lab.max() > 0 else 0,
        threshold=float(t),
        sat_pct=sat,
    )
    return img, mask, m


def main():
    ap = argparse.ArgumentParser(
        description="Measure the area of bacteria in images (one image or a "
                    "directory) and write per-image overlays.")
    ap.add_argument("--version", action="version",
                    version="bacteria-area %s" % __version__)
    ap.add_argument("--input", action="append", required=True,
                    help="image file or directory (repeatable)")
    ap.add_argument("--outdir", default="results/area")
    ap.add_argument("--um-per-px", type=float, default=None,
                    help="pixel size (um/px); overrides the file metadata")
    ap.add_argument("--manifest", default=None,
                    help="manifest.tsv with the pixel size "
                         "(auto-detected when the input directory holds one)")
    ap.add_argument("--channel", type=int, default=PARAMS["channel"])
    ap.add_argument("--no-recursive", action="store_true",
                    help="do not descend into sub-directories")
    ap.add_argument("--no-overlays", action="store_true")
    ap.add_argument("--overlay-dpi", type=int, default=150)
    ap.add_argument("--min-size", type=int, default=None, dest="min_size")
    ap.add_argument("--ball-radius", type=int, default=None, dest="ball_radius")
    ap.add_argument("--cand-frac", type=float, default=None, dest="cand_frac")
    ap.add_argument("--half-max", type=float, default=None, dest="half_max")
    ap.add_argument("--local-win", type=int, default=None, dest="local_win")
    ap.add_argument("--close-radius", type=int, default=None, dest="close_radius")
    args = ap.parse_args()

    params = dict(PARAMS)
    params["channel"] = args.channel
    for k in ("min_size", "ball_radius", "cand_frac", "half_max", "local_win",
              "close_radius"):
        v = getattr(args, k)
        if v is not None:
            params[k] = v

    outdir = os.path.abspath(args.outdir)
    ovdir = os.path.join(outdir, "overlays")
    os.makedirs(outdir, exist_ok=True)
    if not args.no_overlays:
        os.makedirs(ovdir, exist_ok=True)

    rows = []
    used_stems = set()
    for inp in args.input:
        inp = os.path.abspath(inp)
        manifest = args.manifest
        jobs = []                       # (path, group, um/px or None, rec)

        if os.path.isfile(inp):
            base = os.path.basename(inp)
            rec = None
            if base.lower().endswith(".raw"):
                # a raw plane needs the manifest of its export
                man = manifest
                if man is None:
                    for d in (os.path.dirname(inp), os.path.dirname(os.path.dirname(inp))):
                        cand = os.path.join(d, "manifest.tsv")
                        if os.path.isfile(cand):
                            man = cand
                            break
                if man and os.path.isfile(man):
                    for r in read_manifest(man):
                        if r["file"] == base:
                            rec = r
                            break
                if rec is None:
                    raise SystemExit("%s is a raw plane but has no manifest "
                                     "entry (use --manifest)" % inp)
            jobs.append((inp, os.path.splitext(base)[0], None, rec))
        else:
            man = manifest or os.path.join(inp, "manifest.tsv")
            if manifest is None and os.path.isfile(man):
                manifest = man
            recs = read_manifest(manifest) if manifest and os.path.isfile(manifest) else []
            sample = recs[0].get("sample", os.path.basename(inp)) if recs \
                else os.path.basename(inp)
            handled = set()
            if recs:
                # Bio-Formats export: raw planes described by the manifest
                raw_dir = os.path.join(inp, "raw")
                if not os.path.isdir(raw_dir):
                    raw_dir = inp
                for rec in recs:
                    p = os.path.join(raw_dir, rec["file"])
                    if os.path.isfile(p):
                        jobs.append((p, rec.get("sample", sample), None, rec))
                        handled.add(os.path.abspath(p))
            # any other image file below the directory (tif/png/jpg/...)
            for p in find_images(inp, recursive=not args.no_recursive):
                if os.path.abspath(p) in handled:
                    continue
                jobs.append((p, group_of(p, inp), None, None))

        if not jobs:
            print("no images found in", inp)
            continue

        for path, group, um, rec in jobs:
            base = os.path.basename(path)
            stem = os.path.splitext(base)[0]
            ext = os.path.splitext(base)[1].lstrip(".").lower()
            if stem in used_stems:          # e.g. field.tif and field.png
                stem = "%s_%s" % (stem, ext)
            used_stems.add(stem)
            um = args.um_per_px if args.um_per_px is not None else um
            try:
                if rec is not None:
                    img, mask, m = process_raw(path, rec, params, um)
                else:
                    img, mask, m = process_image(path, params, um)
            except Exception as e:
                print("FAILED %s: %s" % (path, e))
                continue
            m.update(group=group, image=stem, path=path)
            rows.append(m)
            if not args.no_overlays:
                draw_overlay(img, mask, m["um_per_px"],
                             os.path.join(ovdir, stem + "_area.png"),
                             "%s  %s" % (m["group"], stem),
                             m["mask_area_um2"], m["mask_area_px"],
                             dpi=args.overlay_dpi)
            if m["um_per_px"]:
                print("%-28s area = %12.1f um^2  (%7.3f %% of field)"
                      % (stem, m["mask_area_um2"], m["coverage_pct"]))
            else:
                print("%-28s area = %12.0f px    (%7.3f %% of field)"
                      "   [no pixel size - area in px only]"
                      % (stem, m["mask_area_px"], m["coverage_pct"]))

    if not rows:
        raise SystemExit("no images processed")

    df = pd.DataFrame(rows)
    df["params"] = json.dumps(params, sort_keys=True)
    cols = ["group", "image", "path", "width", "height", "um_per_px",
            "mask_area_px", "mask_area_um2", "coverage_pct", "n_components",
            "threshold", "sat_pct", "params"]
    df = df[cols]
    df.to_csv(os.path.join(outdir, "areas.csv"), index=False)

    # ---- per-group summary (statistics over the calibrated images only;
    #      images without a pixel size contribute to the pixel-based columns)
    summ = []
    for g, sub in df.groupby("group"):
        cal = sub[sub.um_per_px.notna()]
        v = cal.mask_area_um2
        summ.append(dict(
            group=g, n_images=len(sub), n_calibrated=len(cal),
            mean_area_um2=v.mean() if len(v) else np.nan,
            sd_area_um2=v.std(ddof=1) if len(v) > 1 else np.nan,
            median_area_um2=v.median() if len(v) else np.nan,
            q25_area_um2=v.quantile(0.25) if len(v) else np.nan,
            q75_area_um2=v.quantile(0.75) if len(v) else np.nan,
            mean_mask_area_px=sub.mask_area_px.mean(),
            sd_mask_area_px=sub.mask_area_px.std(ddof=1) if len(sub) > 1 else np.nan,
            mean_coverage_pct=sub.coverage_pct.mean(),
            sd_coverage_pct=sub.coverage_pct.std(ddof=1) if len(sub) > 1 else np.nan,
            mean_n_components=sub.n_components.mean(),
            mean_sat_pct=sub.sat_pct.mean(),
        ))
    summ = pd.DataFrame(summ)
    summ.to_csv(os.path.join(outdir, "area_summary.csv"), index=False)
    with open(os.path.join(outdir, "area_summary.md"), "w") as fh:
        fh.write(summ.to_csv(index=False))
    with open(os.path.join(outdir, "parameters.json"), "w") as fh:
        json.dump(params, fh, indent=2, sort_keys=True)
    try:
        make_figure(df, outdir)
    except Exception as e:
        print("figure skipped:", e)

    print("\n" + summ.to_string(index=False))
    print("\nwrote", os.path.join(outdir, "areas.csv"))
    if not args.no_overlays:
        print("overlays ->", ovdir)


def process_raw(path, rec, params, um=None):
    """Raw plane written by export_lif.py (width/height from the manifest)."""
    w, h = int(rec["width"]), int(rec["height"])
    nch = int(rec["bytes"]) // (w * h)
    data = np.fromfile(path, dtype=np.uint8)
    img = data.reshape(h, w, nch)[..., params["channel"]]
    if um is None:
        um = (float(rec["um_per_px_x"]) + float(rec["um_per_px_y"])) / 2.0
    x = preprocess(img, params)
    mask, t = segment(x, params)
    area_px = float(mask.sum())
    lab = label(mask, connectivity=2)
    m = dict(
        width=w, height=h, um_per_px=um,
        mask_area_px=area_px,
        mask_area_um2=area_px * um ** 2,
        coverage_pct=100.0 * area_px / float(w * h),
        n_components=int(lab.max()) if lab.max() > 0 else 0,
        threshold=float(t),
        sat_pct=100.0 * float((img >= 255).sum()) / img.size,
    )
    return img, mask, m


if __name__ == "__main__":
    main()
