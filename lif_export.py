"""Export every image of a Leica .lif to raw pixel files + a manifest TSV.

Step 1 of the pipeline; the measurement is done by bacteria_area.py, which
reads the raw planes of a directory that holds this manifest.

Env:
  IMG_PATH  input .lif
  OUT_DIR   output directory (raw/ + manifest.tsv)
"""
import os
import re

from loci.formats import ImageReader
from loci.formats import FormatTools
from loci.formats.ome import OMEXMLMetadataImpl

IMG_PATH = os.environ["IMG_PATH"]
OUT_DIR = os.environ["OUT_DIR"]
RAW_DIR = os.path.join(OUT_DIR, "raw")
if not os.path.isdir(RAW_DIR):
    os.makedirs(RAW_DIR)

base = re.sub(r"[^A-Za-z0-9._-]+", "_", os.path.splitext(os.path.basename(IMG_PATH))[0])

r = ImageReader()
store = OMEXMLMetadataImpl()
r.setMetadataStore(store)
r.setId(IMG_PATH)

umx = 0.8711752442996742
umy = 0.8711749755620724
try:
    px = store.getPixelsPhysicalSizeX(0)
    umx = float(px.value())
    umy = float(store.getPixelsPhysicalSizeY(0).value())
except:
    pass

cml = r.getCoreMetadataList()
print("core metadata entries:", len(cml), " um/px:", umx, umy)

manifest = open(os.path.join(OUT_DIR, "manifest.tsv"), "w")
manifest.write("file\tsample\tcore\tchannel\tz\tt\twidth\theight\tbytes\tpixel_type\tum_per_px_x\tum_per_px_y\n")

n_written = 0
for i in range(len(cml)):
    c = cml[i]
    r.setCoreIndex(i)
    n = r.getImageCount()
    print("core %d: %dx%d z=%d c=%d t=%d images=%d rgb=%s inter=%s" % (
        i, c.sizeX, c.sizeY, c.sizeZ, c.sizeC, c.sizeT, n, c.rgb, c.interleaved))
    for p in range(n):
        data = r.openBytes(p)
        w = r.getSizeX()
        h = r.getSizeY()
        cidx = p
        nch = 3 if c.rgb else r.getEffectiveSizeC()
        ch = cidx % nch
        z = (cidx // nch) % max(1, r.getSizeZ())
        t = cidx // (nch * max(1, r.getSizeZ()))
        ptype = FormatTools.getPixelTypeString(c.pixelType)
        raw = os.path.join(RAW_DIR, "%s_core%02d_c%02d_z%02d_t%02d.raw" % (base, i, ch, z, t))
        fh = open(raw, "wb")
        fh.write(data)
        fh.close()
        manifest.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%.6f\t%.6f\n" % (
            os.path.basename(raw), base, i, ch, z, t, w, h, len(data), ptype, umx, umy))
        n_written += 1

manifest.close()
r.close()
print("wrote %d raw planes" % n_written)
