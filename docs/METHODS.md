# Methods (manuscript paragraph — area of bacteria)

The area covered by filamentous cyanobacteria was measured in
wide-field fluorescence micrographs acquired with a Leica SP8 microscope and
stored as Leica `.lif` files. Images were imported with Bio-Formats 8.5.0
(via Fiji/ImageJ); the red detection channel (8-bit, 1536 × 1024 px) was used
for all measurements, the remaining channels containing only bleed-through of
the same fluorophore. Pixel size was taken from the file metadata
(0.8712 µm px⁻¹), giving a field of view of 1.19 mm². The analysed samples
were 1 and 10 days post inoculation (dpi) of one inoculation series, plus a
third sample acquired separately.

Raw images were denoised with a 3 × 3 median filter and the background was
flattened by subtracting a rolling-ball estimate (radius 40 px). Compact
saturated artefacts (≤ 80 px, ≥ 80 % of the image maximum, on a dim
surround) were replaced by the local median before segmentation, because
otherwise a single artefact determines the image maximum and hence any
threshold derived from it. Filaments were then segmented in two steps. First,
all candidate filaments were detected with a generous image-relative
threshold (the larger of a log-transformed Otsu threshold and 15 % of the
image maximum); a plain Otsu threshold detected only one of the three
filaments present in one field of the 1 dpi sample, whereas the rule used
found all three (61.0 µm versus 16.9 µm of segmented area). Second, each pixel
was kept only if its intensity exceeded half of the local maximum intensity
(21 × 21 px neighbourhood), which removes the out-of-focus halo around
saturated filaments while retaining the faint end of filaments that fade along
their length, and makes the segmentation independent of the overall brightness
of the sample. Objects smaller than 20 px were discarded, the binary mask was
closed (radius 2 px) and internal holes were filled. The same parameters were
applied to all images.

The area of bacteria in each field was computed as the number of segmented
mask pixels multiplied by the square of the pixel size,

    A = n_mask pixels × (0.8712 µm)²   [µm²],

and is reported together with the covered fraction of the field
(100 × n_mask / n_pixels) and the number of connected structures in the mask.
No tracing of individual filaments and no length estimate was performed; the
quantity reported is the projected area of the fluorescent biomass, so that
filaments crossing or lying on top of one another in the projection are
counted once.

The segmented area was 669.0 µm² per field (median; interquartile range
479.1–971.8 µm²; mean 751.3 ± 373.3 µm², n = 24 fields, i.e. 0.063 ± 0.031 %
of the field) at 1 dpi and 32 410.0 µm² per field (interquartile range
13 393.7–76 573.1 µm²; mean 51 687.7 ± 47 484.9 µm², n = 24 fields, i.e.
4.33 ± 3.98 % of the field) at 10 dpi. The fields within one file are
technical replicates of a single acquisition, so the biological replication is
one culture per time point and the field-to-field dispersion must not be read
as population variability.

The measurement assumes that all fluorescent structures are bacterial material.
Two limits apply. First, the area is a projection: filaments stacked along
the optical axis occupy the same pixels and are not counted separately. Second,
in the saturated sample (2.8–7.2 % clipped pixels) the flattened intensity
profile widens the segmented region, so the values obtained from it are upper
bounds; this sample should be re-acquired at a lower exposure before it is
compared quantitatively with the others. Objects smaller than 20 px
(≈ 15 µm²) fall below the detection limit. The segmentation parameters were
checked against alternatives: the measured area changes by ≤ 5 % for
background radii of 60–100 px and by ≤ 2 % when the local half-maximum factor
is changed from 0.5 to 0.4.

All analysis was performed in Python (3.9 or newer) with numpy, scipy,
scikit-image, pandas, matplotlib, tifffile and Pillow; the versions used for
the reported values were numpy 2.5.4, scipy 1.18.1, scikit-image 0.26.0,
pandas 3.1.1 and matplotlib 3.11.2. The analysis script, the per-image
measurement table and, for every field, an overlay of the exact pixels that
enter the calculation (the segmented mask drawn over the original image, with
its border marked) are provided as supplementary material. The code is
available at <https://github.com/dugsean/bacteria-area> (repository
`bacteria-area`).
