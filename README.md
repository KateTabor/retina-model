# retina-model

Simple, reproducible pipeline for:

1. **Generating an image set** from a Jupyter notebook + YAML config.
2. **Synthesizing retinal responses** with ISETBio (human optics; midget RGC mosaics with/without S‑cone input to the surround).
3. **Exporting per‑image mRGC responses** for downstream ML/analysis.

> **Note on large files**: images and `.mat` outputs are tracked with **Git LFS**.

---

## 1) Image set generation (Python/Jupyter)

* **Inputs**: a YAML config (paths, categories/splits, sizes, random seed, etc.).
* **Entry point**: a Jupyter notebook that reads the YAML and writes the dataset.
* **Outputs** (per image, written into split folders):

  * Raster images: `.tif`, `.png`, `.jpeg`
  * Arrays/metadata for MATLAB workflows: `.mat`, `.csv`

**Typical layout** (example):

```
imageset/
  train/
  val/
  test/
```

**Key Python deps** (typical): `colour-science`, `numpy`, `pillow`, `matplotlib`, `pyyaml`, `jupyter`.

> The generator produces both raster images (for quick inspection/ML) **and** `.mat` files (float arrays / supporting metadata used by MATLAB + ISETBio steps).

---

## 2) ISETBio model (MATLAB)

We synthesize retinal responses to each dataset image using ISETBio/ISETCam.

* **Optics**: human optics model.
* **Neural front‑end**: on midget RGC mosaic (mRGC).
* **Surround option**: run variants with **S‑cone participation in the surround on** (two comparable mosaics / parameterizations).

**Per image**, we run the mRGC mosaic over the stimulus and save responses.

---

## 3) Output: mRGC responses (MATLAB `.mat`)

Responses are saved in dated run folders (example shown):

```
20251108-mRGCResp/
  train/
  val/
  test/
```

Each file contains the response for a single source image (naming mirrors the image filename). 

---

## Repository layout (representative)

```
notebooks/                   # Jupyter for image generation
configs/                     # YAML for dataset specs
imageset/                    # train/val/test raster + .mat from the generator
2025MMDD-mRGCResp/           # train/val/test mRGC outputs (.mat)
matlMatlab_ISETBio_scripts/  # MATLAB/ISETBio scripts (building optics/mosaics, batch run)
```
---

## How to run (quick sketch)

**Image set (Python):**

1. Create/adjust `configs/imageset.yaml`.
2. Open the generator notebook and run all cells.
3. Verify counts in `imageset/{train,val,test}` and spot‑check images.

**Retina synthesis (MATLAB/ISETBio):**

1. Ensure ISETBio/ISETCam on path; initialize.
2. Point the script to `imageset/` and choose the mRGC model (w/o S‑cone‑surround).
3. Run the batch; outputs go to a dated `*‑mRGCResp/` folder with `{train,val,test}`.

---

## Git LFS

Large binary outputs are tracked via LFS (patterns include `*.tif`, `*.png`, `*.jpeg`, `*.mat`).