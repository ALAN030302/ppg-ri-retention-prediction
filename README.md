# PPG-RI Retention Prediction Supporting Code, Data, and Models

This repository contains the public supporting package for the manuscript:

> Cross-Platform Retention Time Prediction of Environmental Contaminants: A Polypropylene Glycol Retention Index-Driven Graph Neural Network

The package provides processed training and validation data, PPG retention-index conversion scripts, figure-generation scripts, and trained Chemprop model artifacts used to support the manuscript analyses.

## Version

Current release prepared for manuscript submission: `v1.1.0` (2026-08-16).

This release updates the previous public package with the revised SI-0624 code/data/model folder. It also adds the same-condition PPG-RI control model and the revised PPG-RI training table used in the final manuscript and Supporting Information.

## Repository Contents

```text
ppg-ri-retention-prediction/
├── README.md
├── README_English_code_package.md
├── CODE_AVAILABILITY.md
├── CITATION.cff
├── LICENSE_NOTICE.md
├── environment.yml
├── requirements.txt
├── MODEL_AND_DATA_MANIFEST.md
├── Rentention_index_analysis/
│   ├── analysis_PPGRTI.py
│   ├── calculate.py
│   ├── MZ.py
│   ├── Peaktaking.py
│   ├── taking.py
│   ├── trans.py
│   ├── trans_sdf.py
│   ├── fig.py
│   ├── FIG2.py
│   ├── graph1.py ... graph4.py
│   ├── redraw_figure4_condition3_reproducible.py
│   ├── zw_graph3.py
│   ├── model_graph.py
│   ├── model_comp_fig.py
│   ├── chart_settings.json
│   ├── SDF_Output/
│   └── analysis_charts/
├── models/
│   ├── model_1/best.pt
│   ├── model_2/best.pt
│   ├── model_3/best.pt
│   ├── same_condition/best.pt
│   └── smrt_model/best.pt
└── training_data/
    ├── SECOND-CLEAN-TRAIN.csv
    ├── ppg_ri_training.csv
    ├── smrt.csv
    ├── smrt.xlsx
    ├── val_RT.csv
    ├── val_RT.xlsx
    ├── val_smrt_RT.csv
    ├── val_smrt_RT.xlsx
    ├── val_smrt_RI.csv
    ├── val_smrt_RI.xlsx
    └── SMRT/
```

Note: the directory name `Rentention_index_analysis` is retained to preserve the script paths used in the working package.

## Data Files

The processed datasets are stored in `training_data/`.

| File or folder | Records | Description |
|---|---:|---|
| `SECOND-CLEAN-TRAIN.csv` | 650,161 | Public/SMRT-derived RT records used for the large public RT baseline. |
| `SMRT/SMRT_dataset_SMRT_dataset.csv` | 79,955 | Processed SMRT dataset with PubChem CID, RT, SMILES, InChI, formula, molecular descriptors, and metadata. |
| `smrt.csv` / `smrt.xlsx` | 79,955 | SMRT-derived RT table used for model construction and conversion checks. |
| `val_RT.csv` / `val_RT.xlsx` | 269 | Laboratory validation RT table. |
| `val_smrt_RT.csv` / `val_smrt_RT.xlsx` | 277 | Validation compounds represented on the reported SMRT/literature RT scale. |
| `val_smrt_RI.csv` / `val_smrt_RI.xlsx` | 277 | Validation compounds after PPG-RI conversion. |
| `ppg_ri_training.csv` | 265 | Revised PPG-RI training table used for the final PPG-RI model workflow. |

See `MODEL_AND_DATA_MANIFEST.md` for a compact machine-readable-style inventory.

## Trained Model Artifacts

Trained Chemprop model weights are stored as `best.pt` files under `models/`.

| Folder | Artifact | Purpose |
|---|---|---|
| `models/model_1/` | `best.pt` | Direct RT model trained with laboratory RT labels. |
| `models/model_2/` | `best.pt` | Direct RT model trained under the reported SMRT/literature chromatographic condition. |
| `models/model_3/` | `best.pt` | Cross-condition PPG-RI model. |
| `models/same_condition/` | `best.pt` | Same-condition PPG-RI control model used for the Supporting Information analysis. |
| `models/smrt_model/` | `best.pt` | Large public/SMRT-derived RT baseline model. |

The model artifacts are provided for transparency and reuse. Chemprop and PyTorch version compatibility can affect direct loading; when exact reproduction is required, rebuild the environment from `environment.yml` and `requirements.txt`.

## Main Script Groups

| Script group | Main files | Purpose |
|---|---|---|
| PPG-RI calculation and conversion | `analysis_PPGRTI.py`, `calculate.py`, `trans.py` | PPG calibration, RT-to-RI conversion, RI-to-RT decoding, cross-condition transfer, and validation summaries. |
| mzML and peak processing | `MZ.py`, `Peaktaking.py`, `taking.py` | m/z calculation, mzML peak extraction, peak matching, and retention-index preprocessing. |
| Structure conversion | `trans_sdf.py` | Conversion and processing of molecular structure files. |
| Model comparison and plotting | `fig.py`, `FIG2.py`, `graph1.py` to `graph4.py`, `model_graph.py`, `model_comp_fig.py`, `redraw_figure4_condition3_reproducible.py`, `zw_graph3.py` | Figure generation and model-comparison visualization. |
| PubChem helper | `e.py` | Optional PubChem metadata lookup helper. It requires internet access and should be rerun only when external lookup is needed. |

`Rentention_index_analysis/analysis_charts/` contains generated chart outputs from previous runs.

## Software Environment

A practical environment can be created with conda:

```bash
conda env create -f environment.yml
conda activate ppg-ri-retention
pip install -r requirements.txt
```

The package uses common scientific Python packages together with Chemprop, RDKit, PyTorch dependencies installed through Chemprop, and optional LC-MS processing packages such as `pymzml` and `pyopenms`.

## Basic Reuse Workflow

1. Inspect the processed data in `training_data/`.
2. Load a trained model artifact from `models/<model_name>/best.pt` using a compatible Chemprop/PyTorch environment.
3. Use `Rentention_index_analysis/analysis_PPGRTI.py` or related scripts to calculate PPG-RIs, convert between RT and RI, and reproduce validation summaries.
4. Use the plotting scripts in `Rentention_index_analysis/` to regenerate model-comparison figures.

Some scripts were originally used as analysis notebooks/scripts during manuscript preparation and may require local path selection or minor path adaptation for independent reruns.

## Reproducibility Notes

- This repository is a reproducibility and peer-review support archive, not a polished software package.
- Several scripts are figure-specific or exploratory.
- Processed third-party/public datasets are included only to support the reported analysis. Users should cite the original data sources where applicable.
- Raw proprietary instrument files are not included unless explicitly listed in the manuscript or Supporting Information.

## Citation

If you use this package, cite the associated manuscript and this repository release. See `CITATION.cff`.

## Contact

For questions about the manuscript, data package, or model artifacts, contact the corresponding authors listed in the manuscript.
