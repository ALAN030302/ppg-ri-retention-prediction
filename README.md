# PPG-RI Retention Prediction Supporting Code and Data

This repository contains the English-language code, model configuration files, trained model artifacts, processed training/validation data, and analysis scripts supporting the manuscript:

> Cross-Platform Retention Time Prediction of Environmental Contaminants: A Polypropylene Glycol Retention Index-Driven Graph Neural Network

The package supports the study workflow for converting chromatographic retention times (RTs) to polypropylene glycol retention indices (PPG-RIs), training and evaluating graph neural network models, comparing RT- and RI-based prediction strategies, and generating model-comparison figures.

## Repository Status

This package is intended as a reproducibility and peer-review support archive. The code was translated and lightly cleaned for public release. Several scripts remain legacy analysis scripts that were used interactively during the study, so exact reruns may require selecting local input files or adapting paths to the user's environment.

Current public record:

- GitHub repository: <https://github.com/ALAN030302/ppg-ri-retention-prediction>
- Versioned GitHub release: `v1.0.0`

## Directory Structure

```text
ppg-ri-retention-prediction/
├── README.md
├── README_English_code_package.md
├── CODE_AVAILABILITY.md
├── CITATION.cff
├── LICENSE_NOTICE.md
├── environment.yml
├── requirements.txt
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
│   ├── model_graph.py
│   ├── model_comp_fig.py
│   ├── chart_settings.json
│   ├── SDF_Output/
│   └── analysis_charts/
├── models/
│   ├── model_1/config.toml
│   ├── model_2/config.toml
│   ├── model_3/config.toml
│   └── smrt_model/config.toml
└── training_data/
    ├── SECOND-CLEAN-TRAIN.csv
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

Note: the folder name `Rentention_index_analysis` is retained to preserve script paths from the working package.

## File Groups

### `training_data/`

Processed datasets used for model training, validation, and comparison.

| File or folder | Description |
|---|---|
| `SECOND-CLEAN-TRAIN.csv` | Large SMRT-derived RT training dataset with SMILES and retention-time values. |
| `smrt.csv`, `smrt.xlsx` | SMRT-derived RT data converted to minutes. |
| `val_RT.csv`, `val_RT.xlsx` | Validation compounds with experimentally measured RT values. |
| `val_smrt_RT.csv`, `val_smrt_RT.xlsx` | Validation compounds represented on the SMRT/literature RT scale. |
| `val_smrt_RI.csv`, `val_smrt_RI.xlsx` | Validation compounds after conversion to PPG-RI values. |
| `SMRT/` | Source/processed SMRT dataset files and metadata summaries. |

Typical columns include:

| Column | Meaning |
|---|---|
| `SMILES` / `smiles` | Molecular structure in SMILES format. |
| `rt` / `RT` / `Retention Time` | Retention time or retention-index target, depending on the file. |
| `RETENTION_TIME` | Original SMRT retention-time field. |

### `models/`

Chemprop model configuration files. The public-release configs use repository-relative paths.

| Model folder | Input file | Purpose |
|---|---|---|
| `model_1/` | `training_data/val_RT.csv` | Experimental RT model. |
| `model_2/` | `training_data/val_smrt_RT.csv` | SMRT/literature RT model. |
| `model_3/` | `training_data/val_smrt_RI.csv` | PPG-RI model. |
| `smrt_model/` | `training_data/SECOND-CLEAN-TRAIN.csv` | Large SMRT-trained baseline model. |

### `Rentention_index_analysis/`

Analysis and figure-generation scripts.

| Script group | Main files | Purpose |
|---|---|---|
| PPG-RI calculation and conversion | `analysis_PPGRTI.py`, `calculate.py`, `trans.py` | PPG calibration, RT-to-RI conversion, cross-condition transfer, and validation summaries. |
| mzML / peak processing | `MZ.py`, `Peaktaking.py`, `taking.py` | m/z calculation, mzML peak extraction, compound matching, and retention-index preprocessing. |
| Structure conversion | `trans_sdf.py` | Conversion and processing of molecular structure files. |
| Model comparison and plotting | `fig.py`, `FIG2.py`, `graph1.py` to `graph4.py`, `model_graph.py`, `model_comp_fig.py` | Figure generation and model-comparison visualization. |
| PubChem helper | `e.py` | Optional PubChem-based metadata retrieval helper. It requires internet access and should be rerun only when external lookup is needed. |

`analysis_charts/` contains generated chart outputs from previous runs.

## Software Environment

The exact package versions used during model training should be recorded from the final author environment if available. The files below provide a practical audit/rerun starting point:

- `environment.yml` for a conda-based environment.
- `requirements.txt` for pip-installable Python dependencies.

Example setup:

```bash
conda env create -f environment.yml
conda activate ppg-ri-retention
pip install -r requirements.txt
```

Chemprop and PyTorch version compatibility can vary. For exact reproduction, use the Chemprop major version used for the final manuscript analysis.

## Basic Rerun Workflow

### 1. Check the Chemprop configs

The `models/*/config.toml` files now use repository-relative paths. Confirm the `data-path` and `output-dir` fields before training:

```toml
data-path = "training_data/val_smrt_RI.csv"
output-dir = "models/model_3"
```

### 2. Train or rerun Chemprop models

Example command:

```bash
chemprop train --config-path models/model_3/config.toml
```

Repeat for `model_1`, `model_2`, `model_3`, and `smrt_model` as needed.

### 3. Run PPG-RI analysis scripts

Some scripts are interactive and may open a file-picker window. Select inputs from `training_data/` or from the processed workbook generated by the preceding step.

Typical entry points:

```bash
python Rentention_index_analysis/analysis_PPGRTI.py
python Rentention_index_analysis/calculate.py
python Rentention_index_analysis/trans.py
```

### 4. Generate comparison figures

```bash
python Rentention_index_analysis/fig.py
python Rentention_index_analysis/model_comp_fig.py
python Rentention_index_analysis/graph1.py
python Rentention_index_analysis/graph2.py
python Rentention_index_analysis/graph3.py
python Rentention_index_analysis/graph4.py
```

Generated plots may be written to `Rentention_index_analysis/analysis_charts/` or to a path specified inside each script.

## Reproducibility Notes

- The package contains processed data, model configurations, trained model artifacts, and analysis scripts used to support the manuscript.
- Some scripts are exploratory or figure-specific rather than a single automated pipeline.
- Some helper scripts require local intermediate files that are not part of the final manuscript dataset.
- Keep raw, processed, and figure-source data separate in any future repository expansion.

## Contact

For questions about the manuscript and data package, contact the corresponding authors listed in the manuscript.
