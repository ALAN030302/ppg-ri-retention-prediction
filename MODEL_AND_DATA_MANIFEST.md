# Model and Data Manifest

Release: `v1.1.0`  
Prepared: 2026-08-16

## Processed Data

| Path | Records | Key columns | Notes |
|---|---:|---|---|
| `training_data/SECOND-CLEAN-TRAIN.csv` | 650,161 | `Smiles`, `Retention Time` | Large public/SMRT-derived RT training table. |
| `training_data/SMRT/SMRT_dataset_SMRT_dataset.csv` | 79,955 | `PUBCHEM_COMPOUND_CID`, `RETENTION_TIME`, `ID`, `SMILES`, `InChI` | Processed SMRT source table with molecular metadata. |
| `training_data/smrt.csv` | 79,955 | `SMILES`, `rt` | SMRT-derived RT table used in model construction. |
| `training_data/val_RT.csv` | 269 | `smiles`, `rt` | Laboratory validation RT table. |
| `training_data/val_smrt_RT.csv` | 277 | `smiles`, `rt` | Validation compounds on the reported SMRT/literature RT scale. |
| `training_data/val_smrt_RI.csv` | 277 | `smiles`, `RT` | Validation compounds after PPG-RI conversion. |
| `training_data/ppg_ri_training.csv` | 265 | `smiles`, `ri`, `source_count`, `source_compound_ids` | Revised PPG-RI training table used for the final PPG-RI workflow. |

Spreadsheet copies are provided for selected CSV files to support manual inspection.

## Model Artifacts

| Path | Approximate size | Description |
|---|---:|---|
| `models/model_1/best.pt` | 1.22 MB | Direct RT model trained with laboratory RT labels. |
| `models/model_2/best.pt` | 1.22 MB | Direct RT model trained under the reported SMRT/literature chromatographic condition. |
| `models/model_3/best.pt` | 1.23 MB | Cross-condition PPG-RI model. |
| `models/same_condition/best.pt` | 1.22 MB | Same-condition PPG-RI control model. |
| `models/smrt_model/best.pt` | 1.22 MB | Large public/SMRT-derived RT baseline model. |

## Main Analysis Scripts

| Path | Role |
|---|---|
| `Rentention_index_analysis/analysis_PPGRTI.py` | Main PPG-RI calculation, standard-curve fitting, conversion, validation, and visualization workflow. |
| `Rentention_index_analysis/calculate.py` | Supporting PPG-RI and RT calculation workflow. |
| `Rentention_index_analysis/MZ.py`, `Peaktaking.py`, `taking.py` | m/z, peak extraction, and compound-matching utilities. |
| `Rentention_index_analysis/trans.py`, `trans_sdf.py` | Data and molecular-structure conversion helpers. |
| `Rentention_index_analysis/fig.py`, `FIG2.py`, `graph1.py`, `graph2.py`, `graph3.py`, `graph4.py`, `model_graph.py`, `model_comp_fig.py` | Manuscript and Supporting Information figure-generation scripts. |
| `Rentention_index_analysis/redraw_figure4_condition3_reproducible.py`, `zw_graph3.py` | Revised figure-generation utilities retained from the final SI package. |

## Excluded Local Files

The public release excludes local operating-system and IDE metadata, including `.DS_Store`, `.idea/`, and `__pycache__/`.
