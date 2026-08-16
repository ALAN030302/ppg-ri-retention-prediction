# English Code Package

This folder contains the English-language code, model, and processed-data package prepared from the revised SI-0624 materials.

Python comments, docstrings, user-facing messages, and common labels were translated or lightly cleaned for manuscript submission where practical. The original SI folder was not modified.

## Main Script Groups

- `Rentention_index_analysis/analysis_PPGRTI.py`: PPG retention-index calculation, conversion, validation, and visualization workflow.
- `Rentention_index_analysis/calculate.py`, `MZ.py`, `Peaktaking.py`, and `taking.py`: m/z calculation, mzML peak extraction, compound matching, and retention-index calculation.
- `Rentention_index_analysis/fig.py`, `FIG2.py`, `graph*.py`, `model_graph.py`, `model_comp_fig.py`, `redraw_figure4_condition3_reproducible.py`, and `zw_graph3.py`: plotting and model-comparison figure generation.
- `training_data/`: processed model input datasets and validation data.
- `models/`: trained Chemprop model artifacts, including the revised PPG-RI and same-condition control models.

## Reproducibility Note

Some scripts were originally used interactively during manuscript preparation and may contain local path assumptions. When rerunning independently, replace local absolute paths with paths relative to this repository.

For the full public-release description, see `README.md` and `MODEL_AND_DATA_MANIFEST.md`.
