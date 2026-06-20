# English Code Package

This folder is an English-language copy of the SI code and model/data package.
Python comments, docstrings, user-facing messages, and common labels were translated for manuscript submission.
The original SI folder was not modified.

## Main script groups

- `Rentention_index_analysis/analysis_PPGRTI.py`: PPG retention index calculation, conversion, validation, and visualization workflow.
- `Rentention_index_analysis/calculate.py`, `MZ.py`, `Peaktaking.py`, and `taking.py`: m/z calculation, mzML peak extraction, compound matching, and retention-index calculation.
- `Rentention_index_analysis/fig.py`, `FIG2.py`, `graph*.py`, `model_graph.py`, and `model_comp_fig.py`: plotting and model-comparison figure generation.
- `training_data/`: model input datasets and validation data.
- `models/`: trained model folders, configuration files, checkpoints, logs, and prediction outputs.

## Reproducibility note

Some original scripts contain local absolute paths. These paths should be replaced by paths relative to this SI package before independent reruns.
