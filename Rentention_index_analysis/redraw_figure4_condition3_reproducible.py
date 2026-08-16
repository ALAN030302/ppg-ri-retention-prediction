#!/usr/bin/env python3
"""Redraw all six Figure 4 panels from reproducible model outputs.

The three direct-RT comparator columns are read from Table S10. The PPG RI
column is replaced by the locked three-model Chemprop ensemble prediction
converted with the SMRT condition10 PPG ruler. Author-adjusted predictions
are intentionally not accepted by this script.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
BASE_TABLE = ROOT / "SI-0624" / "SI-0624" / "tables" / "TableS10_Prediction_results_summary_for_Section_3_1_3.xlsx"
REPRODUCIBLE_PPG = ROOT / "Model_ppg_ri_tuning" / "28_predicted_rt_smrt_ensemble.csv"
OUT = ROOT / "Model_ppg_ri_reproducible_audit"

COLORS = {
    "Public Database Model": "#012F48",
    "Non Literature Condition Model": "#7A0101",
    "Literature Condition Model": "#035830",
    "PPG RI Model": "#669ABA",
    "Actual": "#4C4C4C",
}
PURPLE = "#5B2A86"
COLS = {
    "Public Database Model": "Literature model (60k)",
    "Non Literature Condition Model": "Non-literature model",
    "Literature Condition Model": "Literature condition model",
    "PPG RI Model": "PPG RI Model reproducible",
}
ORDER = list(COLS)


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "font.size": 10,
            "axes.labelsize": 12,
            "axes.titlesize": 13,
            "legend.fontsize": 9,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "axes.linewidth": 0.75,
            "axes.unicode_minus": False,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def metrics(actual: pd.Series, pred: pd.Series) -> dict[str, float | int]:
    mask = actual.notna() & pred.notna()
    x = actual[mask].to_numpy(float)
    y = pred[mask].to_numpy(float)
    residual = y - x
    return {
        "n": int(len(x)),
        "r2": float(1 - np.sum(residual**2) / np.sum((x - x.mean()) ** 2)),
        "mae_min": float(np.mean(np.abs(residual))),
        "rmse_min": float(np.sqrt(np.mean(residual**2))),
        "mean_error_min": float(np.mean(residual)),
    }


def load_data() -> pd.DataFrame:
    base = pd.read_excel(BASE_TABLE)
    raw = pd.read_csv(REPRODUCIBLE_PPG)
    merged = base.merge(
        raw[["smiles", "predicted_rt", "predicted_ri"]],
        on="smiles",
        how="left",
        validate="one_to_one",
    )
    if merged["predicted_rt"].isna().any():
        raise ValueError("Not all 28 target structures matched the raw Chemprop output.")
    merged[COLS["PPG RI Model"]] = merged["predicted_rt"]
    return merged


def style_axis(ax: plt.Axes) -> None:
    ax.grid(True, color="#B8B8B8", linewidth=0.55, alpha=0.35)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.75)
        spine.set_color("#333333")
    ax.tick_params(length=3, width=0.7, color="#333333")


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.14,
        1.10,
        f"({label})",
        transform=ax.transAxes,
        color=PURPLE,
        fontsize=20,
        fontweight="bold",
        va="top",
        ha="left",
        clip_on=False,
    )


def scatter_all(ax: plt.Axes, df: pd.DataFrame) -> None:
    actual = df["true_rt"]
    lo = min(float(actual.min()), *(float(df[COLS[m]].min()) for m in ORDER)) - 0.2
    hi = max(float(actual.max()), *(float(df[COLS[m]].max()) for m in ORDER)) + 0.2
    for model in ORDER:
        pred = df[COLS[model]]
        mask = actual.notna() & pred.notna()
        x = actual[mask].to_numpy(float)
        y = pred[mask].to_numpy(float)
        ax.scatter(
            x,
            y,
            s=28,
            color=COLORS[model],
            alpha=0.72,
            edgecolor="#263238",
            linewidth=0.25,
            label=model,
            zorder=3,
        )
        fit = np.polyfit(x, y, 1)
        xr = np.linspace(x.min(), x.max(), 120)
        ax.plot(xr, np.polyval(fit, xr), "--", color=COLORS[model], linewidth=1.0, alpha=0.75)
    ax.plot([lo, hi], [lo, hi], ":", color="#9A9A9A", linewidth=1.3, label="Ideal line")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("Actual Retention Time")
    ax.set_ylabel("Predicted Retention Time")
    ax.set_title("Scatter Plots of All Models", pad=4)
    ax.legend(loc="upper left", frameon=False, handlelength=1.0, handletextpad=0.45, labelspacing=0.35, fontsize=8.5)
    style_axis(ax)


def scatter_ppg(ax: plt.Axes, df: pd.DataFrame) -> None:
    x = df["true_rt"].to_numpy(float)
    y = df[COLS["PPG RI Model"]].to_numpy(float)
    summary = metrics(df["true_rt"], df[COLS["PPG RI Model"]])
    lo = min(x.min(), y.min()) - 0.2
    hi = max(x.max(), y.max()) + 0.2
    ax.scatter(
        x,
        y,
        s=42,
        color=COLORS["PPG RI Model"],
        alpha=0.72,
        edgecolor="#54758A",
        linewidth=0.35,
        zorder=3,
    )
    fit = np.polyfit(x, y, 1)
    xr = np.linspace(x.min(), x.max(), 120)
    ax.plot(xr, np.polyval(fit, xr), "--", color="red", linewidth=2.5)
    ax.plot([lo, hi], [lo, hi], ":", color="#A8A8A8", linewidth=1.7)
    ax.text(
        0.05,
        0.95,
        f"R² = {summary['r2']:.3f}\nMAE = {summary['mae_min']:.2f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12,
    )
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("Actual Retention Time")
    ax.set_ylabel("Predicted Retention Time")
    ax.set_title("PPG RI Model", pad=4)
    style_axis(ax)


def boxplot_panel(ax: plt.Axes, df: pd.DataFrame) -> None:
    values = [np.abs(df[COLS[m]] - df["true_rt"]).to_numpy(float) for m in ORDER]
    bp = ax.boxplot(
        values,
        tick_labels=list("ABCD"),
        patch_artist=True,
        widths=0.45,
        boxprops={"linewidth": 0.75},
        medianprops={"color": "#E58C2C", "linewidth": 1.0},
        whiskerprops={"linewidth": 0.75},
        capprops={"linewidth": 0.75},
        flierprops={"marker": "o", "markersize": 3.2, "markerfacecolor": "white", "markeredgecolor": "#333333"},
    )
    for patch, model in zip(bp["boxes"], ORDER):
        patch.set_facecolor(COLORS[model])
        patch.set_alpha(0.72)
    ax.set_ylabel("Absolute Error")
    handles = [mpl.patches.Patch(facecolor=COLORS[m], edgecolor="none", alpha=0.72) for m in ORDER]
    means = [np.mean(v) for v in values]
    labels = [f"{letter}: {model}(Mean={mean:.3f})" for letter, model, mean in zip("ABCD", ORDER, means)]
    ax.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(-0.04, -0.06),
        frameon=False,
        fontsize=8.0,
        handlelength=1.8,
        labelspacing=0.38,
    )
    style_axis(ax)


def trend_panel(ax: plt.Axes, df: pd.DataFrame) -> None:
    d = df.sort_values("true_rt").reset_index(drop=True)
    x = np.arange(len(d))
    ax.plot(x, d["true_rt"], color=COLORS["Actual"], linewidth=2.2, label="Actual", zorder=5)
    for model in ORDER:
        ax.plot(x, d[COLS[model]], "--", color=COLORS[model], alpha=0.72, linewidth=1.5, label=model)
    ax.set_xlabel("Compound (sorted by Actual RT)")
    ax.set_ylabel("Retention Time")
    ax.legend(loc="upper left", frameon=False, fontsize=8.0, handlelength=2.2, labelspacing=0.35)
    style_axis(ax)


def error_distribution(ax: plt.Axes, df: pd.DataFrame) -> None:
    all_errors = []
    for model in ORDER:
        error = (df[COLS[model]] - df["true_rt"]).to_numpy(float)
        all_errors.extend(error.tolist())
        ax.hist(
            error,
            bins=30,
            density=True,
            alpha=0.45,
            color=COLORS[model],
            edgecolor="#222222",
            linewidth=0.35,
            label=model,
        )
    ax.axvline(0, color="red", linestyle="--", linewidth=1.5, label="Zero error")
    ax.set_xlim(min(all_errors) - 0.2, max(all_errors) + 0.2)
    ax.set_xlabel("Prediction Error (Predicted - Actual)")
    ax.set_ylabel("Density")
    ax.set_title("Error Distribution", pad=4)
    ax.legend(loc="upper left", frameon=False, fontsize=7.5, handlelength=1.6, labelspacing=0.30)
    style_axis(ax)


def residual_panel(ax: plt.Axes, df: pd.DataFrame) -> None:
    for model in ORDER:
        residual = df[COLS[model]] - df["true_rt"]
        ax.scatter(
            df["true_rt"],
            residual,
            s=19,
            color=COLORS[model],
            alpha=0.62,
            edgecolor="none",
            label=model,
        )
    ax.axhline(0, color="red", linestyle="--", linewidth=1.5)
    ax.set_xlabel("Actual Retention Time")
    ax.set_ylabel("Residual")
    ax.set_title("Residual Plot", pad=4)
    ax.legend(loc="lower left", frameon=False, fontsize=7.5, handletextpad=0.35, labelspacing=0.30)
    style_axis(ax)


def save_figure(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    fig.savefig(
        stem.with_suffix(".tiff"),
        dpi=600,
        format="tiff",
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


def main() -> None:
    configure_style()
    df = load_data()
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "Figure4_condition3_to_SMRT_ensemble_REPRODUCIBLE_source_data.csv", index=False)

    fig = plt.figure(figsize=(8.3, 11.07))
    grid = fig.add_gridspec(
        3,
        2,
        left=0.09,
        right=0.97,
        bottom=0.06,
        top=0.98,
        wspace=0.28,
        hspace=0.48,
        height_ratios=[1.0, 0.98, 1.0],
    )
    axes = [fig.add_subplot(grid[i, j]) for i in range(3) for j in range(2)]
    scatter_all(axes[0], df)
    scatter_ppg(axes[1], df)
    boxplot_panel(axes[2], df)
    trend_panel(axes[3], df)
    error_distribution(axes[4], df)
    residual_panel(axes[5], df)
    for label, ax in zip("ABCDEF", axes):
        panel_label(ax, label)
    stem = OUT / "Figure4_condition3_to_SMRT_ensemble_REPRODUCIBLE"
    save_figure(fig, stem)

    summary = {model: metrics(df["true_rt"], df[column]) for model, column in COLS.items()}
    summary["provenance"] = {
        "comparator_source": str(BASE_TABLE),
        "ppg_source": str(REPRODUCIBLE_PPG),
        "ppg_note": "Locked three-model Chemprop ensemble RI prediction converted with the SMRT condition10 PPG ruler.",
    }
    (OUT / "Figure4_condition3_to_SMRT_ensemble_REPRODUCIBLE_metrics.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    with Image.open(stem.with_suffix(".tiff")) as image:
        dpi_value = image.info.get("dpi")
        dpi = [float(value) for value in dpi_value] if dpi_value else None
        print(json.dumps({"tiff_pixels": [image.width, image.height], "tiff_dpi": dpi, "metrics": summary}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
