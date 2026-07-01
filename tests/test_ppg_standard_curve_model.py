from __future__ import annotations

import sys
import types
import unittest
from collections import namedtuple
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Rentention_index_analysis"))


def _install_import_stubs() -> None:
    scipy = types.ModuleType("scipy")
    stats = types.ModuleType("scipy.stats")
    result_type = namedtuple("LinregressResult", "slope intercept rvalue pvalue stderr")

    def linregress(x, y):
        x_arr = np.asarray(x, dtype=float)
        y_arr = np.asarray(y, dtype=float)
        slope, intercept = np.polyfit(x_arr, y_arr, 1)
        y_hat = slope * x_arr + intercept
        ss_res = float(np.sum((y_arr - y_hat) ** 2))
        ss_tot = float(np.sum((y_arr - np.mean(y_arr)) ** 2))
        r_squared = 1.0 - ss_res / ss_tot if ss_tot else 1.0
        rvalue = float(np.sqrt(max(r_squared, 0.0)))
        return result_type(float(slope), float(intercept), rvalue, 0.0, 0.0)

    stats.linregress = linregress
    scipy.stats = stats
    sys.modules.setdefault("scipy", scipy)
    sys.modules.setdefault("scipy.stats", stats)

    seaborn = types.ModuleType("seaborn")
    sys.modules.setdefault("seaborn", seaborn)

    matplotlib = types.ModuleType("matplotlib")
    pyplot = types.ModuleType("matplotlib.pyplot")
    backends = types.ModuleType("matplotlib.backends")
    backend_tkagg = types.ModuleType("matplotlib.backends.backend_tkagg")
    figure_mod = types.ModuleType("matplotlib.figure")

    class Figure:
        pass

    pyplot.Figure = Figure
    backend_tkagg.FigureCanvasTkAgg = object
    backend_tkagg.NavigationToolbar2Tk = object
    figure_mod.Figure = Figure
    matplotlib.pyplot = pyplot
    sys.modules.setdefault("matplotlib", matplotlib)
    sys.modules.setdefault("matplotlib.pyplot", pyplot)
    sys.modules.setdefault("matplotlib.backends", backends)
    sys.modules.setdefault("matplotlib.backends.backend_tkagg", backend_tkagg)
    sys.modules.setdefault("matplotlib.figure", figure_mod)


_install_import_stubs()

from analysis_PPGRTI import PPGIndexCalculator  # noqa: E402


class PPGStandardCurveModelTest(unittest.TestCase):
    def setUp(self) -> None:
        self.calculator = PPGIndexCalculator()
        self.calculator.ppg_data["condition1"] = pd.DataFrame(
            {
                "degree_of_polymerization": [3, 4, 5, 6],
                "retention_time": [2.0, 3.0, 4.0, 5.0],
            }
        )
        self.calculator.ppg_data["condition2"] = pd.DataFrame(
            {
                "degree_of_polymerization": [3, 4, 5, 6],
                "retention_time": [4.0, 6.0, 8.0, 10.0],
            }
        )
        self.calculator.compound_data["validation_condition1"] = pd.DataFrame(
            {"compound_name": ["A"], "retention_time": [3.5]}
        )
        self.calculator.compound_data["validation_condition2"] = pd.DataFrame(
            {"compound_name": ["A"], "retention_time": [7.2]}
        )

    def test_standard_curve_defaults_to_linear_n_model(self) -> None:
        success, _ = self.calculator.fit_standard_curve("condition1")
        self.assertTrue(success)
        self.assertEqual(self.calculator.standard_curves["condition1"]["model_type"], "linear")

    def test_regression_index_calculation_uses_linear_curve_when_missing(self) -> None:
        success, _ = self.calculator.calculate_ppg_index("condition1", method="regression")
        self.assertTrue(success)
        self.assertEqual(self.calculator.standard_curves["condition1"]["model_type"], "linear")

    def test_cross_condition_conversion_uses_linear_curve_when_missing(self) -> None:
        self.calculator.ppg_indices["condition1"] = {
            "indices": {
                "validation_condition1": pd.DataFrame(
                    {"compound_name": ["A"], "calculatePPGindex": [450.0]}
                )
            }
        }
        result, _ = self.calculator.convert_ppg_index_to_rt("condition1", "condition2")
        self.assertFalse(result.empty)
        self.assertEqual(self.calculator.standard_curves["condition2"]["model_type"], "linear")
        self.assertAlmostEqual(float(result.loc[0, "condition2_predicted_RT"]), 7.0)


if __name__ == "__main__":
    unittest.main()
