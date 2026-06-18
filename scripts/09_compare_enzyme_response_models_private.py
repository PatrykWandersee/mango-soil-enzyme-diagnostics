from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


DATA_PATH = Path("data/processed/private/enzyme_diagnostics_private.csv")
SUMMARY_PATH = Path("tables/private/main_enzyme_threshold_summary_private.csv")
OUT_PATH = Path("tables/private/enzyme_response_model_comparison_private.csv")

OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

RESPONSE = "Prod_rel_pct"

INDICATORS = [
    "GMea",
    "Beta_glic",
    "Arilsulf",
    "Ativ_Enzimat_Total",
    "GMea_per_Clay",
    "Beta_por_Argila",
]


def linear_plateau(x, plateau, slope, breakpoint):
    return plateau - slope * np.maximum(0, breakpoint - x)


def quadratic_plateau(x, plateau, curvature, breakpoint):
    return plateau - curvature * np.maximum(0, breakpoint - x) ** 2


def model_metrics(y, yhat, k):
    n = len(y)
    residuals = y - yhat
    rss = float(np.sum(residuals ** 2))
    tss = float(np.sum((y - np.mean(y)) ** 2))
    rmse = float(np.sqrt(rss / n))
    r2 = float(1 - rss / tss) if tss > 0 else np.nan

    # AIC for Gaussian residuals up to an additive constant.
    rss_for_aic = max(rss, 1e-12)
    aic = float(n * np.log(rss_for_aic / n) + 2 * k)

    if n > k + 1:
        aicc = float(aic + (2 * k * (k + 1)) / (n - k - 1))
    else:
        aicc = np.nan

    return {
        "n": n,
        "k": k,
        "rss": rss,
        "rmse": rmse,
        "r2": r2,
        "aic": aic,
        "aicc": aicc,
    }


def fit_linear(x, y):
    slope, intercept = np.polyfit(x, y, 1)
    yhat = intercept + slope * x

    row = {
        "model": "linear",
        "intercept_or_plateau": intercept,
        "linear_slope": slope,
        "quadratic_term": np.nan,
        "breakpoint_or_vertex": np.nan,
        "fit_status": "ok",
    }
    row.update(model_metrics(y, yhat, k=2))
    return row


def fit_quadratic(x, y):
    a, b, c = np.polyfit(x, y, 2)
    yhat = a * x ** 2 + b * x + c

    if a < 0:
        vertex = -b / (2 * a)
    else:
        vertex = np.nan

    row = {
        "model": "quadratic",
        "intercept_or_plateau": c,
        "linear_slope": b,
        "quadratic_term": a,
        "breakpoint_or_vertex": vertex,
        "fit_status": "ok",
    }
    row.update(model_metrics(y, yhat, k=3))
    return row


def get_initial_breakpoint(indicator, x, summary):
    if summary is not None and indicator in set(summary["indicator"]):
        value = summary.loc[
            summary["indicator"] == indicator,
            "threshold_tree",
        ].iloc[0]
        if np.isfinite(value):
            return float(np.clip(value, np.min(x), np.max(x)))

    return float(np.median(x))


def fit_linear_plateau(indicator, x, y, summary):
    x_min = float(np.min(x))
    x_max = float(np.max(x))
    breakpoint0 = get_initial_breakpoint(indicator, x, summary)

    plateau0 = float(np.percentile(y, 75))
    slope0 = max((np.percentile(y, 75) - np.percentile(y, 25)) / max(x_max - x_min, 1e-9), 1e-6)

    bounds_lower = [-50.0, 0.0, x_min]
    bounds_upper = [150.0, np.inf, x_max]

    try:
        params, _ = curve_fit(
            linear_plateau,
            x,
            y,
            p0=[plateau0, slope0, breakpoint0],
            bounds=(bounds_lower, bounds_upper),
            maxfev=20000,
        )
        plateau, slope, breakpoint = params
        yhat = linear_plateau(x, plateau, slope, breakpoint)

        row = {
            "model": "linear_plateau",
            "intercept_or_plateau": plateau,
            "linear_slope": slope,
            "quadratic_term": np.nan,
            "breakpoint_or_vertex": breakpoint,
            "fit_status": "ok",
        }
        row.update(model_metrics(y, yhat, k=3))
        return row

    except Exception as exc:
        return {
            "model": "linear_plateau",
            "intercept_or_plateau": np.nan,
            "linear_slope": np.nan,
            "quadratic_term": np.nan,
            "breakpoint_or_vertex": np.nan,
            "fit_status": f"failed: {exc}",
            "n": len(y),
            "k": 3,
            "rss": np.nan,
            "rmse": np.nan,
            "r2": np.nan,
            "aic": np.nan,
            "aicc": np.nan,
        }


def fit_quadratic_plateau(indicator, x, y, summary):
    x_min = float(np.min(x))
    x_max = float(np.max(x))
    breakpoint0 = get_initial_breakpoint(indicator, x, summary)

    plateau0 = float(np.percentile(y, 75))
    curvature0 = max(
        (np.percentile(y, 75) - np.percentile(y, 25))
        / max((breakpoint0 - x_min) ** 2, 1e-9),
        1e-8,
    )

    bounds_lower = [-50.0, 0.0, x_min]
    bounds_upper = [150.0, np.inf, x_max]

    try:
        params, _ = curve_fit(
            quadratic_plateau,
            x,
            y,
            p0=[plateau0, curvature0, breakpoint0],
            bounds=(bounds_lower, bounds_upper),
            maxfev=20000,
        )
        plateau, curvature, breakpoint = params
        yhat = quadratic_plateau(x, plateau, curvature, breakpoint)

        row = {
            "model": "quadratic_plateau",
            "intercept_or_plateau": plateau,
            "linear_slope": np.nan,
            "quadratic_term": curvature,
            "breakpoint_or_vertex": breakpoint,
            "fit_status": "ok",
        }
        row.update(model_metrics(y, yhat, k=3))
        return row

    except Exception as exc:
        return {
            "model": "quadratic_plateau",
            "intercept_or_plateau": np.nan,
            "linear_slope": np.nan,
            "quadratic_term": np.nan,
            "breakpoint_or_vertex": np.nan,
            "fit_status": f"failed: {exc}",
            "n": len(y),
            "k": 3,
            "rss": np.nan,
            "rmse": np.nan,
            "r2": np.nan,
            "aic": np.nan,
            "aicc": np.nan,
        }


def main():
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Input data not found: {DATA_PATH}")

    df = pd.read_csv(DATA_PATH)

    summary = None
    if SUMMARY_PATH.exists():
        summary = pd.read_csv(SUMMARY_PATH)

    rows = []

    for indicator in INDICATORS:
        sub = df[[indicator, RESPONSE]].dropna().copy()

        x = sub[indicator].to_numpy(dtype=float)
        y = sub[RESPONSE].to_numpy(dtype=float)

        fitted_rows = [
            fit_linear(x, y),
            fit_quadratic(x, y),
            fit_linear_plateau(indicator, x, y, summary),
            fit_quadratic_plateau(indicator, x, y, summary),
        ]

        for row in fitted_rows:
            row["indicator"] = indicator
            rows.append(row)

    results = pd.DataFrame(rows)

    results["best_aicc_for_indicator"] = results.groupby("indicator")["aicc"].transform("min")
    results["delta_aicc"] = results["aicc"] - results["best_aicc_for_indicator"]

    results = results[
        [
            "indicator",
            "model",
            "fit_status",
            "n",
            "k",
            "r2",
            "rmse",
            "aic",
            "aicc",
            "delta_aicc",
            "rss",
            "intercept_or_plateau",
            "linear_slope",
            "quadratic_term",
            "breakpoint_or_vertex",
        ]
    ].sort_values(["indicator", "delta_aicc", "aicc"])

    results.to_csv(OUT_PATH, index=False)

    print("Private enzyme response model comparison completed.")
    print(f"Output file: {OUT_PATH}")

    print("\nBest model by AICc for each indicator:")
    best = results.sort_values(["indicator", "delta_aicc"]).groupby("indicator").head(1)
    display_cols = [
        "indicator",
        "model",
        "r2",
        "rmse",
        "aicc",
        "delta_aicc",
        "intercept_or_plateau",
        "breakpoint_or_vertex",
        "fit_status",
    ]
    print(best[display_cols].to_string(index=False))

    print("\nFull model comparison:")
    display_cols = [
        "indicator",
        "model",
        "r2",
        "rmse",
        "aicc",
        "delta_aicc",
        "breakpoint_or_vertex",
        "fit_status",
    ]
    print(results[display_cols].to_string(index=False))


if __name__ == "__main__":
    main()
