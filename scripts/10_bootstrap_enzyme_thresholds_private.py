from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit
from sklearn.tree import DecisionTreeRegressor


DATA_PATH = Path("data/processed/private/enzyme_diagnostics_private.csv")
OUT_DIR = Path("tables/private")
OUT_DIR.mkdir(parents=True, exist_ok=True)

BOOT_VALUES_PATH = OUT_DIR / "enzyme_threshold_bootstrap_values_private.csv"
BOOT_SUMMARY_PATH = OUT_DIR / "enzyme_threshold_bootstrap_summary_private.csv"

RESPONSE = "Prod_rel_pct"
YIELD_CUTOFF = 70.0
MIN_LEAF_SIZE = 8
N_BOOT = 500
RANDOM_SEED = 42

INDICATORS = [
    "GMea",
    "Beta_glic",
    "Arilsulf",
]


def candidate_thresholds(values):
    unique_values = np.sort(pd.Series(values).dropna().unique())
    if len(unique_values) < 2:
        return np.array([])
    return (unique_values[:-1] + unique_values[1:]) / 2


def balanced_accuracy_binary(y_true, y_pred):
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)

    classes = np.unique(y_true)
    if len(classes) < 2:
        return np.nan

    sensitivity = np.mean(y_pred[y_true == 1] == 1)
    specificity = np.mean(y_pred[y_true == 0] == 0)

    return 0.5 * (sensitivity + specificity)


def cate_nelson_threshold(x, y, yield_cutoff):
    sub = pd.DataFrame({"x": x, "y": y}).dropna()

    if sub.shape[0] < 10:
        return np.nan

    y_class = (sub["y"] >= yield_cutoff).astype(int).to_numpy()

    if len(np.unique(y_class)) < 2:
        return np.nan

    thresholds = candidate_thresholds(sub["x"])
    if len(thresholds) == 0:
        return np.nan

    best_threshold = np.nan
    best_balanced_accuracy = -np.inf
    best_accuracy = -np.inf

    for threshold in thresholds:
        x_class = (sub["x"].to_numpy() >= threshold).astype(int)

        balanced_accuracy = balanced_accuracy_binary(y_class, x_class)
        accuracy = np.mean(x_class == y_class)

        if np.isnan(balanced_accuracy):
            continue

        if (
            balanced_accuracy > best_balanced_accuracy
            or (
                balanced_accuracy == best_balanced_accuracy
                and accuracy > best_accuracy
            )
        ):
            best_balanced_accuracy = balanced_accuracy
            best_accuracy = accuracy
            best_threshold = threshold

    return float(best_threshold)


def tree_threshold(x, y):
    sub = pd.DataFrame({"x": x, "y": y}).dropna()

    if sub.shape[0] < 2 * MIN_LEAF_SIZE:
        return np.nan

    model = DecisionTreeRegressor(
        max_depth=1,
        min_samples_leaf=MIN_LEAF_SIZE,
        random_state=42,
    )

    model.fit(sub[["x"]].to_numpy(), sub["y"].to_numpy())
    threshold = model.tree_.threshold[0]

    if threshold == -2:
        return np.nan

    return float(threshold)


def linear_plateau(x, plateau, slope, breakpoint):
    return plateau - slope * np.maximum(0, breakpoint - x)


def linear_plateau_breakpoint(x, y):
    sub = pd.DataFrame({"x": x, "y": y}).dropna()

    if sub.shape[0] < 12:
        return np.nan

    x_arr = sub["x"].to_numpy(dtype=float)
    y_arr = sub["y"].to_numpy(dtype=float)

    x_min = float(np.min(x_arr))
    x_max = float(np.max(x_arr))

    if x_max <= x_min:
        return np.nan

    plateau0 = float(np.percentile(y_arr, 75))
    breakpoint0 = float(np.median(x_arr))
    slope0 = max(
        (np.percentile(y_arr, 75) - np.percentile(y_arr, 25))
        / max(x_max - x_min, 1e-9),
        1e-6,
    )

    try:
        params, _ = curve_fit(
            linear_plateau,
            x_arr,
            y_arr,
            p0=[plateau0, slope0, breakpoint0],
            bounds=([-50.0, 0.0, x_min], [150.0, np.inf, x_max]),
            maxfev=20000,
        )
        return float(params[2])

    except Exception:
        return np.nan


def high_yield_iqr(x, y, yield_cutoff):
    sub = pd.DataFrame({"x": x, "y": y}).dropna()
    high = sub[sub["y"] >= yield_cutoff]

    if high.shape[0] < 3:
        return {
            "n_high_yield": int(high.shape[0]),
            "high_yield_p25": np.nan,
            "high_yield_median": np.nan,
            "high_yield_p75": np.nan,
        }

    q = high["x"].quantile([0.25, 0.50, 0.75])

    return {
        "n_high_yield": int(high.shape[0]),
        "high_yield_p25": float(q.loc[0.25]),
        "high_yield_median": float(q.loc[0.50]),
        "high_yield_p75": float(q.loc[0.75]),
    }


def spearman_rho(x, y):
    sub = pd.DataFrame({"x": x, "y": y}).dropna()

    if sub.shape[0] < 3:
        return np.nan

    result = stats.spearmanr(sub["x"], sub["y"])

    return float(result.statistic)


def summarize_bootstrap(values, original_estimate):
    values = pd.Series(values).dropna()

    if values.empty:
        return {
            "original_estimate": original_estimate,
            "valid_bootstrap_n": 0,
            "bootstrap_median": np.nan,
            "ci_2_5": np.nan,
            "ci_97_5": np.nan,
            "bootstrap_mean": np.nan,
            "bootstrap_sd": np.nan,
        }

    return {
        "original_estimate": original_estimate,
        "valid_bootstrap_n": int(values.shape[0]),
        "bootstrap_median": float(values.median()),
        "ci_2_5": float(values.quantile(0.025)),
        "ci_97_5": float(values.quantile(0.975)),
        "bootstrap_mean": float(values.mean()),
        "bootstrap_sd": float(values.std(ddof=1)),
    }


def compute_metrics(df, indicator):
    x = df[indicator].to_numpy(dtype=float)
    y = df[RESPONSE].to_numpy(dtype=float)

    iqr = high_yield_iqr(x, y, YIELD_CUTOFF)

    return {
        "spearman_rho": spearman_rho(x, y),
        "cate_nelson_threshold_70": cate_nelson_threshold(x, y, YIELD_CUTOFF),
        "tree_threshold": tree_threshold(x, y),
        "linear_plateau_breakpoint": linear_plateau_breakpoint(x, y),
        **iqr,
    }


def main():
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Input data not found: {DATA_PATH}")

    df = pd.read_csv(DATA_PATH)

    missing = [col for col in [RESPONSE] + INDICATORS if col not in df.columns]
    if missing:
        raise ValueError("Missing columns: " + ", ".join(missing))

    rng = np.random.default_rng(RANDOM_SEED)

    original_rows = {}
    for indicator in INDICATORS:
        original_rows[indicator] = compute_metrics(df, indicator)

    boot_rows = []

    for boot_id in range(1, N_BOOT + 1):
        sample_indices = rng.integers(0, df.shape[0], df.shape[0])
        boot_df = df.iloc[sample_indices].reset_index(drop=True)

        for indicator in INDICATORS:
            metrics = compute_metrics(boot_df, indicator)
            boot_rows.append(
                {
                    "bootstrap_id": boot_id,
                    "indicator": indicator,
                    **metrics,
                }
            )

    boot_values = pd.DataFrame(boot_rows)
    boot_values.to_csv(BOOT_VALUES_PATH, index=False)

    summary_rows = []
    metric_columns = [
        "spearman_rho",
        "cate_nelson_threshold_70",
        "tree_threshold",
        "linear_plateau_breakpoint",
        "n_high_yield",
        "high_yield_p25",
        "high_yield_median",
        "high_yield_p75",
    ]

    for indicator in INDICATORS:
        indicator_boot = boot_values[boot_values["indicator"] == indicator]

        for metric in metric_columns:
            original_estimate = original_rows[indicator][metric]
            summary = summarize_bootstrap(
                indicator_boot[metric],
                original_estimate=original_estimate,
            )
            summary_rows.append(
                {
                    "indicator": indicator,
                    "metric": metric,
                    **summary,
                }
            )

    boot_summary = pd.DataFrame(summary_rows)
    boot_summary.to_csv(BOOT_SUMMARY_PATH, index=False)

    print("Bootstrap analysis completed.")
    print(f"Bootstrap iterations: {N_BOOT}")
    print(f"Bootstrap values: {BOOT_VALUES_PATH}")
    print(f"Bootstrap summary: {BOOT_SUMMARY_PATH}")

    print("\nKey bootstrap summary:")
    key_metrics = [
        "spearman_rho",
        "cate_nelson_threshold_70",
        "tree_threshold",
        "linear_plateau_breakpoint",
        "high_yield_p25",
        "high_yield_p75",
    ]
    key = boot_summary[boot_summary["metric"].isin(key_metrics)].copy()

    display_cols = [
        "indicator",
        "metric",
        "original_estimate",
        "valid_bootstrap_n",
        "bootstrap_median",
        "ci_2_5",
        "ci_97_5",
    ]

    print(key[display_cols].to_string(index=False))


if __name__ == "__main__":
    main()
