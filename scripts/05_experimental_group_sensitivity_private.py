from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import balanced_accuracy_score
from sklearn.tree import DecisionTreeRegressor


INPUT_PATH = Path("data/processed/private/enzyme_diagnostics_private.csv")
TABLE_DIR = Path("tables/private")
TABLE_DIR.mkdir(parents=True, exist_ok=True)

RESPONSE = "Prod_rel_pct"
MAIN_YIELD_CUTOFF = 70
MIN_LEAF_SIZE = 8

INDICATORS = [
    "GMea",
    "Beta_glic",
    "Arilsulf",
    "Ativ_Enzimat_Total",
    "GMea_per_Clay",
    "log_GMea_per_Clay",
    "Beta_por_Argila",
    "Aril_por_Argila",
    "Ativ_Enzim_por_MO",
    "qGMea",
]


def detect_experimental_mask(df):
    mask = pd.Series(False, index=df.index)

    for col in ["Fazenda", "Area"]:
        if col in df.columns:
            values = df[col].astype(str).str.lower()
            mask = mask | values.str.contains("exp", na=False)
            mask = mask | values.str.contains("experimental", na=False)

    return mask


def candidate_thresholds(values):
    unique_values = np.sort(pd.Series(values).dropna().unique())
    if len(unique_values) < 2:
        return np.array([])
    return (unique_values[:-1] + unique_values[1:]) / 2


def cate_nelson_fixed_y(df, indicator, response, yield_cutoff):
    sub = df[[indicator, response]].dropna().copy()

    if sub.shape[0] < 10:
        return None

    y_class = (sub[response] >= yield_cutoff).astype(int)
    thresholds = candidate_thresholds(sub[indicator])

    best = None

    for threshold in thresholds:
        x_class = (sub[indicator] >= threshold).astype(int)

        accuracy = (x_class == y_class).mean()
        balanced_accuracy = balanced_accuracy_score(y_class, x_class)

        current = {
            "indicator": indicator,
            "n": int(sub.shape[0]),
            "yield_cutoff": yield_cutoff,
            "threshold_cate_nelson": float(threshold),
            "accuracy": float(accuracy),
            "balanced_accuracy": float(balanced_accuracy),
            "low_indicator_low_yield": int(((x_class == 0) & (y_class == 0)).sum()),
            "low_indicator_high_yield": int(((x_class == 0) & (y_class == 1)).sum()),
            "high_indicator_low_yield": int(((x_class == 1) & (y_class == 0)).sum()),
            "high_indicator_high_yield": int(((x_class == 1) & (y_class == 1)).sum()),
        }

        if best is None:
            best = current
        elif (
            current["balanced_accuracy"] > best["balanced_accuracy"]
            or (
                current["balanced_accuracy"] == best["balanced_accuracy"]
                and current["accuracy"] > best["accuracy"]
            )
        ):
            best = current

    return best


def tree_stump_threshold(df, indicator, response):
    sub = df[[indicator, response]].dropna().copy()

    if sub.shape[0] < 2 * MIN_LEAF_SIZE:
        return None

    x = sub[[indicator]].to_numpy()
    y = sub[response].to_numpy()

    model = DecisionTreeRegressor(
        max_depth=1,
        min_samples_leaf=MIN_LEAF_SIZE,
        random_state=42,
    )
    model.fit(x, y)

    threshold = model.tree_.threshold[0]

    if threshold == -2:
        return None

    left_mask = sub[indicator] <= threshold
    right_mask = sub[indicator] > threshold

    return {
        "indicator": indicator,
        "threshold_tree": float(threshold),
        "tree_left_n": int(left_mask.sum()),
        "tree_right_n": int(right_mask.sum()),
        "tree_left_mean_yield": float(sub.loc[left_mask, response].mean()),
        "tree_right_mean_yield": float(sub.loc[right_mask, response].mean()),
        "tree_yield_difference": float(
            sub.loc[right_mask, response].mean()
            - sub.loc[left_mask, response].mean()
        ),
        "tree_r2_in_sample": float(model.score(x, y)),
    }


def high_yield_range(df, indicator, response, yield_cutoff):
    sub = df[[indicator, response]].dropna().copy()
    high = sub[sub[response] >= yield_cutoff]

    if high.empty:
        return None

    q = high[indicator].quantile([0.10, 0.25, 0.50, 0.75, 0.90])

    return {
        "indicator": indicator,
        "n_high_yield": int(high.shape[0]),
        "high_yield_p10": float(q.loc[0.10]),
        "high_yield_p25": float(q.loc[0.25]),
        "high_yield_median": float(q.loc[0.50]),
        "high_yield_p75": float(q.loc[0.75]),
        "high_yield_p90": float(q.loc[0.90]),
    }


def spearman_summary(df, indicator, response):
    sub = df[[indicator, response]].dropna()

    if sub.shape[0] < 3:
        return {
            "spearman_rho": None,
            "spearman_p": None,
        }

    result = stats.spearmanr(sub[indicator], sub[response])

    return {
        "spearman_rho": float(result.statistic),
        "spearman_p": float(result.pvalue),
    }


def analyze_subset(df, subset_name):
    rows = []

    for indicator in INDICATORS:
        cate = cate_nelson_fixed_y(df, indicator, RESPONSE, MAIN_YIELD_CUTOFF)
        tree = tree_stump_threshold(df, indicator, RESPONSE)
        high_range = high_yield_range(df, indicator, RESPONSE, MAIN_YIELD_CUTOFF)
        spearman = spearman_summary(df, indicator, RESPONSE)

        row = {
            "subset": subset_name,
            "indicator": indicator,
            "n_subset": int(df.shape[0]),
            "yield_cutoff": MAIN_YIELD_CUTOFF,
        }

        if cate is not None:
            row.update(cate)
        if tree is not None:
            row.update(tree)
        if high_range is not None:
            row.update(high_range)
        row.update(spearman)

        rows.append(row)

    return pd.DataFrame(rows)


def main():
    if not INPUT_PATH.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_PATH}")

    df = pd.read_csv(INPUT_PATH)

    experimental_mask = detect_experimental_mask(df)

    print("Experimental detection:")
    print(f"- Total rows: {df.shape[0]}")
    print(f"- Experimental-like rows: {int(experimental_mask.sum())}")
    print(f"- Non-experimental rows: {int((~experimental_mask).sum())}")

    all_data = df.copy()
    without_experimental = df.loc[~experimental_mask].copy()

    results = pd.concat(
        [
            analyze_subset(all_data, "all_data"),
            analyze_subset(without_experimental, "without_experimental"),
        ],
        ignore_index=True,
    )

    out_path = TABLE_DIR / "experimental_group_threshold_sensitivity_private.csv"
    results.to_csv(out_path, index=False)

    print(f"\nOutput file: {out_path}")

    print("\nMain 70% threshold sensitivity:")
    display_cols = [
        "subset",
        "indicator",
        "n_subset",
        "n_high_yield",
        "spearman_rho",
        "threshold_cate_nelson",
        "balanced_accuracy",
        "threshold_tree",
        "tree_yield_difference",
        "high_yield_p25",
        "high_yield_median",
        "high_yield_p75",
    ]

    key = results[
        results["indicator"].isin(
            ["GMea", "Beta_glic", "Arilsulf", "GMea_per_Clay", "Beta_por_Argila"]
        )
    ].copy()

    print(key[display_cols].to_string(index=False))


if __name__ == "__main__":
    main()
