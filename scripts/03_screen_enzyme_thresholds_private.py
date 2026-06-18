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

INDICATORS = [
    "GMea",
    "Ativ_Enzimat_Total",
    "Beta_glic",
    "Arilsulf",
    "GMea_por_Argila",
    "GMea_per_Clay",
    "log_GMea_per_Clay",
    "Beta_por_Argila",
    "Aril_por_Argila",
    "Ativ_Enzim_por_MO",
    "qGMea",
]

YIELD_CLASS_CUTOFFS = [70, 80, 90]
MIN_LEAF_SIZE = 8


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

        low_low = int(((x_class == 0) & (y_class == 0)).sum())
        low_high = int(((x_class == 0) & (y_class == 1)).sum())
        high_low = int(((x_class == 1) & (y_class == 0)).sum())
        high_high = int(((x_class == 1) & (y_class == 1)).sum())

        spearman = stats.spearmanr(sub[indicator], sub[response])

        current = {
            "method": "cate_nelson_fixed_y",
            "indicator": indicator,
            "response": response,
            "yield_cutoff": yield_cutoff,
            "n": int(sub.shape[0]),
            "threshold": float(threshold),
            "accuracy": float(accuracy),
            "balanced_accuracy": float(balanced_accuracy),
            "low_indicator_low_yield": low_low,
            "low_indicator_high_yield": low_high,
            "high_indicator_low_yield": high_low,
            "high_indicator_high_yield": high_high,
            "spearman_rho": float(spearman.statistic),
            "spearman_p": float(spearman.pvalue),
        }

        if best is None:
            best = current
        else:
            if (
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
    if sub.shape[0] < (2 * MIN_LEAF_SIZE):
        return None

    x = sub[[indicator]].to_numpy()
    y = sub[response].to_numpy()

    model = DecisionTreeRegressor(
        max_depth=1,
        min_samples_leaf=MIN_LEAF_SIZE,
        random_state=42,
    )
    model.fit(x, y)

    tree = model.tree_
    threshold = tree.threshold[0]

    if threshold == -2:
        return None

    left_mask = sub[indicator] <= threshold
    right_mask = sub[indicator] > threshold

    left_mean = sub.loc[left_mask, response].mean()
    right_mean = sub.loc[right_mask, response].mean()

    return {
        "method": "regression_tree_stump",
        "indicator": indicator,
        "response": response,
        "n": int(sub.shape[0]),
        "threshold": float(threshold),
        "left_n": int(left_mask.sum()),
        "right_n": int(right_mask.sum()),
        "left_mean_yield": float(left_mean),
        "right_mean_yield": float(right_mean),
        "yield_difference_right_minus_left": float(right_mean - left_mean),
        "r2_in_sample": float(model.score(x, y)),
    }


def high_yield_reference_range(df, indicator, response, yield_cutoff):
    sub = df[[indicator, response]].dropna().copy()
    high = sub[sub[response] >= yield_cutoff]

    if high.empty:
        return None

    quantiles = high[indicator].quantile([0.10, 0.25, 0.50, 0.75, 0.90])

    return {
        "indicator": indicator,
        "response": response,
        "yield_cutoff": yield_cutoff,
        "n_total": int(sub.shape[0]),
        "n_high_yield": int(high.shape[0]),
        "p10": float(quantiles.loc[0.10]),
        "p25": float(quantiles.loc[0.25]),
        "median": float(quantiles.loc[0.50]),
        "p75": float(quantiles.loc[0.75]),
        "p90": float(quantiles.loc[0.90]),
    }


def main():
    if not INPUT_PATH.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_PATH}")

    df = pd.read_csv(INPUT_PATH)

    missing = [col for col in [RESPONSE] + INDICATORS if col not in df.columns]
    if missing:
        raise ValueError("Missing columns:\n" + "\n".join(f"- {col}" for col in missing))

    cate_rows = []
    range_rows = []
    tree_rows = []

    for indicator in INDICATORS:
        tree_result = tree_stump_threshold(df, indicator, RESPONSE)
        if tree_result is not None:
            tree_rows.append(tree_result)

        for yield_cutoff in YIELD_CLASS_CUTOFFS:
            cate_result = cate_nelson_fixed_y(df, indicator, RESPONSE, yield_cutoff)
            if cate_result is not None:
                cate_rows.append(cate_result)

            range_result = high_yield_reference_range(
                df, indicator, RESPONSE, yield_cutoff
            )
            if range_result is not None:
                range_rows.append(range_result)

    cate_table = pd.DataFrame(cate_rows)
    tree_table = pd.DataFrame(tree_rows)
    range_table = pd.DataFrame(range_rows)

    cate_path = TABLE_DIR / "enzyme_cate_nelson_thresholds_private.csv"
    tree_path = TABLE_DIR / "enzyme_tree_stump_thresholds_private.csv"
    range_path = TABLE_DIR / "enzyme_high_yield_reference_ranges_private.csv"

    cate_table.to_csv(cate_path, index=False)
    tree_table.to_csv(tree_path, index=False)
    range_table.to_csv(range_path, index=False)

    print("Private enzyme threshold screening completed.")
    print(f"Cate-Nelson table: {cate_path}")
    print(f"Regression tree stump table: {tree_path}")
    print(f"High-yield reference range table: {range_path}")

    print("\nCate-Nelson-style thresholds for yield cutoff = 80%:")
    cols = [
        "indicator",
        "threshold",
        "balanced_accuracy",
        "accuracy",
        "low_indicator_low_yield",
        "low_indicator_high_yield",
        "high_indicator_low_yield",
        "high_indicator_high_yield",
        "spearman_rho",
    ]
    print(
        cate_table[cate_table["yield_cutoff"] == 80]
        .sort_values(["balanced_accuracy", "accuracy"], ascending=False)[cols]
        .head(10)
        .to_string(index=False)
    )

    print("\nRegression tree stump thresholds:")
    tree_cols = [
        "indicator",
        "threshold",
        "left_n",
        "right_n",
        "left_mean_yield",
        "right_mean_yield",
        "yield_difference_right_minus_left",
        "r2_in_sample",
    ]
    print(
        tree_table.sort_values("yield_difference_right_minus_left", ascending=False)[tree_cols]
        .head(10)
        .to_string(index=False)
    )

    print("\nHigh-yield reference ranges for yield cutoff = 80%:")
    range_cols = [
        "indicator",
        "n_high_yield",
        "p10",
        "p25",
        "median",
        "p75",
        "p90",
    ]
    print(
        range_table[range_table["yield_cutoff"] == 80]
        .sort_values("indicator")[range_cols]
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
