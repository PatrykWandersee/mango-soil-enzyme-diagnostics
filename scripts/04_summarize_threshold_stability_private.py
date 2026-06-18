from pathlib import Path

import pandas as pd


TABLE_DIR = Path("tables/private")

CATE_PATH = TABLE_DIR / "enzyme_cate_nelson_thresholds_private.csv"
TREE_PATH = TABLE_DIR / "enzyme_tree_stump_thresholds_private.csv"
RANGE_PATH = TABLE_DIR / "enzyme_high_yield_reference_ranges_private.csv"

OUT_PATH = TABLE_DIR / "enzyme_threshold_stability_summary_private.csv"

KEY_INDICATORS = [
    "GMea",
    "Beta_glic",
    "Arilsulf",
    "Ativ_Enzimat_Total",
    "GMea_per_Clay",
    "log_GMea_per_Clay",
    "Beta_por_Argila",
    "Aril_por_Argila",
]


def main():
    cate = pd.read_csv(CATE_PATH)
    tree = pd.read_csv(TREE_PATH)
    ranges = pd.read_csv(RANGE_PATH)

    cate_key = cate[cate["indicator"].isin(KEY_INDICATORS)].copy()
    ranges_key = ranges[ranges["indicator"].isin(KEY_INDICATORS)].copy()
    tree_key = tree[tree["indicator"].isin(KEY_INDICATORS)].copy()

    merged = cate_key.merge(
        ranges_key[
            [
                "indicator",
                "yield_cutoff",
                "n_high_yield",
                "p10",
                "p25",
                "median",
                "p75",
                "p90",
            ]
        ],
        on=["indicator", "yield_cutoff"],
        how="left",
    )

    tree_key = tree_key[
        [
            "indicator",
            "threshold",
            "left_n",
            "right_n",
            "left_mean_yield",
            "right_mean_yield",
            "yield_difference_right_minus_left",
            "r2_in_sample",
        ]
    ].rename(
        columns={
            "threshold": "tree_minimum_threshold",
            "left_n": "tree_left_n",
            "right_n": "tree_right_n",
            "left_mean_yield": "tree_left_mean_yield",
            "right_mean_yield": "tree_right_mean_yield",
            "yield_difference_right_minus_left": "tree_yield_difference",
            "r2_in_sample": "tree_r2_in_sample",
        }
    )

    merged = merged.merge(tree_key, on="indicator", how="left")

    cols = [
        "indicator",
        "yield_cutoff",
        "n",
        "n_high_yield",
        "threshold",
        "balanced_accuracy",
        "accuracy",
        "low_indicator_low_yield",
        "low_indicator_high_yield",
        "high_indicator_low_yield",
        "high_indicator_high_yield",
        "spearman_rho",
        "p10",
        "p25",
        "median",
        "p75",
        "p90",
        "tree_minimum_threshold",
        "tree_left_n",
        "tree_right_n",
        "tree_left_mean_yield",
        "tree_right_mean_yield",
        "tree_yield_difference",
        "tree_r2_in_sample",
    ]

    merged = merged[cols].sort_values(["indicator", "yield_cutoff"])
    merged.to_csv(OUT_PATH, index=False)

    print("Threshold stability summary created.")
    print(f"Output file: {OUT_PATH}")

    print("\nKey threshold stability summary:")
    display_cols = [
        "indicator",
        "yield_cutoff",
        "n_high_yield",
        "threshold",
        "balanced_accuracy",
        "p25",
        "median",
        "tree_minimum_threshold",
        "tree_yield_difference",
    ]

    for indicator in KEY_INDICATORS:
        sub = merged[merged["indicator"] == indicator]
        if sub.empty:
            continue
        print(f"\n{indicator}:")
        print(sub[display_cols].to_string(index=False))


if __name__ == "__main__":
    main()
