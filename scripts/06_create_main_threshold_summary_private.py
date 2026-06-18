from pathlib import Path

import pandas as pd


TABLE_DIR = Path("tables/private")
OUT_PATH = TABLE_DIR / "main_enzyme_threshold_summary_private.csv"

CORR_PATH = TABLE_DIR / "enzyme_response_correlations_private.csv"
SENS_PATH = TABLE_DIR / "experimental_group_threshold_sensitivity_private.csv"

MAIN_INDICATORS = [
    "GMea",
    "Beta_glic",
    "Arilsulf",
    "Ativ_Enzimat_Total",
    "GMea_per_Clay",
    "Beta_por_Argila",
]

MAIN_YIELD_CUTOFF = 70


def main():
    corr = pd.read_csv(CORR_PATH)
    sens = pd.read_csv(SENS_PATH)

    corr_main = corr[
        (corr["response"] == "Prod_rel_pct")
        & (corr["indicator"].isin(MAIN_INDICATORS))
        & (corr["indicator_group"] == "enzyme_candidate")
    ][
        [
            "indicator",
            "n",
            "spearman_rho",
            "spearman_p",
            "pearson_r",
            "pearson_p",
        ]
    ].copy()

    sens_main = sens[
        (sens["subset"] == "all_data")
        & (sens["yield_cutoff"] == MAIN_YIELD_CUTOFF)
        & (sens["indicator"].isin(MAIN_INDICATORS))
    ].copy()

    sens_cols = [
        "indicator",
        "n_high_yield",
        "threshold_cate_nelson",
        "balanced_accuracy",
        "threshold_tree",
        "tree_left_mean_yield",
        "tree_right_mean_yield",
        "tree_yield_difference",
        "high_yield_p25",
        "high_yield_median",
        "high_yield_p75",
    ]

    summary = corr_main.merge(
        sens_main[sens_cols],
        on="indicator",
        how="left",
    )

    interpretation = {
        "GMea": "Main integrated biological indicator",
        "Beta_glic": "Sensitive minimum-level indicator",
        "Arilsulf": "Complementary sulfur-cycling enzyme indicator",
        "Ativ_Enzimat_Total": "Alternative aggregate enzyme indicator",
        "GMea_per_Clay": "Clay-normalized sensitivity indicator",
        "Beta_por_Argila": "Clay-normalized beta-glucosidase sensitivity indicator",
    }

    summary["proposed_role"] = summary["indicator"].map(interpretation)

    summary = summary[
        [
            "indicator",
            "proposed_role",
            "n",
            "n_high_yield",
            "spearman_rho",
            "spearman_p",
            "pearson_r",
            "pearson_p",
            "threshold_tree",
            "threshold_cate_nelson",
            "balanced_accuracy",
            "tree_left_mean_yield",
            "tree_right_mean_yield",
            "tree_yield_difference",
            "high_yield_p25",
            "high_yield_median",
            "high_yield_p75",
        ]
    ].sort_values("spearman_rho", ascending=False)

    summary.to_csv(OUT_PATH, index=False)

    print("Main private enzyme threshold summary created.")
    print(f"Output file: {OUT_PATH}")

    print("\nMain enzyme threshold summary:")
    display_cols = [
        "indicator",
        "proposed_role",
        "spearman_rho",
        "threshold_tree",
        "threshold_cate_nelson",
        "balanced_accuracy",
        "high_yield_p25",
        "high_yield_median",
        "high_yield_p75",
    ]
    print(summary[display_cols].to_string(index=False))


if __name__ == "__main__":
    main()
