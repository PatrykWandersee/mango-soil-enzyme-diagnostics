from pathlib import Path

import numpy as np
import pandas as pd


TABLE_DIR = Path("tables/private")
TABLE_DIR.mkdir(parents=True, exist_ok=True)

MAIN_PATH = TABLE_DIR / "main_enzyme_threshold_summary_private.csv"
BOOT_PATH = TABLE_DIR / "enzyme_threshold_bootstrap_summary_private.csv"
MODEL_PATH = TABLE_DIR / "enzyme_response_model_comparison_private.csv"
STRUCTURE_MODEL_PATH = TABLE_DIR / "enzyme_orchard_structure_adjusted_models_private.csv"
FARM_MODEL_PATH = TABLE_DIR / "enzyme_farm_effect_models_private.csv"

OUT_SYNTHESIS_CSV = TABLE_DIR / "article_enzyme_diagnostic_synthesis_private.csv"
OUT_SYNTHESIS_MD = TABLE_DIR / "article_enzyme_diagnostic_synthesis_private.md"

OUT_MODEL_CSV = TABLE_DIR / "article_enzyme_model_support_private.csv"
OUT_MODEL_MD = TABLE_DIR / "article_enzyme_model_support_private.md"

INDICATORS = {
    "GMea": {
        "label": "GMea",
        "role": "Main integrated biological indicator",
        "interpretation": (
            "Best overall integrated biological indicator. Strong positive association with yield, "
            "clear diagnostic threshold, and support for non-linear/plateau-like behavior. "
            "Signal is partly integrated with soil organic matter, soil structure, farm context, "
            "orchard age, and planting density."
        ),
    },
    "Beta_glic": {
        "label": "β-glucosidase",
        "role": "Most robust individual enzyme indicator",
        "interpretation": (
            "Most consistent individual enzyme. Maintains useful signal after adjustment for "
            "orchard structure and soil context. Suitable as a minimum biological-level indicator."
        ),
    },
    "Arilsulf": {
        "label": "Arylsulfatase",
        "role": "Complementary sulfur-cycling enzyme indicator",
        "interpretation": (
            "Useful complementary component, but less independent as an isolated diagnostic indicator. "
            "Strongly associated with organic matter, exchangeable bases, soil structure, and farm context."
        ),
    },
}


def read_csv_or_empty(path):
    if path.exists():
        return pd.read_csv(path)
    print(f"Warning: missing file: {path}")
    return pd.DataFrame()


def get_value(df, indicator, column, default=np.nan):
    if df.empty or column not in df.columns or "indicator" not in df.columns:
        return default

    subset = df[df["indicator"] == indicator]

    if subset.empty:
        return default

    value = subset[column].iloc[0]
    return value


def get_model_value(df, indicator, model_name, column, default=np.nan):
    if df.empty or column not in df.columns:
        return default

    if "enzyme" not in df.columns or "model_name" not in df.columns:
        return default

    subset = df[(df["enzyme"] == indicator) & (df["model_name"] == model_name)]

    if subset.empty:
        return default

    return subset[column].iloc[0]


def get_boot_metric(boot_df, indicator, metric):
    if boot_df.empty:
        return {}

    required = {"indicator", "metric"}
    if not required.issubset(set(boot_df.columns)):
        return {}

    subset = boot_df[
        (boot_df["indicator"] == indicator)
        & (boot_df["metric"] == metric)
    ]

    if subset.empty:
        return {}

    return subset.iloc[0].to_dict()


def fmt_number(value, digits=3):
    if pd.isna(value):
        return ""

    return f"{float(value):.{digits}f}"


def fmt_p(value):
    if pd.isna(value):
        return ""

    value = float(value)

    if value < 0.001:
        return "<0.001"

    return f"{value:.3f}"


def fmt_range(low, high, digits=1):
    if pd.isna(low) or pd.isna(high):
        return ""

    return f"{float(low):.{digits}f}–{float(high):.{digits}f}"


def fmt_boot_ci(metric_row, digits=1):
    if not metric_row:
        return ""

    median = metric_row.get("bootstrap_median", np.nan)
    low = metric_row.get("ci_2_5", np.nan)
    high = metric_row.get("ci_97_5", np.nan)

    if pd.isna(median) or pd.isna(low) or pd.isna(high):
        return ""

    return f"{float(median):.{digits}f} ({float(low):.{digits}f}–{float(high):.{digits}f})"


def fmt_coef_p(coef, p, digits=2):
    if pd.isna(coef) or pd.isna(p):
        return ""

    p_text = fmt_p(p)
    if p_text.startswith("<"):
        return f"{float(coef):.{digits}f}; p{p_text}"
    return f"{float(coef):.{digits}f}; p={p_text}"


def competitive_models(model_df, indicator, delta_cutoff=2.0):
    if model_df.empty:
        return ""

    if "indicator" not in model_df.columns:
        return ""

    subset = model_df[
        (model_df["indicator"] == indicator)
        & (model_df["delta_aicc"] <= delta_cutoff)
    ].copy()

    if subset.empty:
        return ""

    subset = subset.sort_values("delta_aicc")

    parts = []
    for _, row in subset.iterrows():
        parts.append(f"{row['model']} (ΔAICc={row['delta_aicc']:.2f})")

    return "; ".join(parts)


def best_model(model_df, indicator):
    if model_df.empty or "indicator" not in model_df.columns:
        return "", np.nan, np.nan, np.nan

    subset = model_df[model_df["indicator"] == indicator].copy()

    if subset.empty:
        return "", np.nan, np.nan, np.nan

    subset = subset.sort_values("delta_aicc")
    row = subset.iloc[0]

    return (
        row.get("model", ""),
        row.get("r2", np.nan),
        row.get("rmse", np.nan),
        row.get("breakpoint_or_vertex", np.nan),
    )


def to_markdown_table(df, path):
    text_cols = [str(col) for col in df.columns]

    lines = []
    lines.append("| " + " | ".join(text_cols) + " |")
    lines.append("| " + " | ".join(["---"] * len(text_cols)) + " |")

    for _, row in df.iterrows():
        values = []
        for col in df.columns:
            value = row[col]
            if pd.isna(value):
                values.append("")
            else:
                values.append(str(value).replace("\n", " "))
        lines.append("| " + " | ".join(values) + " |")

    path.write_text("\n".join(lines) + "\n")


def main():
    main_df = read_csv_or_empty(MAIN_PATH)
    boot_df = read_csv_or_empty(BOOT_PATH)
    model_df = read_csv_or_empty(MODEL_PATH)
    structure_df = read_csv_or_empty(STRUCTURE_MODEL_PATH)
    farm_df = read_csv_or_empty(FARM_MODEL_PATH)

    synthesis_rows = []
    model_rows = []

    for indicator, meta in INDICATORS.items():
        threshold_cate = get_value(main_df, indicator, "threshold_cate_nelson")
        threshold_tree = get_value(main_df, indicator, "threshold_tree")
        spearman_rho = get_value(main_df, indicator, "spearman_rho")
        spearman_p = get_value(main_df, indicator, "spearman_p")
        balanced_accuracy = get_value(main_df, indicator, "balanced_accuracy")

        high_p25 = get_value(main_df, indicator, "high_yield_p25")
        high_median = get_value(main_df, indicator, "high_yield_median")
        high_p75 = get_value(main_df, indicator, "high_yield_p75")

        cate_boot = get_boot_metric(
            boot_df,
            indicator,
            "cate_nelson_threshold_70",
        )
        tree_boot = get_boot_metric(
            boot_df,
            indicator,
            "tree_threshold",
        )
        plateau_boot = get_boot_metric(
            boot_df,
            indicator,
            "linear_plateau_breakpoint",
        )
        iqr_p25_boot = get_boot_metric(
            boot_df,
            indicator,
            "high_yield_p25",
        )
        iqr_p75_boot = get_boot_metric(
            boot_df,
            indicator,
            "high_yield_p75",
        )

        best_model_name, best_model_r2, best_model_rmse, best_break = best_model(
            model_df,
            indicator,
        )

        structure_p = get_model_value(
            structure_df,
            indicator,
            "orchard_structure",
            "enzyme_p_hc3",
        )
        structure_coef = get_model_value(
            structure_df,
            indicator,
            "orchard_structure",
            "enzyme_coef",
        )

        soil_core_p = get_model_value(
            structure_df,
            indicator,
            "soil_core",
            "enzyme_p_hc3",
        )
        soil_core_coef = get_model_value(
            structure_df,
            indicator,
            "soil_core",
            "enzyme_coef",
        )

        farm_fixed_p = get_model_value(
            farm_df,
            indicator,
            "farm_fixed",
            "enzyme_p_hc3",
        )
        farm_fixed_coef = get_model_value(
            farm_df,
            indicator,
            "farm_fixed",
            "enzyme_coef",
        )

        farm_structure_p = get_model_value(
            farm_df,
            indicator,
            "farm_fixed_plus_structure",
            "enzyme_p_hc3",
        )
        farm_structure_coef = get_model_value(
            farm_df,
            indicator,
            "farm_fixed_plus_structure",
            "enzyme_coef",
        )

        farm_all_p = get_model_value(
            farm_df,
            indicator,
            "farm_fixed_plus_structure_and_key_soil",
            "enzyme_p_hc3",
        )
        farm_all_coef = get_model_value(
            farm_df,
            indicator,
            "farm_fixed_plus_structure_and_key_soil",
            "enzyme_coef",
        )

        synthesis_rows.append(
            {
                "indicator": meta["label"],
                "role": meta["role"],
                "Spearman rho": fmt_number(spearman_rho, 3),
                "Spearman p": fmt_p(spearman_p),
                "Cate-Nelson threshold": fmt_number(threshold_cate, 1),
                "Cate-Nelson bootstrap median (95% interval)": fmt_boot_ci(cate_boot, 1),
                "Tree threshold": fmt_number(threshold_tree, 1),
                "Tree bootstrap median (95% interval)": fmt_boot_ci(tree_boot, 1),
                "High-yield IQR": fmt_range(high_p25, high_p75, 1),
                "High-yield IQR bootstrap support": (
                    f"P25 {fmt_boot_ci(iqr_p25_boot, 1)}; "
                    f"P75 {fmt_boot_ci(iqr_p75_boot, 1)}"
                ),
                "Balanced accuracy": fmt_number(balanced_accuracy, 3),
                "Interpretation": meta["interpretation"],
            }
        )

        model_rows.append(
            {
                "indicator": meta["label"],
                "best response model": best_model_name,
                "best model R2": fmt_number(best_model_r2, 3),
                "best model RMSE": fmt_number(best_model_rmse, 2),
                "best model breakpoint or vertex": fmt_number(best_break, 1),
                "competitive models (delta AICc <= 2)": competitive_models(model_df, indicator),
                "linear-plateau bootstrap median (95% interval)": fmt_boot_ci(plateau_boot, 1),
                "orchard-structure adjusted enzyme effect": fmt_coef_p(structure_coef, structure_p),
                "soil-core adjusted enzyme effect": fmt_coef_p(soil_core_coef, soil_core_p),
                "farm fixed-effect enzyme effect": fmt_coef_p(farm_fixed_coef, farm_fixed_p),
                "farm + structure enzyme effect": fmt_coef_p(farm_structure_coef, farm_structure_p),
                "farm + structure + key soil enzyme effect": fmt_coef_p(farm_all_coef, farm_all_p),
            }
        )

    synthesis_df = pd.DataFrame(synthesis_rows)
    model_support_df = pd.DataFrame(model_rows)

    synthesis_df.to_csv(OUT_SYNTHESIS_CSV, index=False)
    model_support_df.to_csv(OUT_MODEL_CSV, index=False)

    to_markdown_table(synthesis_df, OUT_SYNTHESIS_MD)
    to_markdown_table(model_support_df, OUT_MODEL_MD)

    print("Article-level enzyme synthesis created.")
    print(f"Diagnostic synthesis CSV: {OUT_SYNTHESIS_CSV}")
    print(f"Diagnostic synthesis MD:  {OUT_SYNTHESIS_MD}")
    print(f"Model support CSV:        {OUT_MODEL_CSV}")
    print(f"Model support MD:         {OUT_MODEL_MD}")

    print("\nDiagnostic synthesis:")
    print(synthesis_df.to_string(index=False))

    print("\nModel support synthesis:")
    print(model_support_df.to_string(index=False))


if __name__ == "__main__":
    main()
