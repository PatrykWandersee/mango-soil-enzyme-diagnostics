from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats


INPUT_PATH = Path("data/processed/private/enzyme_diagnostics_private.csv")

TABLE_DIR = Path("tables/private")
FIGURE_DIR = Path("figures/private")

TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)


RESPONSE_COLUMNS = [
    "Prod_rel_pct",
    "Prod_rel_ha_pct",
]

ENZYME_CANDIDATES = [
    "Beta_glic",
    "Arilsulf",
    "GMea",
    "Ativ_Enzimat_Total",
    "Rel_Beta_Aril",
    "Beta_por_Argila",
    "Aril_por_Argila",
    "GMea_por_Argila",
    "Ativ_Enzim_por_MO",
    "GMea_per_Clay",
    "log_GMea_per_Clay",
    "qBeta",
    "qAril",
    "qGMea",
]

CONTEXT_COLUMNS = [
    "MO_g_dm3",
    "C_org_g_d3",
    "Argila_g_kg",
    "Floculacao_pct",
    "Ds_g_cm3",
    "pH",
    "CE_dS_m",
    "PST",
    "PM1_mg_dm3",
    "SB_cmolc_Kg",
]


def correlation_pair(df, x_col, y_col):
    sub = df[[x_col, y_col]].dropna()

    if sub.shape[0] < 3:
        return {
            "n": sub.shape[0],
            "spearman_rho": None,
            "spearman_p": None,
            "pearson_r": None,
            "pearson_p": None,
        }

    spearman = stats.spearmanr(sub[x_col], sub[y_col])
    pearson = stats.pearsonr(sub[x_col], sub[y_col])

    return {
        "n": sub.shape[0],
        "spearman_rho": spearman.statistic,
        "spearman_p": spearman.pvalue,
        "pearson_r": pearson.statistic,
        "pearson_p": pearson.pvalue,
    }


def save_scatter(df, x_col, y_col, output_path):
    sub = df[[x_col, y_col]].dropna()

    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    ax.scatter(sub[x_col], sub[y_col], alpha=0.8)
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)

    if sub.shape[0] >= 3:
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            sub[x_col], sub[y_col]
        )
        x_min = sub[x_col].min()
        x_max = sub[x_col].max()
        x_line = [x_min, x_max]
        y_line = [intercept + slope * x_min, intercept + slope * x_max]
        ax.plot(x_line, y_line)
        ax.text(
            0.04,
            0.96,
            f"r = {r_value:.3f}\np = {p_value:.3f}",
            transform=ax.transAxes,
            va="top",
            ha="left",
        )

    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def main():
    if not INPUT_PATH.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_PATH}")

    df = pd.read_csv(INPUT_PATH)

    required_columns = RESPONSE_COLUMNS + ENZYME_CANDIDATES + CONTEXT_COLUMNS
    missing = [col for col in required_columns if col not in df.columns]
    if missing:
        raise ValueError(
            "Missing required columns:\n" + "\n".join(f"- {col}" for col in missing)
        )

    rows = []

    for response in RESPONSE_COLUMNS:
        for indicator in ENZYME_CANDIDATES:
            result = correlation_pair(df, indicator, response)
            rows.append(
                {
                    "indicator_group": "enzyme_candidate",
                    "indicator": indicator,
                    "response": response,
                    **result,
                }
            )

        for indicator in CONTEXT_COLUMNS:
            result = correlation_pair(df, indicator, response)
            rows.append(
                {
                    "indicator_group": "soil_context",
                    "indicator": indicator,
                    "response": response,
                    **result,
                }
            )

    corr_table = pd.DataFrame(rows)
    corr_table["abs_spearman_rho"] = corr_table["spearman_rho"].abs()
    corr_table = corr_table.sort_values(
        ["response", "indicator_group", "abs_spearman_rho"],
        ascending=[True, True, False],
    )

    corr_path = TABLE_DIR / "enzyme_response_correlations_private.csv"
    corr_table.to_csv(corr_path, index=False)

    enzyme_only = corr_table[
        (corr_table["indicator_group"] == "enzyme_candidate")
        & (corr_table["response"] == "Prod_rel_pct")
    ].copy()
    top_enzyme = enzyme_only.sort_values("abs_spearman_rho", ascending=False).head(8)

    top_path = TABLE_DIR / "top_enzyme_indicators_prod_rel_pct_private.csv"
    top_enzyme.to_csv(top_path, index=False)

    # Save scatter plots for the main raw/derived enzyme indicators.
    main_plot_indicators = [
        "Beta_glic",
        "Arilsulf",
        "GMea",
        "Beta_por_Argila",
        "Aril_por_Argila",
        "GMea_per_Clay",
        "Ativ_Enzim_por_MO",
        "qGMea",
    ]

    for indicator in main_plot_indicators:
        output_path = FIGURE_DIR / f"scatter_{indicator}_vs_Prod_rel_pct_private.png"
        save_scatter(df, indicator, "Prod_rel_pct", output_path)

    print("Private enzyme relationship exploration completed.")
    print(f"Input shape: {df.shape[0]} rows x {df.shape[1]} columns")
    print(f"Correlation table: {corr_path}")
    print(f"Top enzyme table: {top_path}")
    print(f"Scatter plots saved to: {FIGURE_DIR}")

    print("\nTop enzyme indicators for relative yield per plant:")
    display_cols = [
        "indicator",
        "n",
        "spearman_rho",
        "spearman_p",
        "pearson_r",
        "pearson_p",
    ]
    print(top_enzyme[display_cols].to_string(index=False))


if __name__ == "__main__":
    main()
