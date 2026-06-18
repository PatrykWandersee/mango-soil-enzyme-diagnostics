from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm


DATA_PATH = Path("data/processed/private/enzyme_diagnostics_private.csv")
THRESHOLD_PATH = Path("tables/private/main_enzyme_threshold_summary_private.csv")

OUT_DIR = Path("tables/private")
OUT_DIR.mkdir(parents=True, exist_ok=True)

CORR_OUT = OUT_DIR / "enzyme_soil_context_correlations_private.csv"
ADJUSTED_OUT = OUT_DIR / "enzyme_context_adjusted_models_private.csv"
GROUP_OUT = OUT_DIR / "enzyme_threshold_context_group_summary_private.csv"

RESPONSE = "Prod_rel_pct"
YIELD_CUTOFF = 70.0

MAIN_INDICATORS = [
    "GMea",
    "Beta_glic",
    "Arilsulf",
]

CONTEXT_CANDIDATES = [
    "Argila_g_kg",
    "Areia_g_kg",
    "MO_g_dm3",
    "C_org_g_d3",
    "pH",
    "CE_dS_m",
    "PST",
    "Ds_g_cm3",
    "Floculacao_pct",
    "Dispercao_pct",
    "PM1_mg_dm3",
    "SB_cmolc_Kg",
    "Ca_Troc_cmolc_Kg",
    "Mg_Troc_cmolc_Kg",
    "K_Troc_cmolc_Kg",
    "Na_Troc_cmolc_Kg",
    "RAS",
    "CROSS",
    "Idade_anos",
    "Densidade_pl_ha",
]


def spearman_pair(df, x_col, y_col):
    sub = df[[x_col, y_col]].dropna()

    if sub.shape[0] < 4:
        return {
            "n": sub.shape[0],
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
        }

    result = stats.spearmanr(sub[x_col], sub[y_col])

    return {
        "n": sub.shape[0],
        "spearman_rho": float(result.statistic),
        "spearman_p": float(result.pvalue),
    }


def zscore(series):
    values = series.astype(float)
    sd = values.std(ddof=0)

    if sd == 0 or np.isnan(sd):
        return values * np.nan

    return (values - values.mean()) / sd


def fit_adjusted_model(df, enzyme_col, context_cols):
    cols = [RESPONSE, enzyme_col] + context_cols
    sub = df[cols].dropna().copy()

    if sub.shape[0] < len(context_cols) + 8:
        return None

    y = sub[RESPONSE].astype(float)

    x_data = pd.DataFrame(index=sub.index)
    x_data[f"{enzyme_col}_z"] = zscore(sub[enzyme_col])

    for col in context_cols:
        x_data[f"{col}_z"] = zscore(sub[col])

    x_data = sm.add_constant(x_data, has_constant="add")

    try:
        model = sm.OLS(y, x_data).fit()
    except Exception:
        return None

    row = {
        "enzyme": enzyme_col,
        "adjustment": "+".join(context_cols) if context_cols else "none",
        "n": int(sub.shape[0]),
        "r2": float(model.rsquared),
        "adj_r2": float(model.rsquared_adj),
        "aic": float(model.aic),
        "enzyme_coef": float(model.params.get(f"{enzyme_col}_z", np.nan)),
        "enzyme_p": float(model.pvalues.get(f"{enzyme_col}_z", np.nan)),
    }

    for col in context_cols:
        row[f"{col}_coef"] = float(model.params.get(f"{col}_z", np.nan))
        row[f"{col}_p"] = float(model.pvalues.get(f"{col}_z", np.nan))

    return row


def summarize_threshold_groups(df, thresholds, enzyme_col, context_cols):
    threshold = thresholds.get(enzyme_col, np.nan)

    if not np.isfinite(threshold):
        return []

    sub = df[[RESPONSE, enzyme_col] + context_cols].dropna(subset=[RESPONSE, enzyme_col]).copy()
    sub["diagnostic_group"] = np.where(
        sub[enzyme_col] >= threshold,
        "at_or_above_threshold",
        "below_threshold",
    )

    rows = []

    for group_name, group_df in sub.groupby("diagnostic_group"):
        row = {
            "enzyme": enzyme_col,
            "threshold_cate_nelson": threshold,
            "group": group_name,
            "n": int(group_df.shape[0]),
            "yield_mean": float(group_df[RESPONSE].mean()),
            "yield_median": float(group_df[RESPONSE].median()),
            "yield_p25": float(group_df[RESPONSE].quantile(0.25)),
            "yield_p75": float(group_df[RESPONSE].quantile(0.75)),
        }

        for col in context_cols:
            valid = group_df[col].dropna()

            row[f"{col}_median"] = float(valid.median()) if len(valid) else np.nan
            row[f"{col}_p25"] = float(valid.quantile(0.25)) if len(valid) else np.nan
            row[f"{col}_p75"] = float(valid.quantile(0.75)) if len(valid) else np.nan

        rows.append(row)

    return rows


def main():
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Input data not found: {DATA_PATH}")

    df = pd.read_csv(DATA_PATH)

    indicators = [col for col in MAIN_INDICATORS if col in df.columns]
    context_cols = [col for col in CONTEXT_CANDIDATES if col in df.columns]

    missing_indicators = [col for col in MAIN_INDICATORS if col not in df.columns]
    if missing_indicators:
        print("Missing enzyme indicators:", missing_indicators)

    if RESPONSE not in df.columns:
        raise ValueError(f"Response column not found: {RESPONSE}")

    print("Available enzyme indicators:", indicators)
    print("Available context variables:", context_cols)

    # 1) Spearman correlations: enzymes vs soil context; soil context vs yield.
    corr_rows = []

    for enzyme in indicators:
        for context in context_cols:
            result = spearman_pair(df, enzyme, context)
            corr_rows.append(
                {
                    "relationship": "enzyme_vs_context",
                    "x": enzyme,
                    "y": context,
                    **result,
                }
            )

    for context in context_cols:
        result = spearman_pair(df, context, RESPONSE)
        corr_rows.append(
            {
                "relationship": "context_vs_yield",
                "x": context,
                "y": RESPONSE,
                **result,
            }
        )

    corr_df = pd.DataFrame(corr_rows)
    corr_df.to_csv(CORR_OUT, index=False)

    # 2) One-context adjusted models: yield ~ enzyme + context.
    adjusted_rows = []

    for enzyme in indicators:
        model = fit_adjusted_model(df, enzyme, [])
        if model is not None:
            adjusted_rows.append(model)

        for context in context_cols:
            model = fit_adjusted_model(df, enzyme, [context])
            if model is not None:
                adjusted_rows.append(model)

        core_context = [
            col
            for col in [
                "Argila_g_kg",
                "MO_g_dm3",
                "pH",
                "CE_dS_m",
                "PST",
                "Ds_g_cm3",
            ]
            if col in context_cols
        ]

        model = fit_adjusted_model(df, enzyme, core_context)
        if model is not None:
            adjusted_rows.append(model)

    adjusted_df = pd.DataFrame(adjusted_rows)
    adjusted_df.to_csv(ADJUSTED_OUT, index=False)

    # 3) Context summary below vs above Cate-Nelson thresholds.
    thresholds = {}

    if THRESHOLD_PATH.exists():
        th = pd.read_csv(THRESHOLD_PATH)
        if "indicator" in th.columns and "threshold_cate_nelson" in th.columns:
            thresholds = dict(zip(th["indicator"], th["threshold_cate_nelson"]))

    group_rows = []
    for enzyme in indicators:
        group_rows.extend(
            summarize_threshold_groups(
                df=df,
                thresholds=thresholds,
                enzyme_col=enzyme,
                context_cols=context_cols,
            )
        )

    group_df = pd.DataFrame(group_rows)
    group_df.to_csv(GROUP_OUT, index=False)

    print("\nSaved outputs:")
    print(f"- {CORR_OUT}")
    print(f"- {ADJUSTED_OUT}")
    print(f"- {GROUP_OUT}")

    print("\nStrongest enzyme-context Spearman correlations:")
    enzyme_context = corr_df[corr_df["relationship"] == "enzyme_vs_context"].copy()
    enzyme_context["abs_rho"] = enzyme_context["spearman_rho"].abs()
    print(
        enzyme_context
        .sort_values("abs_rho", ascending=False)
        .head(20)[["x", "y", "n", "spearman_rho", "spearman_p"]]
        .to_string(index=False)
    )

    print("\nStrongest context-yield Spearman correlations:")
    context_yield = corr_df[corr_df["relationship"] == "context_vs_yield"].copy()
    context_yield["abs_rho"] = context_yield["spearman_rho"].abs()
    print(
        context_yield
        .sort_values("abs_rho", ascending=False)
        .head(15)[["x", "y", "n", "spearman_rho", "spearman_p"]]
        .to_string(index=False)
    )

    print("\nEnzyme-only and core-context adjusted models:")
    if not adjusted_df.empty:
        view = adjusted_df[
            adjusted_df["adjustment"].isin(
                [
                    "none",
                    "Argila_g_kg+MO_g_dm3+pH+CE_dS_m+PST+Ds_g_cm3",
                ]
            )
        ].copy()

        display_cols = [
            "enzyme",
            "adjustment",
            "n",
            "r2",
            "adj_r2",
            "enzyme_coef",
            "enzyme_p",
        ]

        print(view[display_cols].to_string(index=False))
    else:
        print("No adjusted models were fitted.")


if __name__ == "__main__":
    main()
