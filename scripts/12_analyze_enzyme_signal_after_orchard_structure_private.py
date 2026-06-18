from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm


DATA_PATH = Path("data/processed/private/enzyme_diagnostics_private.csv")

OUT_DIR = Path("tables/private")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_MODELS = OUT_DIR / "enzyme_orchard_structure_adjusted_models_private.csv"
OUT_PARTIAL = OUT_DIR / "enzyme_partial_correlations_private.csv"

RESPONSE = "Prod_rel_pct"

ENZYMES = [
    "GMea",
    "Beta_glic",
    "Arilsulf",
]

MODEL_SETS = {
    "enzyme_only": [],
    "orchard_structure": ["Idade_anos", "Densidade_pl_ha"],
    "soil_core": ["Argila_g_kg", "MO_g_dm3", "pH", "CE_dS_m", "PST", "Ds_g_cm3"],
    "structure_plus_key_soil": ["Idade_anos", "Densidade_pl_ha", "MO_g_dm3", "PST", "Ds_g_cm3"],
    "structure_plus_fertility": ["Idade_anos", "Densidade_pl_ha", "MO_g_dm3", "SB_cmolc_Kg", "PST"],
}


def zscore(series):
    values = series.astype(float)
    sd = values.std(ddof=0)

    if sd == 0 or np.isnan(sd):
        return values * np.nan

    return (values - values.mean()) / sd


def fit_model(df, enzyme, covariates):
    required = [RESPONSE, enzyme] + covariates
    required = [col for col in required if col in df.columns]

    sub = df[required].dropna().copy()

    if sub.shape[0] < len(covariates) + 8:
        return None

    y = sub[RESPONSE].astype(float)

    x = pd.DataFrame(index=sub.index)
    x[f"{enzyme}_z"] = zscore(sub[enzyme])

    used_covariates = []
    for cov in covariates:
        if cov in sub.columns:
            x[f"{cov}_z"] = zscore(sub[cov])
            used_covariates.append(cov)

    x = sm.add_constant(x, has_constant="add")

    try:
        model = sm.OLS(y, x).fit(cov_type="HC3")
    except Exception:
        return None

    row = {
        "enzyme": enzyme,
        "model_set": "+".join(used_covariates) if used_covariates else "enzyme_only",
        "n": int(sub.shape[0]),
        "r2": float(model.rsquared),
        "adj_r2": float(model.rsquared_adj),
        "aic": float(model.aic),
        "enzyme_coef": float(model.params.get(f"{enzyme}_z", np.nan)),
        "enzyme_p_hc3": float(model.pvalues.get(f"{enzyme}_z", np.nan)),
    }

    for cov in used_covariates:
        row[f"{cov}_coef"] = float(model.params.get(f"{cov}_z", np.nan))
        row[f"{cov}_p_hc3"] = float(model.pvalues.get(f"{cov}_z", np.nan))

    return row


def residualize(series, covariate_df):
    data = pd.concat([series, covariate_df], axis=1).dropna()

    y = data.iloc[:, 0].astype(float)

    if covariate_df.shape[1] == 0:
        return pd.Series(y - y.mean(), index=data.index)

    x = pd.DataFrame(index=data.index)
    for col in data.columns[1:]:
        x[col] = zscore(data[col])

    x = sm.add_constant(x, has_constant="add")
    model = sm.OLS(y, x).fit()

    return pd.Series(model.resid, index=data.index)


def partial_spearman(df, enzyme, covariates):
    available_covariates = [cov for cov in covariates if cov in df.columns]
    required = [RESPONSE, enzyme] + available_covariates

    sub = df[required].dropna().copy()

    if sub.shape[0] < len(available_covariates) + 8:
        return None

    cov_df = sub[available_covariates].copy()

    y_resid = residualize(sub[RESPONSE], cov_df)
    x_resid = residualize(sub[enzyme], cov_df)

    aligned = pd.concat(
        [
            y_resid.rename("yield_resid"),
            x_resid.rename("enzyme_resid"),
        ],
        axis=1,
    ).dropna()

    if aligned.shape[0] < 4:
        return None

    result = stats.spearmanr(aligned["enzyme_resid"], aligned["yield_resid"])

    return {
        "enzyme": enzyme,
        "adjustment": "+".join(available_covariates) if available_covariates else "none",
        "n": int(aligned.shape[0]),
        "partial_spearman_rho": float(result.statistic),
        "partial_spearman_p": float(result.pvalue),
    }


def main():
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Input data not found: {DATA_PATH}")

    df = pd.read_csv(DATA_PATH)

    model_rows = []

    for enzyme in ENZYMES:
        if enzyme not in df.columns:
            print(f"Skipping missing enzyme: {enzyme}")
            continue

        for model_name, covariates in MODEL_SETS.items():
            row = fit_model(df, enzyme, covariates)
            if row is not None:
                row["model_name"] = model_name
                model_rows.append(row)

    model_df = pd.DataFrame(model_rows)
    model_df = model_df[
        [
            "enzyme",
            "model_name",
            "model_set",
            "n",
            "r2",
            "adj_r2",
            "aic",
            "enzyme_coef",
            "enzyme_p_hc3",
        ]
        + [
            col for col in model_df.columns
            if col.endswith("_coef") or col.endswith("_p_hc3")
        ]
    ]

    model_df.to_csv(OUT_MODELS, index=False)

    partial_rows = []

    for enzyme in ENZYMES:
        if enzyme not in df.columns:
            continue

        for model_name, covariates in MODEL_SETS.items():
            row = partial_spearman(df, enzyme, covariates)
            if row is not None:
                row["model_name"] = model_name
                partial_rows.append(row)

    partial_df = pd.DataFrame(partial_rows)
    partial_df.to_csv(OUT_PARTIAL, index=False)

    print("Orchard-structure-adjusted enzyme analysis completed.")
    print(f"Model output: {OUT_MODELS}")
    print(f"Partial correlation output: {OUT_PARTIAL}")

    print("\nAdjusted model summary:")
    display_cols = [
        "enzyme",
        "model_name",
        "n",
        "r2",
        "adj_r2",
        "enzyme_coef",
        "enzyme_p_hc3",
    ]
    print(model_df[display_cols].to_string(index=False))

    print("\nPartial Spearman summary:")
    display_cols = [
        "enzyme",
        "model_name",
        "n",
        "partial_spearman_rho",
        "partial_spearman_p",
    ]
    print(partial_df[display_cols].to_string(index=False))


if __name__ == "__main__":
    main()
