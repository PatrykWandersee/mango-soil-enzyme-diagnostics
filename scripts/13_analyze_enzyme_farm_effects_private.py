from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm


DATA_PATH = Path("data/processed/private/enzyme_diagnostics_private.csv")

OUT_DIR = Path("tables/private")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_MODELS = OUT_DIR / "enzyme_farm_effect_models_private.csv"
OUT_PARTIAL = OUT_DIR / "enzyme_farm_partial_correlations_private.csv"
OUT_FARM_SUMMARY = OUT_DIR / "farm_level_enzyme_summary_private.csv"

RESPONSE = "Prod_rel_pct"
FARM = "Fazenda"

ENZYMES = [
    "GMea",
    "Beta_glic",
    "Arilsulf",
]

STRUCTURE_COVARIATES = [
    "Idade_anos",
    "Densidade_pl_ha",
]

KEY_SOIL_COVARIATES = [
    "MO_g_dm3",
    "PST",
    "Ds_g_cm3",
]


def zscore(series):
    values = series.astype(float)
    sd = values.std(ddof=0)

    if sd == 0 or np.isnan(sd):
        return values * np.nan

    return (values - values.mean()) / sd


def build_design_matrix(df, enzyme, include_farm=False, covariates=None):
    if covariates is None:
        covariates = []

    x = pd.DataFrame(index=df.index)
    x[f"{enzyme}_z"] = zscore(df[enzyme])

    used_covariates = []
    for cov in covariates:
        if cov in df.columns:
            x[f"{cov}_z"] = zscore(df[cov])
            used_covariates.append(cov)

    if include_farm:
        farm_dummies = pd.get_dummies(df[FARM].astype(str), prefix="farm", drop_first=True).astype(float)
        x = pd.concat([x, farm_dummies], axis=1)

    x = sm.add_constant(x, has_constant="add")
    x = x.apply(pd.to_numeric, errors="coerce").astype(float)

    return x, used_covariates


def fit_model(df, enzyme, model_name, include_farm=False, covariates=None):
    if covariates is None:
        covariates = []

    required = [RESPONSE, enzyme]
    if include_farm:
        required.append(FARM)
    required += [cov for cov in covariates if cov in df.columns]

    sub = df[required].dropna().copy()

    if sub.shape[0] < 10:
        return None

    y = sub[RESPONSE].astype(float)
    x, used_covariates = build_design_matrix(
        sub,
        enzyme=enzyme,
        include_farm=include_farm,
        covariates=covariates,
    )

    try:
        model = sm.OLS(y, x).fit(cov_type="HC3")
    except Exception as exc:
        return {
            "enzyme": enzyme,
            "model_name": model_name,
            "fit_status": f"failed: {exc}",
        }

    enzyme_term = f"{enzyme}_z"

    return {
        "enzyme": enzyme,
        "model_name": model_name,
        "fit_status": "ok",
        "n": int(sub.shape[0]),
        "n_farms": int(sub[FARM].nunique()) if include_farm else np.nan,
        "r2": float(model.rsquared),
        "adj_r2": float(model.rsquared_adj),
        "aic": float(model.aic),
        "enzyme_coef": float(model.params.get(enzyme_term, np.nan)),
        "enzyme_p_hc3": float(model.pvalues.get(enzyme_term, np.nan)),
        "include_farm": include_farm,
        "covariates": "+".join(used_covariates) if used_covariates else "none",
    }


def residualize_against_covariates(df, target, covariates, include_farm=False):
    required = [target]
    if include_farm:
        required.append(FARM)
    required += [cov for cov in covariates if cov in df.columns]

    sub = df[required].dropna().copy()

    y = sub[target].astype(float)

    x_parts = []

    for cov in covariates:
        if cov in sub.columns:
            x_parts.append(zscore(sub[cov]).rename(f"{cov}_z"))

    if include_farm:
        farm_dummies = pd.get_dummies(sub[FARM].astype(str), prefix="farm", drop_first=True).astype(float)
        x_parts.append(farm_dummies)

    if x_parts:
        x = pd.concat(x_parts, axis=1)
        x = sm.add_constant(x, has_constant="add")
        x = x.apply(pd.to_numeric, errors="coerce").astype(float)
        model = sm.OLS(y, x).fit()
        resid = pd.Series(model.resid, index=sub.index)
    else:
        resid = pd.Series(y - y.mean(), index=sub.index)

    return resid


def partial_spearman(df, enzyme, adjustment_name, include_farm=False, covariates=None):
    if covariates is None:
        covariates = []

    required = [RESPONSE, enzyme]
    if include_farm:
        required.append(FARM)
    required += [cov for cov in covariates if cov in df.columns]

    sub = df[required].dropna().copy()

    if sub.shape[0] < 10:
        return None

    y_resid = residualize_against_covariates(
        sub,
        target=RESPONSE,
        covariates=covariates,
        include_farm=include_farm,
    )

    x_resid = residualize_against_covariates(
        sub,
        target=enzyme,
        covariates=covariates,
        include_farm=include_farm,
    )

    aligned = pd.concat(
        [
            y_resid.rename("yield_resid"),
            x_resid.rename("enzyme_resid"),
        ],
        axis=1,
    ).dropna()

    if aligned.shape[0] < 10:
        return None

    result = stats.spearmanr(aligned["enzyme_resid"], aligned["yield_resid"])

    return {
        "enzyme": enzyme,
        "adjustment_name": adjustment_name,
        "n": int(aligned.shape[0]),
        "partial_spearman_rho": float(result.statistic),
        "partial_spearman_p": float(result.pvalue),
        "include_farm": include_farm,
        "covariates": "+".join([cov for cov in covariates if cov in df.columns]) if covariates else "none",
    }


def farm_level_summary(df):
    available_cols = [FARM, RESPONSE] + [col for col in ENZYMES if col in df.columns]

    extra_cols = [
        "Idade_anos",
        "Densidade_pl_ha",
        "MO_g_dm3",
        "PST",
        "Ds_g_cm3",
    ]

    available_cols += [col for col in extra_cols if col in df.columns]
    available_cols = list(dict.fromkeys(available_cols))

    sub = df[available_cols].dropna(subset=[FARM]).copy()

    agg = {}
    for col in available_cols:
        if col == FARM:
            continue
        agg[f"{col}_mean"] = (col, "mean")
        agg[f"{col}_median"] = (col, "median")

    out = sub.groupby(FARM).agg(
        n=(RESPONSE, "size"),
        **agg,
    ).reset_index()

    return out


def main():
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Input data not found: {DATA_PATH}")

    df = pd.read_csv(DATA_PATH)

    required = [RESPONSE, FARM]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError("Missing required columns: " + ", ".join(missing))

    enzymes = [enzyme for enzyme in ENZYMES if enzyme in df.columns]

    model_specs = [
        ("enzyme_only", False, []),
        ("farm_fixed", True, []),
        ("farm_fixed_plus_structure", True, STRUCTURE_COVARIATES),
        ("farm_fixed_plus_key_soil", True, KEY_SOIL_COVARIATES),
        ("farm_fixed_plus_structure_and_key_soil", True, STRUCTURE_COVARIATES + KEY_SOIL_COVARIATES),
    ]

    model_rows = []

    for enzyme in enzymes:
        for model_name, include_farm, covariates in model_specs:
            row = fit_model(
                df,
                enzyme=enzyme,
                model_name=model_name,
                include_farm=include_farm,
                covariates=covariates,
            )
            if row is not None:
                model_rows.append(row)

    model_df = pd.DataFrame(model_rows)
    model_df.to_csv(OUT_MODELS, index=False)

    partial_specs = [
        ("none", False, []),
        ("farm", True, []),
        ("farm_plus_structure", True, STRUCTURE_COVARIATES),
        ("farm_plus_key_soil", True, KEY_SOIL_COVARIATES),
        ("farm_plus_structure_and_key_soil", True, STRUCTURE_COVARIATES + KEY_SOIL_COVARIATES),
    ]

    partial_rows = []

    for enzyme in enzymes:
        for adjustment_name, include_farm, covariates in partial_specs:
            row = partial_spearman(
                df,
                enzyme=enzyme,
                adjustment_name=adjustment_name,
                include_farm=include_farm,
                covariates=covariates,
            )
            if row is not None:
                partial_rows.append(row)

    partial_df = pd.DataFrame(partial_rows)
    partial_df.to_csv(OUT_PARTIAL, index=False)

    farm_summary = farm_level_summary(df)
    farm_summary.to_csv(OUT_FARM_SUMMARY, index=False)

    print("Farm-effect enzyme analysis completed.")
    print(f"Model output: {OUT_MODELS}")
    print(f"Partial correlation output: {OUT_PARTIAL}")
    print(f"Farm-level summary: {OUT_FARM_SUMMARY}")

    print("\nFarm counts:")
    print(df[FARM].value_counts().sort_index().to_string())

    print("\nFarm fixed-effect model summary:")
    display_cols = [
        "enzyme",
        "model_name",
        "n",
        "n_farms",
        "r2",
        "adj_r2",
        "enzyme_coef",
        "enzyme_p_hc3",
        "covariates",
        "fit_status",
    ]
    print(model_df[display_cols].to_string(index=False))

    print("\nFarm-adjusted partial Spearman summary:")
    display_cols = [
        "enzyme",
        "adjustment_name",
        "n",
        "partial_spearman_rho",
        "partial_spearman_p",
        "covariates",
    ]
    print(partial_df[display_cols].to_string(index=False))


if __name__ == "__main__":
    main()
