from pathlib import Path

import pandas as pd


INPUT_PATH = Path("data/private/soil_quality_processed_private.csv")
OUTPUT_PATH = Path("data/processed/private/enzyme_diagnostics_private.csv")
SUMMARY_PATH = Path("tables/private/enzyme_diagnostics_summary_private.csv")

OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
SUMMARY_PATH.parent.mkdir(parents=True, exist_ok=True)


ID_COLUMNS = [
    "Area",
    "Fazenda",
]

STRUCTURAL_COLUMNS = [
    "Idade_anos",
    "Densidade_pl_ha",
]

PRIMARY_RESPONSE_COLUMNS = [
    "Prod_rel_pct",
]

SENSITIVITY_RESPONSE_COLUMNS = [
    "Prod_rel_ha_pct",
]

RAW_ENZYME_COLUMNS = [
    "Beta_glic",
    "Arilsulf",
    "GMea",
    "Ativ_Enzimat_Total",
    "Rel_Beta_Aril",
]

NORMALIZED_ENZYME_COLUMNS = [
    "Beta_por_Argila",
    "Aril_por_Argila",
    "GMea_por_Argila",
    "Ativ_Enzim_por_MO",
    "GMea_per_Clay",
    "log_GMea_per_Clay",
]

BIOLOGICAL_QUOTIENT_COLUMNS = [
    "qBeta",
    "qAril",
    "qGMea",
]

ORGANIC_MATTER_COLUMNS = [
    "MO_g_dm3",
    "C_org_g_d3",
]

TEXTURE_STRUCTURE_COLUMNS = [
    "Areia_g_kg",
    "Argila_g_kg",
    "Silte_g_kg",
    "Floculacao_pct",
    "Dispercao_pct",
    "Ds_g_cm3",
]

CHEMICAL_CONTEXT_COLUMNS = [
    "pH",
    "CE_dS_m",
    "PST",
    "PM1_mg_dm3",
    "SB_cmolc_Kg",
    "Ca_Troc_cmolc_Kg",
    "Mg_Troc_cmolc_Kg",
    "K_Troc_cmolc_Kg",
    "Na_Troc_cmolc_Kg",
    "Na_pct_sol",
    "RAS",
    "CROSS",
]

COLUMN_GROUPS = {
    "id": ID_COLUMNS,
    "structural": STRUCTURAL_COLUMNS,
    "primary_response": PRIMARY_RESPONSE_COLUMNS,
    "sensitivity_response": SENSITIVITY_RESPONSE_COLUMNS,
    "raw_enzymes": RAW_ENZYME_COLUMNS,
    "normalized_enzymes": NORMALIZED_ENZYME_COLUMNS,
    "biological_quotients": BIOLOGICAL_QUOTIENT_COLUMNS,
    "organic_matter": ORGANIC_MATTER_COLUMNS,
    "texture_structure": TEXTURE_STRUCTURE_COLUMNS,
    "chemical_context": CHEMICAL_CONTEXT_COLUMNS,
}


def flatten_column_groups():
    columns = []
    for group_columns in COLUMN_GROUPS.values():
        columns.extend(group_columns)
    return list(dict.fromkeys(columns))


def main():
    if not INPUT_PATH.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_PATH}")

    df = pd.read_csv(INPUT_PATH)
    selected_columns = flatten_column_groups()

    missing_columns = [col for col in selected_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(
            "The following selected columns are missing from the input dataset:\n"
            + "\n".join(f"- {col}" for col in missing_columns)
        )

    out = df[selected_columns].copy()
    out.to_csv(OUTPUT_PATH, index=False)

    rows = []
    for group_name, columns in COLUMN_GROUPS.items():
        for col in columns:
            series = out[col]
            rows.append(
                {
                    "group": group_name,
                    "column": col,
                    "dtype": str(series.dtype),
                    "non_missing": int(series.notna().sum()),
                    "missing": int(series.isna().sum()),
                    "missing_pct": float(series.isna().mean() * 100),
                    "min": float(series.min()) if pd.api.types.is_numeric_dtype(series) else None,
                    "mean": float(series.mean()) if pd.api.types.is_numeric_dtype(series) else None,
                    "median": float(series.median()) if pd.api.types.is_numeric_dtype(series) else None,
                    "max": float(series.max()) if pd.api.types.is_numeric_dtype(series) else None,
                }
            )

    summary = pd.DataFrame(rows)
    summary.to_csv(SUMMARY_PATH, index=False)

    print("Private enzyme diagnostics dataset created.")
    print(f"Input shape: {df.shape[0]} rows x {df.shape[1]} columns")
    print(f"Output shape: {out.shape[0]} rows x {out.shape[1]} columns")
    print(f"Output file: {OUTPUT_PATH}")
    print(f"Summary file: {SUMMARY_PATH}")

    print("\nSelected column groups:")
    for group_name, columns in COLUMN_GROUPS.items():
        print(f"\n{group_name}:")
        for col in columns:
            print(f"- {col}")


if __name__ == "__main__":
    main()
