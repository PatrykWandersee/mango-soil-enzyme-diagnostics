from pathlib import Path

import pandas as pd


DATA_PATH = Path("data/private/soil_quality_processed_private.csv")
OUT_DIR = Path("tables/private")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def find_columns(columns, keywords):
    matches = []
    for col in columns:
        col_lower = col.lower()
        if any(keyword.lower() in col_lower for keyword in keywords):
            matches.append(col)
    return matches


def main():
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Private dataset not found: {DATA_PATH}")

    df = pd.read_csv(DATA_PATH)

    print("\nDataset loaded successfully.")
    print(f"Rows: {df.shape[0]}")
    print(f"Columns: {df.shape[1]}")

    print("\nColumns:")
    for col in df.columns:
        print(f"- {col}")

    overview = pd.DataFrame({
        "column": df.columns,
        "dtype": [str(df[col].dtype) for col in df.columns],
        "non_missing": [df[col].notna().sum() for col in df.columns],
        "missing": [df[col].isna().sum() for col in df.columns],
        "missing_pct": [df[col].isna().mean() * 100 for col in df.columns],
        "unique_values": [df[col].nunique(dropna=True) for col in df.columns],
    })

    overview.to_csv(OUT_DIR / "dataset_column_overview_private.csv", index=False)

    candidate_groups = {
        "biological": ["beta", "glic", "gluc", "aril", "aryl", "gmea", "enzyme", "enzim"],
        "texture_structure": ["argila", "clay", "areia", "sand", "floc", "ds", "density", "dens"],
        "organic_matter": ["mo", "som", "organic"],
        "salinity_sodicity": ["ce", "ece", "cond", "na", "sod", "pst", "esp"],
        "fertility": ["ca", "mg", "k", "p_", "pm1", "mehlich", "sb", "ph"],
        "yield": ["prod", "yield", "rel"],
        "farm_group": ["fazenda", "farm", "area", "grupo"],
    }

    rows = []
    for group, keywords in candidate_groups.items():
        for col in find_columns(df.columns, keywords):
            rows.append({"group": group, "column": col})

    candidates = pd.DataFrame(rows).drop_duplicates()
    candidates.to_csv(OUT_DIR / "candidate_columns_private.csv", index=False)

    numeric_cols = df.select_dtypes(include="number").columns
    numeric_summary = df[numeric_cols].describe().T
    numeric_summary.to_csv(OUT_DIR / "numeric_summary_private.csv")

    print("\nPrivate inspection outputs saved to:")
    print(f"- {OUT_DIR / 'dataset_column_overview_private.csv'}")
    print(f"- {OUT_DIR / 'candidate_columns_private.csv'}")
    print(f"- {OUT_DIR / 'numeric_summary_private.csv'}")

    print("\nCandidate column groups:")
    if candidates.empty:
        print("No candidate columns detected by keyword search.")
    else:
        for group in candidates["group"].unique():
            cols = candidates.loc[candidates["group"] == group, "column"].tolist()
            print(f"\n{group}:")
            for col in cols:
                print(f"- {col}")


if __name__ == "__main__":
    main()
