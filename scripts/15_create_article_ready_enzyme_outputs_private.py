from pathlib import Path
import shutil

import pandas as pd


TABLE_DIR = Path("tables/private")
FIGURE_DIR = Path("figures/private")

OUT_TABLE_DIR = TABLE_DIR / "manuscript"
OUT_FIGURE_DIR = FIGURE_DIR / "manuscript"

OUT_TABLE_DIR.mkdir(parents=True, exist_ok=True)
OUT_FIGURE_DIR.mkdir(parents=True, exist_ok=True)


INPUT_DIAGNOSTIC_TABLE = TABLE_DIR / "article_enzyme_diagnostic_synthesis_private.csv"
INPUT_MODEL_TABLE = TABLE_DIR / "article_enzyme_model_support_private.csv"

OUTPUT_DIAGNOSTIC_TABLE_CSV = OUT_TABLE_DIR / "table_enzyme_diagnostic_thresholds_private.csv"
OUTPUT_DIAGNOSTIC_TABLE_MD = OUT_TABLE_DIR / "table_enzyme_diagnostic_thresholds_private.md"

OUTPUT_MODEL_TABLE_CSV = OUT_TABLE_DIR / "table_enzyme_model_support_private.csv"
OUTPUT_MODEL_TABLE_MD = OUT_TABLE_DIR / "table_enzyme_model_support_private.md"

MANIFEST_PATH = OUT_TABLE_DIR / "article_ready_enzyme_outputs_manifest_private.md"

FIGURE_SOURCES = {
    "main_vertical_png": (
        FIGURE_DIR / "enzyme_threshold_panel_vertical.png",
        OUT_FIGURE_DIR / "figure_enzyme_diagnostic_panel_vertical.png",
    ),
    "main_vertical_pdf": (
        FIGURE_DIR / "enzyme_threshold_panel_vertical.pdf",
        OUT_FIGURE_DIR / "figure_enzyme_diagnostic_panel_vertical.pdf",
    ),
    "supplement_horizontal_png": (
        FIGURE_DIR / "enzyme_threshold_panel_horizontal.png",
        OUT_FIGURE_DIR / "figure_enzyme_diagnostic_panel_horizontal.png",
    ),
    "supplement_horizontal_pdf": (
        FIGURE_DIR / "enzyme_threshold_panel_horizontal.pdf",
        OUT_FIGURE_DIR / "figure_enzyme_diagnostic_panel_horizontal.pdf",
    ),
}


def to_markdown_table(df, path):
    lines = []
    cols = [str(col) for col in df.columns]

    lines.append("| " + " | ".join(cols) + " |")
    lines.append("| " + " | ".join(["---"] * len(cols)) + " |")

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


def copy_if_exists(src, dst):
    if src.exists():
        shutil.copy2(src, dst)
        return True
    return False


def main():
    generated = []
    missing = []

    if not INPUT_DIAGNOSTIC_TABLE.exists():
        missing.append(str(INPUT_DIAGNOSTIC_TABLE))
    else:
        diagnostic = pd.read_csv(INPUT_DIAGNOSTIC_TABLE)

        # Keep the main table focused on the diagnostic message.
        diagnostic_cols = [
            "indicator",
            "role",
            "Spearman rho",
            "Spearman p",
            "Cate-Nelson threshold",
            "Cate-Nelson bootstrap median (95% interval)",
            "High-yield IQR",
            "Balanced accuracy",
            "Interpretation",
        ]
        diagnostic_cols = [col for col in diagnostic_cols if col in diagnostic.columns]
        diagnostic_out = diagnostic[diagnostic_cols].copy()

        diagnostic_out.to_csv(OUTPUT_DIAGNOSTIC_TABLE_CSV, index=False)
        to_markdown_table(diagnostic_out, OUTPUT_DIAGNOSTIC_TABLE_MD)

        generated.extend(
            [
                OUTPUT_DIAGNOSTIC_TABLE_CSV,
                OUTPUT_DIAGNOSTIC_TABLE_MD,
            ]
        )

    if not INPUT_MODEL_TABLE.exists():
        missing.append(str(INPUT_MODEL_TABLE))
    else:
        model = pd.read_csv(INPUT_MODEL_TABLE)

        # Keep model support as a secondary table, not the main diagnostic table.
        model_cols = [
            "indicator",
            "best response model",
            "best model R2",
            "best model RMSE",
            "best model breakpoint or vertex",
            "competitive models (delta AICc <= 2)",
            "linear-plateau bootstrap median (95% interval)",
            "orchard-structure adjusted enzyme effect",
            "soil-core adjusted enzyme effect",
            "farm fixed-effect enzyme effect",
            "farm + structure enzyme effect",
            "farm + structure + key soil enzyme effect",
        ]
        model_cols = [col for col in model_cols if col in model.columns]
        model_out = model[model_cols].copy()

        model_out.to_csv(OUTPUT_MODEL_TABLE_CSV, index=False)
        to_markdown_table(model_out, OUTPUT_MODEL_TABLE_MD)

        generated.extend(
            [
                OUTPUT_MODEL_TABLE_CSV,
                OUTPUT_MODEL_TABLE_MD,
            ]
        )

    for _, (src, dst) in FIGURE_SOURCES.items():
        if copy_if_exists(src, dst):
            generated.append(dst)
        else:
            missing.append(str(src))

    manifest_lines = [
        "# Article-ready private enzyme outputs",
        "",
        "This folder contains private article-ready outputs generated from the private enzyme diagnostics workflow.",
        "",
        "## Intended use",
        "",
        "- `figure_enzyme_diagnostic_panel_vertical.*`: candidate main diagnostic figure.",
        "- `figure_enzyme_diagnostic_panel_horizontal.*`: candidate supplementary or presentation figure.",
        "- `table_enzyme_diagnostic_thresholds_private.*`: candidate main diagnostic table.",
        "- `table_enzyme_model_support_private.*`: candidate supplementary model-support table.",
        "",
        "## Interpretation notes",
        "",
        "- The vertical diagnostic line represents the Cate–Nelson threshold using the 70% relative yield criterion.",
        "- The shaded band represents the interquartile range of enzyme values among orchards with relative yield per plant ≥70%.",
        "- The shaded band is a dataset-specific high-yield-associated range, not a universal optimum range.",
        "- Plateau breakpoints support non-linear interpretation but should not be treated as primary diagnostic thresholds unless stability is acceptable.",
        "",
        "## Generated files",
        "",
    ]

    for path in generated:
        manifest_lines.append(f"- `{path}`")

    if missing:
        manifest_lines.extend(
            [
                "",
                "## Missing expected source files",
                "",
            ]
        )
        for path in missing:
            manifest_lines.append(f"- `{path}`")

    MANIFEST_PATH.write_text("\n".join(manifest_lines) + "\n")
    generated.append(MANIFEST_PATH)

    print("Article-ready private enzyme outputs created.")
    print("\nGenerated files:")
    for path in generated:
        print(f"- {path}")

    if missing:
        print("\nMissing expected source files:")
        for path in missing:
            print(f"- {path}")


if __name__ == "__main__":
    main()
