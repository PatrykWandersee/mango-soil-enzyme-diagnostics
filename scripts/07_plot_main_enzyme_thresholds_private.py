from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats


DATA_PATH = Path("data/processed/private/enzyme_diagnostics_private.csv")
SUMMARY_PATH = Path("tables/private/main_enzyme_threshold_summary_private.csv")
FIGURE_DIR = Path("figures/private")
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

RESPONSE = "Prod_rel_pct"

PLOT_INDICATORS = {
    "GMea": "Geometric mean of enzyme activities (GMea)",
    "Beta_glic": "Beta-glucosidase activity",
    "Arilsulf": "Arylsulfatase activity",
    "GMea_per_Clay": "Clay-normalized GMea",
}


def make_plot(df, summary, indicator, label):
    sub = df[[indicator, RESPONSE]].dropna().copy()
    row = summary.loc[summary["indicator"] == indicator].iloc[0]

    threshold_tree = row["threshold_tree"]
    threshold_cate = row["threshold_cate_nelson"]
    p25 = row["high_yield_p25"]
    p75 = row["high_yield_p75"]

    fig, ax = plt.subplots(figsize=(6, 4.5))

    ax.scatter(sub[indicator], sub[RESPONSE], alpha=0.85)

    slope, intercept, r_value, p_value, std_err = stats.linregress(
        sub[indicator], sub[RESPONSE]
    )
    x_min = sub[indicator].min()
    x_max = sub[indicator].max()
    x_line = [x_min, x_max]
    y_line = [intercept + slope * x_min, intercept + slope * x_max]
    ax.plot(x_line, y_line)

    ax.axhline(70, linestyle="--", linewidth=1)
    ax.axvline(threshold_tree, linestyle="--", linewidth=1)
    ax.axvline(threshold_cate, linestyle=":", linewidth=1)
    ax.axvspan(p25, p75, alpha=0.15)

    ax.set_xlabel(label)
    ax.set_ylabel("Relative yield per plant (%)")

    ax.text(
        0.04,
        0.96,
        (
            f"Spearman rho = {row['spearman_rho']:.3f}\n"
            f"Tree threshold = {threshold_tree:.3f}\n"
            f"Cate-Nelson threshold = {threshold_cate:.3f}\n"
            f"High-yield IQR = {p25:.3f}–{p75:.3f}"
        ),
        transform=ax.transAxes,
        va="top",
        ha="left",
    )

    fig.tight_layout()

    output_png = FIGURE_DIR / f"{indicator}_threshold_diagnostic_private.png"
    output_pdf = FIGURE_DIR / f"{indicator}_threshold_diagnostic_private.pdf"

    fig.savefig(output_png, dpi=300)
    fig.savefig(output_pdf)
    plt.close(fig)

    return output_png, output_pdf


def main():
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Input data not found: {DATA_PATH}")
    if not SUMMARY_PATH.exists():
        raise FileNotFoundError(f"Summary table not found: {SUMMARY_PATH}")

    df = pd.read_csv(DATA_PATH)
    summary = pd.read_csv(SUMMARY_PATH)

    print("Creating private enzyme threshold diagnostic plots.")

    for indicator, label in PLOT_INDICATORS.items():
        png, pdf = make_plot(df, summary, indicator, label)
        print(f"- {indicator}: {png}, {pdf}")

    print("\nDone. Private figures were saved to figures/private/.")


if __name__ == "__main__":
    main()
