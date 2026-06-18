from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats


DATA_PATH = Path("data/processed/private/enzyme_diagnostics_private.csv")
SUMMARY_PATH = Path("tables/private/main_enzyme_threshold_summary_private.csv")
FIGURE_DIR = Path("figures/private")
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

RESPONSE = "Prod_rel_pct"
YIELD_CUTOFF = 70

PLOT_ORDER = [
    ("Beta_glic", "Beta-glucosidase activity"),
    ("Arilsulf", "Arylsulfatase activity"),
    ("GMea", "GMea"),
]

PANEL_LABELS = ["A", "B", "C"]


def format_value(value):
    if abs(value) >= 10:
        return f"{value:.1f}"
    return f"{value:.3f}"


def main():
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Input data not found: {DATA_PATH}")
    if not SUMMARY_PATH.exists():
        raise FileNotFoundError(f"Summary table not found: {SUMMARY_PATH}")

    df = pd.read_csv(DATA_PATH)
    summary = pd.read_csv(SUMMARY_PATH)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharey=True)

    for ax, panel_label, (indicator, x_label) in zip(axes, PANEL_LABELS, PLOT_ORDER):
        sub = df[[indicator, RESPONSE]].dropna().copy()
        row = summary.loc[summary["indicator"] == indicator].iloc[0]

        cate_threshold = row["threshold_cate_nelson"]
        tree_threshold = row["threshold_tree"]
        high_p25 = row["high_yield_p25"]
        high_p75 = row["high_yield_p75"]
        rho = row["spearman_rho"]

        ax.scatter(sub[indicator], sub[RESPONSE], alpha=0.82)

        ax.axhline(YIELD_CUTOFF, linestyle="--", linewidth=1)
        ax.axvline(cate_threshold, linestyle="--", linewidth=1)
        ax.axvspan(high_p25, high_p75, alpha=0.15)

        # Show tree threshold only when it is materially different from Cate-Nelson.
        if abs(tree_threshold - cate_threshold) > 0.05 * max(abs(cate_threshold), 1):
            ax.axvline(tree_threshold, linestyle=":", linewidth=1)

        ax.set_xlabel(x_label)

        ax.text(
            0.04,
            0.96,
            panel_label,
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontweight="bold",
        )

        ax.text(
            0.04,
            0.86,
            (
                f"ρ = {rho:.3f}\n"
                f"Threshold = {format_value(cate_threshold)}\n"
                f"IQR ≥70% = {format_value(high_p25)}–{format_value(high_p75)}"
            ),
            transform=ax.transAxes,
            va="top",
            ha="left",
        )

    axes[0].set_ylabel("Relative yield per plant (%)")

    fig.tight_layout()

    output_png = FIGURE_DIR / "three_enzyme_diagnostic_panel_private.png"
    output_pdf = FIGURE_DIR / "three_enzyme_diagnostic_panel_private.pdf"

    fig.savefig(output_png, dpi=300)
    fig.savefig(output_pdf)
    plt.close(fig)

    print("Three-enzyme diagnostic panel created.")
    print(f"PNG: {output_png}")
    print(f"PDF: {output_pdf}")


if __name__ == "__main__":
    main()
