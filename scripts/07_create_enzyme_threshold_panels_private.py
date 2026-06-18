"""Create manuscript-style enzyme threshold panel figures.

This script creates article-ready panel figures for the three main enzyme
indicators using all observations:
- beta-glucosidase
- arylsulfatase
- GMea

Displayed threshold:
- regression-tree threshold (used here as the critical level)

Displayed range:
- interquartile range (IQR) of observations with relative yield >= 70%

Inputs
------
data/processed/private/enzyme_diagnostics_private.csv
tables/private/main_enzyme_threshold_summary_private.csv

Outputs
-------
figures/private/enzyme_threshold_panel_horizontal.png
figures/private/enzyme_threshold_panel_horizontal.pdf
figures/private/enzyme_threshold_panel_vertical.png
figures/private/enzyme_threshold_panel_vertical.pdf
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


DATA_PATH = Path("data/processed/private/enzyme_diagnostics_private.csv")
SUMMARY_PATH = Path("tables/private/main_enzyme_threshold_summary_private.csv")
OUTPUT_DIR = Path("figures/private")

RESPONSE_COLUMN = "Prod_rel_pct"
YIELD_CUTOFF = 70.0

PANEL_CONFIG = [
    {
        "indicator": "Beta_glic",
        "panel_label": "A",
        "title": r"$\beta$-glucosidase",
        "xlabel": r"$\beta$-glucosidase activity ($\mu$g p-nit. g$^{-1}$ h$^{-1}$)",
    },
    {
        "indicator": "Arilsulf",
        "panel_label": "B",
        "title": "Arylsulfatase",
        "xlabel": r"Arylsulfatase activity ($\mu$g p-nit. g$^{-1}$ h$^{-1}$)",
    },
    {
        "indicator": "GMea",
        "panel_label": "C",
        "title": "GMea",
        "xlabel": r"GMea ($\mu$g p-nit. g$^{-1}$ h$^{-1}$)",
    },
]

REQUIRED_SUMMARY_COLUMNS = [
    "indicator",
    "spearman_rho",
    "threshold_tree",
    "high_yield_p25",
    "high_yield_p75",
]

POINT_COLOR = "#4C78A8"
THRESHOLD_COLOR = "#D62728"
YIELD_LINE_COLOR = "#1F77B4"
IQR_COLOR = "#B0BEC5"


def validate_columns(df: pd.DataFrame, required_columns: list[str], df_name: str) -> None:
    """Validate that a dataframe contains all required columns."""
    missing_columns = [column for column in required_columns if column not in df.columns]
    if missing_columns:
        missing_text = ", ".join(missing_columns)
        raise ValueError(f"Missing columns in {df_name}: {missing_text}")


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load dataset and threshold summary."""
    data = pd.read_csv(DATA_PATH)
    summary = pd.read_csv(SUMMARY_PATH)

    validate_columns(
        data,
        [RESPONSE_COLUMN, *[config["indicator"] for config in PANEL_CONFIG]],
        "enzyme dataset",
    )
    validate_columns(summary, REQUIRED_SUMMARY_COLUMNS, "threshold summary")

    return data, summary


def get_summary_lookup(summary: pd.DataFrame) -> dict[str, pd.Series]:
    """Create a lookup dictionary keyed by indicator."""
    lookup = {}
    for config in PANEL_CONFIG:
        indicator = config["indicator"]
        subset = summary[summary["indicator"] == indicator]
        if subset.empty:
            raise ValueError(f"No summary row found for indicator: {indicator}")
        lookup[indicator] = subset.iloc[0]
    return lookup


def format_annotation(row: pd.Series) -> str:
    """Format compact annotation text for one panel."""
    return (
        f"\u03C1 = {row['spearman_rho']:.3f}\n"
        f"Diagnostic threshold = {row['threshold_cate_nelson']:.1f}\n"
        f"IQR of ≥70% yield group = {row['high_yield_p25']:.1f}\u2013{row['high_yield_p75']:.1f}"
    )


def set_x_limits(ax, x_values: pd.Series, row: pd.Series) -> None:
    """Set x-limits with padding around data and summary landmarks."""
    x_min = min(
        x_values.min(),
        row["threshold_cate_nelson"],
        row["high_yield_p25"],
        row["high_yield_p75"],
    )
    x_max = max(
        x_values.max(),
        row["threshold_cate_nelson"],
        row["high_yield_p25"],
        row["high_yield_p75"],
    )
    x_range = x_max - x_min
    padding = x_range * 0.10 if x_range > 0 else 1.0
    ax.set_xlim(x_min - padding, x_max + padding)


def create_single_panel(
    ax,
    data: pd.DataFrame,
    summary_row: pd.Series,
    config: dict,
    show_ylabel: bool,
) -> None:
    """Draw one enzyme threshold panel."""
    indicator = config["indicator"]
    x = data[indicator]
    y = data[RESPONSE_COLUMN]

    ax.axvspan(
        summary_row["high_yield_p25"],
        summary_row["high_yield_p75"],
        color=IQR_COLOR,
        alpha=0.45,
        zorder=0,
    )
    ax.axhline(
        YIELD_CUTOFF,
        color=YIELD_LINE_COLOR,
        linestyle=":",
        linewidth=1.6,
        zorder=1,
    )
    ax.axvline(
        summary_row["threshold_cate_nelson"],
        color=THRESHOLD_COLOR,
        linestyle="--",
        linewidth=1.6,
        zorder=1,
    )
    ax.scatter(
        x,
        y,
        s=70,
        color=POINT_COLOR,
        alpha=0.90,
        edgecolors="white",
        linewidth=0.5,
        zorder=2,
    )

    ax.set_title(config["title"], fontsize=13, pad=10)
    ax.set_xlabel(config["xlabel"], fontsize=11)

    if show_ylabel:
        ax.set_ylabel("Relative yield per plant (%)", fontsize=11)

    ax.set_ylim(5, 105)
    set_x_limits(ax, x, summary_row)
    ax.grid(True, alpha=0.30)

    ax.text(
        0.04,
        0.96,
        config["panel_label"],
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
    )

    ax.text(
        0.06,
        0.84,
        format_annotation(summary_row),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10.5,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.70, "edgecolor": "none"},
    )


def create_panel_figure(
    data: pd.DataFrame,
    summary_lookup: dict[str, pd.Series],
    layout: str,
) -> None:
    """Create and save one panel figure layout."""
    if layout == "horizontal":
        fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.9), sharey=True)
        output_stem = "enzyme_threshold_panel_horizontal"
        legend_ncol = 4
    elif layout == "vertical":
        fig, axes = plt.subplots(3, 1, figsize=(6.8, 11.8), sharey=True)
        output_stem = "enzyme_threshold_panel_vertical"
        legend_ncol = 2
    else:
        raise ValueError(f"Unsupported layout: {layout}")

    if not isinstance(axes, (list, tuple)):
        axes = axes.flatten()

    for index, (ax, config) in enumerate(zip(axes, PANEL_CONFIG)):
        summary_row = summary_lookup[config["indicator"]]
        create_single_panel(
            ax=ax,
            data=data,
            summary_row=summary_row,
            config=config,
            show_ylabel=(index == 0),
        )

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="None",
            markerfacecolor=POINT_COLOR,
            markeredgecolor="white",
            markeredgewidth=0.6,
            markersize=8,
            label="Observation",
        ),
        Line2D(
            [0],
            [0],
            color=THRESHOLD_COLOR,
            linestyle="--",
            linewidth=1.6,
            label="Diagnostic threshold",
        ),
        Line2D(
            [0],
            [0],
            color=YIELD_LINE_COLOR,
            linestyle=":",
            linewidth=1.6,
            label=f"Relative yield cutoff ({YIELD_CUTOFF:.0f}%)",
        ),
        Patch(
            facecolor=IQR_COLOR,
            alpha=0.45,
            label=f"IQR of ≥{YIELD_CUTOFF:.0f}% yield group",
        ),
    ]

    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=legend_ncol,
        frameon=False,
        bbox_to_anchor=(0.5, -0.01),
    )

    fig.tight_layout(rect=[0, 0.06, 1, 1])

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    png_path = OUTPUT_DIR / f"{output_stem}.png"
    pdf_path = OUTPUT_DIR / f"{output_stem}.pdf"

    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    print(f"Created: {png_path}")
    print(f"Created: {pdf_path}")


def main() -> None:
    """Run the full figure-generation workflow."""
    data, summary = load_inputs()
    summary_lookup = get_summary_lookup(summary)

    create_panel_figure(data=data, summary_lookup=summary_lookup, layout="horizontal")
    create_panel_figure(data=data, summary_lookup=summary_lookup, layout="vertical")

    print()
    print("Enzyme threshold panel figures created successfully.")
    print(f"Input dataset: {DATA_PATH}")
    print(f"Input summary: {SUMMARY_PATH}")
    print(f"Output directory: {OUTPUT_DIR}")
    print()
    print("Thresholds used in the figures:")
    for config in PANEL_CONFIG:
        row = summary_lookup[config["indicator"]]
        print(
            f"- {config['indicator']}: "
            f"tree threshold = {row['threshold_cate_nelson']:.3f}; "
            f"high-yield IQR = {row['high_yield_p25']:.3f} to {row['high_yield_p75']:.3f}; "
            f"Spearman rho = {row['spearman_rho']:.3f}"
        )


if __name__ == "__main__":
    main()
