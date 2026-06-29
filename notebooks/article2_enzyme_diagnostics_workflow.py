# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Article 2 workflow: soil enzyme diagnostics
#
# This notebook is a narrative orchestrator for the Article 2 workflow in this
# repository. It does not reimplement the analyses. Instead, it documents the
# objective, checks required private inputs, summarizes existing methodological
# documents, runs the existing scripts in order, and displays the final
# article-ready tables and figures.
#
# The workflow supports the manuscript described in the repository as:
#
# **Soil enzyme diagnostics for irrigated Palmer mango orchards in the Brazilian semiarid region**
#
# Main indicators and interpretation rules are defined in the project documents
# under `docs/`.

# %% [markdown]
# ## Execution policy
#
# This notebook is intended to be run from the repository root or from the
# `notebooks/` directory.
#
# It does not modify raw data. The analysis scripts read private inputs and
# write processed data, tables, and figures under the project output folders.
# The original private dataset remains under `data/private/`.
#
# If an expected file or script is missing, the notebook reports a clear TODO
# instead of assuming an alternative path.

# %%
from __future__ import annotations

import platform
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
from IPython.display import Image, Markdown, display


def find_repo_root(start: Path | None = None) -> Path:
    """Find the repository root from the current working directory."""
    current = (start or Path.cwd()).resolve()
    candidates = [current, *current.parents]
    for candidate in candidates:
        if (candidate / "README.md").exists() and (candidate / "scripts").is_dir():
            return candidate
    raise RuntimeError("TODO: repository root was not found from the current directory.")


REPO_ROOT = find_repo_root()
print(f"Repository root: {REPO_ROOT}")
print(f"Python: {sys.version.split()[0]}")
print(f"Platform: {platform.platform()}")
print(f"Execution timestamp: {datetime.now().isoformat(timespec='seconds')}")

# %% [markdown]
# ## Required private inputs
#
# The repository documentation states that the original field dataset is private.
# This workflow expects the private input dataset at the path used by the
# existing preparation script.

# %%
PRIVATE_INPUTS = [
    REPO_ROOT / "data/private/soil_quality_processed_private.csv",
]

EXPECTED_PREPARED_OUTPUTS = [
    REPO_ROOT / "data/processed/private/enzyme_diagnostics_private.csv",
]

missing_private_inputs = [path for path in PRIVATE_INPUTS if not path.exists()]

if missing_private_inputs:
    display(Markdown("### TODO: missing private input files"))
    for path in missing_private_inputs:
        display(Markdown(f"- `{path.relative_to(REPO_ROOT)}`"))
    raise FileNotFoundError("Private input files are missing. See TODO list above.")

display(Markdown("All required private input files were found:"))
for path in PRIVATE_INPUTS:
    display(Markdown(f"- `{path.relative_to(REPO_ROOT)}`"))

display(Markdown("Expected product after data preparation:"))
for path in EXPECTED_PREPARED_OUTPUTS:
    status = "found" if path.exists() else "TODO: not generated yet"
    display(Markdown(f"- `{path.relative_to(REPO_ROOT)}`: {status}"))

# %% [markdown]
# ## Methodological documents
#
# The cells below display concise excerpts from existing project documents. These
# documents define the article scope, analysis strategy, methodological
# decisions, and manuscript outline. The notebook uses these files as the source
# of the narrative context instead of duplicating their content.

# %%
DOCS = [
    REPO_ROOT / "README.md",
    REPO_ROOT / "docs/article_scope.md",
    REPO_ROOT / "docs/analysis_strategy.md",
    REPO_ROOT / "docs/methodological_decisions.md",
    REPO_ROOT / "docs/manuscript_outline.md",
]


def show_markdown_excerpt(path: Path, max_chars: int = 3000) -> None:
    if not path.exists():
        display(Markdown(f"### TODO: missing document `{path.relative_to(REPO_ROOT)}`"))
        return

    text = path.read_text(encoding="utf-8")
    excerpt = text[:max_chars]
    if len(text) > max_chars:
        excerpt = excerpt.rstrip() + "\n\n..."

    display(Markdown(f"## `{path.relative_to(REPO_ROOT)}`"))
    display(Markdown(excerpt))


for doc_path in DOCS:
    show_markdown_excerpt(doc_path)

# %% [markdown]
# ## Pipeline scripts
#
# The workflow below executes the existing scripts in their numeric order. The
# list is explicit so the notebook remains auditable. The two scripts numbered
# `07` are both retained because both exist in the repository.
#
# The scripts are run with the repository root as the working directory.

# %%
PIPELINE_SCRIPTS = [
    REPO_ROOT / "scripts/00_inspect_private_dataset.py",
    REPO_ROOT / "scripts/01_prepare_enzyme_dataset_private.py",
    REPO_ROOT / "scripts/02_explore_enzyme_relationships_private.py",
    REPO_ROOT / "scripts/03_screen_enzyme_thresholds_private.py",
    REPO_ROOT / "scripts/04_summarize_threshold_stability_private.py",
    REPO_ROOT / "scripts/05_experimental_group_sensitivity_private.py",
    REPO_ROOT / "scripts/06_create_main_threshold_summary_private.py",
    REPO_ROOT / "scripts/07_create_enzyme_threshold_panels_private.py",
    REPO_ROOT / "scripts/07_plot_main_enzyme_thresholds_private.py",
    REPO_ROOT / "scripts/08_plot_three_enzyme_diagnostic_panel_private.py",
    REPO_ROOT / "scripts/09_compare_enzyme_response_models_private.py",
    REPO_ROOT / "scripts/10_bootstrap_enzyme_thresholds_private.py",
    REPO_ROOT / "scripts/11_analyze_enzyme_soil_context_private.py",
    REPO_ROOT / "scripts/12_analyze_enzyme_signal_after_orchard_structure_private.py",
    REPO_ROOT / "scripts/13_analyze_enzyme_farm_effects_private.py",
    REPO_ROOT / "scripts/14_create_article_enzyme_synthesis_private.py",
    REPO_ROOT / "scripts/15_create_article_ready_enzyme_outputs_private.py",
]

missing_scripts = [path for path in PIPELINE_SCRIPTS if not path.exists()]

if missing_scripts:
    display(Markdown("### TODO: missing pipeline scripts"))
    for path in missing_scripts:
        display(Markdown(f"- `{path.relative_to(REPO_ROOT)}`"))
    raise FileNotFoundError("Pipeline scripts are missing. See TODO list above.")

display(Markdown("Pipeline scripts found:"))
for index, script_path in enumerate(PIPELINE_SCRIPTS, start=1):
    display(Markdown(f"{index}. `{script_path.relative_to(REPO_ROOT)}`"))

# %% [markdown]
# ## Run the existing pipeline
#
# By default, this notebook does not run the pipeline. Set `RUN_PIPELINE = True`
# in the next cell to execute the complete Article 2 workflow. When enabled, it
# calls the scripts above in order and stops immediately if any script fails.

# %%
RUN_PIPELINE = False

if not RUN_PIPELINE:
    display(Markdown("Pipeline execution is disabled because `RUN_PIPELINE = False`."))
    display(Markdown("To execute the workflow, set `RUN_PIPELINE = True` and rerun this cell."))
    display(Markdown("Scripts that would be executed:"))
    for index, script_path in enumerate(PIPELINE_SCRIPTS, start=1):
        display(Markdown(f"{index}. `{script_path.relative_to(REPO_ROOT)}`"))
else:
    for script_path in PIPELINE_SCRIPTS:
        relative_script = script_path.relative_to(REPO_ROOT)
        print(f"\n>>> Running {relative_script}")
        subprocess.run(
            [sys.executable, str(relative_script)],
            cwd=REPO_ROOT,
            check=True,
        )

    print("\nPipeline completed.")

# %% [markdown]
# ## Final article-ready tables
#
# The final tables are produced by the existing synthesis and article-ready
# output scripts. This section loads the generated table files if they exist and
# reports TODOs otherwise.

# %%
FINAL_TABLES = [
    REPO_ROOT / "tables/private/manuscript/table_enzyme_diagnostic_thresholds_private.md",
    REPO_ROOT / "tables/private/manuscript/table_enzyme_model_support_private.md",
    REPO_ROOT / "tables/private/manuscript/table_enzyme_diagnostic_thresholds_private.csv",
    REPO_ROOT / "tables/private/manuscript/table_enzyme_model_support_private.csv",
]

for table_path in FINAL_TABLES:
    if not table_path.exists():
        display(Markdown(f"### TODO: missing final table `{table_path.relative_to(REPO_ROOT)}`"))
        continue

    display(Markdown(f"## `{table_path.relative_to(REPO_ROOT)}`"))
    if table_path.suffix == ".md":
        display(Markdown(table_path.read_text(encoding="utf-8")))
    elif table_path.suffix == ".csv":
        display(pd.read_csv(table_path))

# %% [markdown]
# ## Final article-ready figures
#
# The vertical diagnostic panel is the main article-ready figure produced by the
# existing workflow. The PNG is displayed inline when available. The PDF path is
# listed for manuscript use.

# %%
FINAL_FIGURES = [
    REPO_ROOT / "figures/private/manuscript/figure_enzyme_diagnostic_panel_vertical.png",
    REPO_ROOT / "figures/private/manuscript/figure_enzyme_diagnostic_panel_vertical.pdf",
]

for figure_path in FINAL_FIGURES:
    if not figure_path.exists():
        display(Markdown(f"### TODO: missing final figure `{figure_path.relative_to(REPO_ROOT)}`"))
        continue

    display(Markdown(f"## `{figure_path.relative_to(REPO_ROOT)}`"))
    if figure_path.suffix.lower() == ".png":
        display(Image(filename=str(figure_path)))
    else:
        display(Markdown(f"Available at `{figure_path.relative_to(REPO_ROOT)}`"))

# %% [markdown]
# ## Reproducibility checklist
#
# This checklist records the files used by the workflow, the scripts executed,
# the final outputs expected by the manuscript, and the execution environment.

# %%
CHECKLIST_PATHS = {
    "private_inputs": PRIVATE_INPUTS,
    "expected_prepared_outputs": EXPECTED_PREPARED_OUTPUTS,
    "method_documents": DOCS,
    "pipeline_scripts": PIPELINE_SCRIPTS,
    "final_tables": FINAL_TABLES,
    "final_figures": FINAL_FIGURES,
}

display(Markdown("## Checklist"))

for group_name, paths in CHECKLIST_PATHS.items():
    display(Markdown(f"### {group_name}"))
    for path in paths:
        status = "found" if path.exists() else "TODO: missing"
        display(Markdown(f"- `{path.relative_to(REPO_ROOT)}`: {status}"))

display(Markdown("### environment"))
display(Markdown(f"- Python: `{sys.version.split()[0]}`"))
display(Markdown(f"- Platform: `{platform.platform()}`"))
display(Markdown(f"- Repository root: `{REPO_ROOT}`"))
display(Markdown(f"- Timestamp: `{datetime.now().isoformat(timespec='seconds')}`"))

# %% [markdown]
# ## Notes
#
# This notebook intentionally keeps the analysis logic inside the existing
# scripts. Any methodological change should be made in the appropriate script or
# document, then rerun here as an orchestrated workflow.
