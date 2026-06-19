# Mango Soil Enzyme Diagnostics

This repository documents a reproducible analytical workflow for soil enzyme diagnostics in irrigated Palmer mango orchards cultivated under semiarid conditions in Brazil.

The project supports a manuscript focused on soil biological indicators, diagnostic enzyme thresholds, and high-yield-associated reference ranges derived from the relationship between enzyme activity and relative orchard yield.

## Research focus

This workflow evaluates three main soil biological indicators:

* β-glucosidase activity;
* arylsulfatase activity;
* GMea, defined as the geometric mean of β-glucosidase and arylsulfatase activity.

The main analytical goal is to assess whether these indicators can support biological diagnostic interpretation in irrigated mango orchards, with emphasis on:

* enzyme-yield associations;
* Cate–Nelson diagnostic thresholds;
* high-yield-associated reference ranges;
* nonlinear or threshold/plateau-like response behavior;
* bootstrap stability;
* soil, orchard-structure, and farm-context sensitivity.

## Manuscript context

This repository is associated with the manuscript draft:

**Soil Enzyme Diagnostics for Irrigated Palmer Mango Orchards in the Brazilian Semiarid Region**

The current working interpretation is that:

* GMea is the main integrated biological indicator;
* β-glucosidase is the most robust individual enzyme indicator;
* arylsulfatase is a complementary sulfur-cycling indicator;
* proposed thresholds are dataset-specific diagnostic references, not universal optimum values.

## Data availability

The original field dataset is private and is not included in this repository.

The private dataset contains soil biological, chemical, physical, orchard-structure, and productivity variables from irrigated Palmer mango orchards. It is excluded from version control because manuscripts derived from the dataset are still in preparation and data sharing requires coauthor and partner approval.

Public examples will use synthetic or anonymized data products with the same structure as the private analytical dataset. These files will be used only for workflow demonstration and reproducibility of the analytical structure.

## Repository structure

```text
mango-soil-enzyme-diagnostics/
├── data/
│   └── synthetic/              # Future synthetic example datasets
├── docs/                       # Project scope, analysis strategy, and methodological notes
├── figures/
│   └── synthetic/              # Future public example figures
├── manuscript/                 # Manuscript-related placeholder files
├── notebooks/                  # Exploratory notebooks, if needed
├── scripts/                    # Reproducible analysis scripts
├── tables/
│   └── synthetic/              # Future public example tables
├── .gitignore
├── README.md
└── requirements.txt
```

Private analytical inputs and outputs are intentionally excluded from GitHub, including:

```text
data/private/
data/processed/private/
figures/private/
tables/private/
```

## Analytical workflow

The private workflow currently includes scripts for:

1. inspecting the private dataset structure;
2. preparing the enzyme diagnostics dataset;
3. exploring enzyme-yield relationships;
4. screening diagnostic thresholds;
5. summarizing threshold stability;
6. testing sensitivity to the experimental group;
7. creating main diagnostic summary tables;
8. generating diagnostic threshold figures;
9. comparing linear, quadratic, and plateau-type response models;
10. evaluating bootstrap stability;
11. analyzing enzyme relationships with soil context;
12. testing enzyme signals after orchard-structure adjustment;
13. evaluating farm-effect sensitivity;
14. synthesizing article-ready enzyme diagnostic results;
15. organizing private article-ready tables and figures.

The public version of the repository documents the workflow and code structure. Full reproducibility with real data is restricted until data release is approved.

## Main methods

The workflow uses:

* Spearman correlation for monotonic enzyme-yield associations;
* Cate–Nelson analysis for primary diagnostic threshold estimation;
* balanced accuracy to summarize diagnostic classification performance;
* interquartile ranges of the adequate-yield group to define high-yield-associated reference ranges;
* linear, quadratic, linear-plateau, and quadratic-plateau models to assess response shape;
* bootstrap resampling to evaluate threshold and model stability;
* context-adjusted models and farm-effect sensitivity analyses to avoid overinterpreting enzyme activity as an isolated causal driver.

## Interpretation principles

The repository follows a cautious diagnostic interpretation:

* enzyme activity is treated as an integrated biological-edaphic indicator;
* low enzyme activity may identify biologically constrained orchard conditions;
* high enzyme activity does not guarantee high yield;
* Cate–Nelson thresholds are treated as primary diagnostic thresholds;
* plateau breakpoints are treated as supporting evidence for nonlinear behavior;
* high-yield-associated ranges are dataset-specific reference ranges, not recommendation ranges.

## Software

The analyses were developed in Python 3.12.3.

Main Python packages:

* pandas;
* NumPy;
* SciPy;
* statsmodels;
* scikit-learn;
* Matplotlib.

Install the basic requirements with:

```bash
pip install -r requirements.txt
```

## Current status

The private analytical workflow has been implemented and used to generate manuscript-oriented diagnostic tables and figures. The public repository currently documents the workflow structure, methodological decisions, and reproducible scripts.

Next public-facing development steps include:

* adding a synthetic example dataset;
* generating synthetic demonstration tables and figures;
* adding a short notebook showing how the workflow runs with public synthetic data;
* improving package version pinning for reproducibility.

## License and citation

A formal license and citation instructions will be added after manuscript and data-sharing decisions are finalized.
