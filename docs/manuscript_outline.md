# Manuscript outline and working brief

## Working title

**Soil enzyme diagnostics for irrigated Palmer mango orchards in the Brazilian semiarid region**

Alternative titles:

1. **Diagnostic thresholds for soil enzyme activity in irrigated Palmer mango orchards**
2. **Soil biological indicators associated with relative yield in semiarid mango orchards**
3. **β-glucosidase, arylsulfatase, and GMea as biological indicators of orchard performance in irrigated mango systems**

## Project context

This manuscript is part of a broader research effort on soil quality and soil health assessment in irrigated Palmer mango orchards cultivated under semiarid conditions in Brazil.

A previous manuscript focused on the construction and validation of an integrated Soil Quality Index (SQI). The present manuscript has a narrower biological focus: soil enzyme activity as a diagnostic component of orchard soil functioning.

The analysis uses private field data from irrigated Palmer mango orchards. Public repository files should use synthetic or anonymized data only. Original data and private outputs should remain outside version control until publication and coauthor approval.

## Central message

Soil enzyme activity, especially GMea and β-glucosidase, was consistently associated with relative mango yield and provided useful dataset-specific diagnostic thresholds for identifying biologically constrained orchards.

The strongest interpretation is not that enzyme activity independently controls yield by itself, but that enzyme activity integrates relevant biological, chemical, physical, and management-related conditions of the orchard system.

## Main research question

Can soil enzyme activity be used to define diagnostic biological thresholds and high-yield-associated reference ranges for irrigated Palmer mango orchards in the Brazilian semiarid region?

## Specific objectives

1. Evaluate the association between soil enzyme indicators and relative mango yield.
2. Identify diagnostic enzyme thresholds associated with adequate relative yield.
3. Compare the diagnostic performance of β-glucosidase, arylsulfatase, and GMea.
4. Assess whether enzyme-yield relationships are purely linear or show threshold/plateau-like behavior.
5. Evaluate the stability of diagnostic thresholds using bootstrap resampling.
6. Examine whether enzyme signals remain meaningful after considering soil context, orchard structure, and farm effects.
7. Propose dataset-specific biological reference ranges associated with adequate orchard performance.

## Hypotheses

1. Low soil enzyme activity is associated with lower relative mango yield.
2. GMea provides the best integrated biological diagnostic signal because it combines β-glucosidase and arylsulfatase activities.
3. β-glucosidase provides a more robust individual enzyme signal than arylsulfatase.
4. Arylsulfatase is biologically relevant but more dependent on organic matter, soil fertility, structure, and farm context.
5. Enzyme-yield relationships are not adequately represented as purely linear across the full range of enzyme activity.
6. Diagnostic thresholds are more stable and agronomically defensible than plateau breakpoints.
7. High-yield-associated enzyme ranges should be interpreted as dataset-specific reference ranges, not universal optimum values.

## Core variables

### Main response variable

* `Prod_rel_pct`: relative yield per plant (%)

This is the main response variable for diagnostic threshold screening.

### Secondary response variable

* `Prod_rel_ha_pct`: relative yield per hectare (%)

This can be used as a sensitivity response if needed.

### Main enzyme indicators

* `Beta_glic`: β-glucosidase activity
* `Arilsulf`: arylsulfatase activity
* `GMea`: geometric mean of β-glucosidase and arylsulfatase activity

In this manuscript, GMea specifically means:

**GMea = sqrt(β-glucosidase × arylsulfatase)**

### Soil and orchard context variables

Important context variables include:

* soil organic matter (`MO_g_dm3`)
* organic carbon (`C_org_g_d3`)
* clay content (`Argila_g_kg`)
* bulk density (`Ds_g_cm3`)
* clay flocculation (`Floculacao_pct`)
* clay dispersion (`Dispercao_pct`)
* pH
* electrical conductivity (`CE_dS_m`)
* exchangeable bases
* exchangeable sodium percentage (`PST`)
* orchard age (`Idade_anos`)
* planting density (`Densidade_pl_ha`)
* farm identity (`Fazenda`)

## Diagnostic criterion

A relative yield cutoff of **70%** is used to classify orchards into adequate and low relative yield groups.

This cutoff is used for Cate–Nelson threshold screening and for defining the high-yield-associated group used to calculate reference ranges.

The 70% cutoff should be described as an operational diagnostic criterion, not as a universal biological boundary.

## Main methodological decisions

### Primary diagnostic threshold method

Cate–Nelson analysis is used as the main method for estimating diagnostic enzyme thresholds.

Rationale:

* It has a direct agronomic interpretation.
* It identifies a practical threshold separating low-yielding and adequate-yielding orchards.
* It is more appropriate for diagnostic threshold definition than purely linear correlation.
* It was more stable than plateau breakpoints in bootstrap resampling.

### Complementary sensitivity methods

The following methods are used as complementary analyses:

1. **Spearman correlation**

   * Used to assess monotonic association between enzymes and relative yield.

2. **Shallow recursive partitioning**

   * Used as a complementary data-driven split.
   * Not treated as the main diagnostic threshold method.

3. **Linear, quadratic, linear-plateau, and quadratic-plateau models**

   * Used to assess whether enzyme-yield relationships are purely linear or better described by nonlinear/plateau-like responses.
   * Plateau breakpoints are interpreted cautiously because their bootstrap intervals are wide.

4. **Bootstrap resampling**

   * Used to assess stability of correlations, thresholds, high-yield-associated ranges, and plateau breakpoints.

5. **Context-adjusted models**

   * Used to assess whether enzyme indicators retain signal after adjustment for soil context, orchard structure, and farm effects.

## Main results to emphasize

### 1. General enzyme-yield associations

All three enzyme indicators were positively associated with relative yield per plant.

Main Spearman correlations:

* GMea: ρ = 0.608, p < 0.001
* β-glucosidase: ρ = 0.565, p < 0.001
* Arylsulfatase: ρ = 0.559, p < 0.001

GMea had the strongest overall association with relative yield.

### 2. Diagnostic thresholds

Cate–Nelson diagnostic thresholds using the 70% relative yield criterion were:

* GMea: 87.0
* β-glucosidase: 81.4
* Arylsulfatase: 91.6

These values should be interpreted as minimum diagnostic levels associated with lower risk of poor relative yield.

They should not be interpreted as universal optimum values.

### 3. Bootstrap support

Bootstrap results suggested that Cate–Nelson thresholds were more stable than plateau breakpoints.

Bootstrap median and 95% interval for Cate–Nelson thresholds:

* GMea: 96.1 (77.3–132.8)
* β-glucosidase: 81.4 (79.1–120.6)
* Arylsulfatase: 91.6 (82.9–151.4)

β-glucosidase showed the most stable diagnostic threshold.

GMea showed strong overall diagnostic value but a broader bootstrap interval.

Arylsulfatase showed useful biological signal but weaker stability and less independence as an isolated diagnostic indicator.

### 4. High-yield-associated reference ranges

The interquartile range of enzyme values among orchards with relative yield per plant ≥70% was used as a descriptive high-yield-associated reference range.

Observed high-yield-associated IQRs:

* GMea: 114.0–144.0
* β-glucosidase: 105.4–144.4
* Arylsulfatase: 113.2–162.1

These ranges should be described as dataset-specific reference ranges associated with adequate relative yield, not as universal optimal ranges.

### 5. Nonlinear and plateau-like behavior

Model comparison indicated that enzyme-yield relationships were not purely linear.

Best or competitive response models included linear-plateau, quadratic, and quadratic-plateau forms.

For GMea:

* best model: linear-plateau
* R² = 0.432
* RMSE = 16.31
* estimated breakpoint: 140.8
* competitive models also included quadratic and quadratic-plateau forms

For β-glucosidase:

* best model: quadratic
* R² = 0.398
* RMSE = 16.80
* estimated vertex: 142.3
* linear-plateau was also competitive

For arylsulfatase:

* best model: linear-plateau
* R² = 0.344
* RMSE = 17.53
* estimated breakpoint: 140.2
* quadratic and quadratic-plateau were also competitive

However, plateau breakpoints had wide bootstrap intervals and should not be treated as primary diagnostic thresholds.

Suggested interpretation:

Low enzyme activity was associated with lower relative yield, whereas increases above an intermediate level did not necessarily translate into proportional yield gains. This supports a threshold/plateau-like interpretation, where enzyme activity becomes less limiting and other soil, orchard, or management factors may control yield.

### 6. Soil context

Enzyme indicators were strongly associated with soil organic matter and organic carbon.

The strongest enzyme-context relationships included:

* arylsulfatase with soil organic matter
* GMea with soil organic matter
* β-glucosidase with soil organic matter

This supports the interpretation that enzyme activity reflects an integrated biological and edaphic condition.

### 7. Orchard structure and farm effects

Orchard age and planting density were strongly associated with relative yield.

Farm effects absorbed substantial variation in enzyme-yield relationships, indicating that part of the enzyme signal is structured at the orchard/farm-system level.

When only farm fixed effects were considered, enzyme effects became weak. However, after adjusting for orchard structure, GMea and β-glucosidase retained useful signal.

This supports the interpretation that GMea and β-glucosidase are more consistent than arylsulfatase, while enzyme activity as a whole should be interpreted as an integrated system-level indicator rather than a fully independent causal driver.

## Proposed hierarchy of indicators

### GMea

Role: **main integrated biological indicator**

Justification:

* strongest Spearman association with relative yield
* clear diagnostic threshold
* good overall diagnostic performance
* supported by nonlinear/plateau-like model behavior
* integrates the behavior of both β-glucosidase and arylsulfatase

Interpretation:

GMea is the best overall indicator of soil biological functioning among the tested enzyme variables. It should be used as the main biological diagnostic indicator.

### β-glucosidase

Role: **most robust individual enzyme indicator**

Justification:

* consistent positive association with yield
* stable diagnostic threshold
* maintained useful signal after adjustment for orchard structure and soil context
* suitable as a minimum biological-level indicator

Interpretation:

β-glucosidase is the most robust individual enzyme indicator and may be especially useful for identifying orchards below a minimum biological functioning level.

### Arylsulfatase

Role: **complementary sulfur-cycling enzyme indicator**

Justification:

* positive association with relative yield
* biologically relevant
* strongly associated with organic matter, bases, and soil structure
* weaker independent signal after adjustment

Interpretation:

Arylsulfatase should be interpreted as a complementary enzyme indicator. It contributes biological information but is less reliable as a standalone diagnostic indicator than GMea or β-glucosidase.

## Candidate main figure

**Figure 1. Soil enzyme diagnostic thresholds and high-yield-associated reference ranges.**

Candidate file:

* `figures/private/manuscript/figure_enzyme_diagnostic_panel_vertical.png`
* `figures/private/manuscript/figure_enzyme_diagnostic_panel_vertical.pdf`

Figure content:

* Panel A: β-glucosidase
* Panel B: arylsulfatase
* Panel C: GMea
* y-axis: relative yield per plant (%)
* x-axis: enzyme activity
* red dashed vertical line: Cate–Nelson diagnostic threshold
* blue dotted horizontal line: 70% relative yield cutoff
* shaded band: IQR of enzyme values among orchards with relative yield ≥70%

Interpretation note:

The shaded band represents a dataset-specific high-yield-associated range and should not be interpreted as a universal optimum range.

Potential figure caption:

**Figure 1. Diagnostic thresholds and high-yield-associated reference ranges for soil enzyme activity in irrigated Palmer mango orchards. Points represent individual observations. The horizontal dotted line indicates the 70% relative yield cutoff used to define adequate relative performance. The vertical dashed line indicates the Cate–Nelson diagnostic threshold for each enzyme indicator. The shaded band represents the interquartile range of enzyme values observed among orchards with relative yield per plant ≥70%. This band represents a dataset-specific high-yield-associated range, not a universal optimum range.**

## Candidate main table

**Table 1. Diagnostic thresholds and high-yield-associated reference ranges for soil enzyme indicators.**

Candidate file:

* `tables/private/manuscript/table_enzyme_diagnostic_thresholds_private.csv`
* `tables/private/manuscript/table_enzyme_diagnostic_thresholds_private.md`

Suggested table columns:

* enzyme indicator
* role
* Spearman correlation
* Cate–Nelson threshold
* bootstrap median and 95% interval
* high-yield-associated IQR
* balanced accuracy
* interpretation

## Candidate supplementary table

**Supplementary Table S1. Model support for nonlinear and context-adjusted enzyme-yield relationships.**

Candidate file:

* `tables/private/manuscript/table_enzyme_model_support_private.csv`
* `tables/private/manuscript/table_enzyme_model_support_private.md`

Purpose:

This table supports the interpretation that enzyme-yield relationships are not purely linear and that enzyme signals are partly integrated with orchard structure, soil context, and farm-level effects.

## Suggested manuscript structure

### Abstract

The abstract should emphasize:

1. The need for biological diagnostic indicators in irrigated mango orchards.
2. The use of β-glucosidase, arylsulfatase, and GMea as enzyme indicators.
3. The association with relative yield.
4. The use of Cate–Nelson thresholds and high-yield-associated IQRs.
5. The main finding: GMea and β-glucosidase were the most useful indicators, while arylsulfatase was complementary.
6. The caution: ranges are dataset-specific and require external validation.

### Introduction

The introduction should cover:

1. Importance of mango production in irrigated semiarid systems.
2. Soil biological functioning as a component of orchard sustainability.
3. Enzyme activity as an early and sensitive indicator of soil biological processes.
4. β-glucosidase as an indicator of carbon cycling.
5. Arylsulfatase as an indicator of sulfur cycling.
6. GMea as an integrated enzyme index.
7. Need for diagnostic thresholds or reference ranges rather than only descriptive enzyme values.
8. Gap: lack of enzyme diagnostic ranges for irrigated mango orchards in semiarid conditions.
9. Objective and hypotheses.

### Materials and methods

Suggested subsections:

1. Study area and orchard selection
2. Soil sampling and laboratory analyses
3. Enzyme activity measurements
4. Relative yield calculation
5. GMea calculation
6. Diagnostic threshold estimation
7. Model comparison
8. Bootstrap stability analysis
9. Soil context, orchard structure, and farm-effect analyses
10. Software and reproducibility

Important wording:

* Use “diagnostic threshold” instead of “optimal threshold”.
* Use “high-yield-associated reference range” instead of “optimal range”.
* Use “dataset-specific” where appropriate.
* State that plateau breakpoints were used as interpretive support, not as primary diagnostic thresholds.

### Results

Suggested results flow:

#### 1. Enzyme activity and relative yield

Report Spearman correlations and identify GMea as the strongest integrated indicator.

#### 2. Diagnostic thresholds

Report Cate–Nelson thresholds and balanced accuracy.

Emphasize:

* GMea threshold: 87.0
* β-glucosidase threshold: 81.4
* arylsulfatase threshold: 91.6

#### 3. High-yield-associated ranges

Report IQRs for the ≥70% relative yield group.

Emphasize that these are descriptive reference ranges.

#### 4. Nonlinear response and plateau-like behavior

Report model comparison results.

Emphasize that nonlinear and plateau-like models were competitive but plateau breakpoints were less stable than Cate–Nelson thresholds.

#### 5. Soil and orchard context

Report that enzyme activity was strongly associated with organic matter, organic carbon, exchangeable bases, structure, orchard age, density, and farm context.

#### 6. Indicator hierarchy

State that:

* GMea is the main integrated biological indicator.
* β-glucosidase is the most robust individual enzyme indicator.
* arylsulfatase is complementary.

### Discussion

Suggested discussion points:

#### 1. Why GMea performed well

GMea integrates carbon and sulfur cycling enzyme activity and may better represent overall soil biological functioning than either enzyme alone.

#### 2. Why β-glucosidase was robust

β-glucosidase is directly linked to carbon cycling and may respond strongly to biologically active organic matter and root-zone microbial activity.

#### 3. Why arylsulfatase was complementary

Arylsulfatase is relevant to sulfur cycling but appears more dependent on soil organic matter, exchangeable bases, structure, and farm-level context. Its diagnostic value is stronger as part of an integrated biological assessment than as an isolated indicator.

#### 4. Threshold versus optimum

The Cate–Nelson thresholds should be interpreted as minimum diagnostic levels, not optimum values.

The high-yield-associated IQRs should be interpreted as observed reference ranges within this dataset, not as universal recommendations.

#### 5. Nonlinear response

The enzyme-yield relationship resembles a threshold/plateau response rather than a strictly linear response. This is agronomically reasonable because once minimum biological functioning is reached, other factors such as salinity, sodicity, soil structure, fertility balance, orchard age, planting density, irrigation, or management may become limiting.

#### 6. System-level interpretation

Enzyme activity should be interpreted as an integrated biological-edaphic signal. The goal is not to claim that enzyme activity alone causes yield variation, but that enzyme activity helps diagnose orchard biological condition associated with yield performance.

#### 7. Limitations

Important limitations:

* dataset-specific thresholds
* limited number of high-yield observations
* farm-level structure
* no external validation yet
* possible confounding by orchard age, density, and management
* enzyme activity is an indicator, not a direct management prescription

#### 8. Practical implications

Potential practical use:

* identify orchards with low biological functioning
* support soil health monitoring in irrigated mango systems
* complement chemical and physical soil indicators
* guide further investigation of organic matter management, structure, sodicity, and biological functioning

### Conclusions

Suggested conclusion points:

1. Soil enzyme activity was consistently associated with relative yield in irrigated Palmer mango orchards.
2. GMea was the best integrated biological indicator.
3. β-glucosidase was the most robust individual enzyme indicator.
4. Arylsulfatase was useful as a complementary sulfur-cycling indicator but less independent as a standalone diagnostic variable.
5. Cate–Nelson thresholds provided practical diagnostic values associated with adequate relative yield.
6. High-yield-associated IQRs provided dataset-specific reference ranges.
7. Enzyme-yield relationships were not purely linear, supporting a threshold/plateau-like interpretation.
8. Further validation is required before using these values as general recommendations.

## Terms to use

Preferred terms:

* diagnostic threshold
* minimum biological level
* high-yield-associated range
* dataset-specific reference range
* biological indicator
* integrated biological signal
* soil biological functioning
* threshold/plateau-like response
* adequate relative yield group

## Terms to avoid or use cautiously

Avoid or use cautiously:

* optimal range
* ideal range
* universal threshold
* recommendation range
* enzyme sufficiency range
* enzyme causes yield increase
* high enzyme activity reduces yield

If “optimal” is used, it should be qualified as dataset-specific and preliminary.

## Suggested short narrative for the manuscript

Soil enzyme activity was positively associated with relative yield in irrigated Palmer mango orchards, but the relationships were not adequately described as purely linear responses. GMea showed the strongest integrated diagnostic behavior, while β-glucosidase was the most robust individual enzyme indicator. Arylsulfatase contributed complementary information related to sulfur cycling but was more dependent on organic matter, soil structure, exchangeable bases, and farm context. Cate–Nelson thresholds provided practical minimum biological levels associated with adequate relative yield, whereas the interquartile range of enzyme values among orchards with relative yield ≥70% provided dataset-specific high-yield-associated reference ranges. These values should be interpreted as preliminary diagnostic references requiring external validation, not as universal optimum ranges.

## Remaining tasks before manuscript drafting

1. Decide whether to include relative yield per hectare as a sensitivity result.
2. Decide whether to include clay-normalized enzyme indicators in supplementary material only.
3. Create final manuscript-ready figure caption.
4. Create final manuscript-ready table caption.
5. Prepare a methods paragraph for Cate–Nelson analysis.
6. Prepare a methods paragraph for bootstrap resampling.
7. Prepare a methods paragraph for nonlinear/plateau model comparison.
8. Prepare a limitations paragraph.
9. Decide target journal and formatting requirements.
10. Create the initial LaTeX manuscript structure.

## Current article-ready private outputs

Candidate main figure:

* `figures/private/manuscript/figure_enzyme_diagnostic_panel_vertical.png`
* `figures/private/manuscript/figure_enzyme_diagnostic_panel_vertical.pdf`

Candidate supplementary/presentation figure:

* `figures/private/manuscript/figure_enzyme_diagnostic_panel_horizontal.png`
* `figures/private/manuscript/figure_enzyme_diagnostic_panel_horizontal.pdf`

Candidate main table:

* `tables/private/manuscript/table_enzyme_diagnostic_thresholds_private.csv`
* `tables/private/manuscript/table_enzyme_diagnostic_thresholds_private.md`

Candidate supplementary model-support table:

* `tables/private/manuscript/table_enzyme_model_support_private.csv`
* `tables/private/manuscript/table_enzyme_model_support_private.md`

## Working interpretation to preserve

The manuscript should not overclaim causality. The safest and strongest interpretation is that enzyme activity acts as an integrated biological indicator of orchard soil condition. Low enzyme activity identifies biologically constrained systems, while high enzyme activity does not guarantee high yield because other factors may become limiting. GMea and β-glucosidase are the most defensible indicators for diagnostic use, while arylsulfatase is best interpreted as complementary.
