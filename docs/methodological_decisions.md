# Methodological decisions

This document records the main methodological decisions adopted for the private exploratory workflow of the `mango-soil-enzyme-diagnostics` project.

## Scope

This analysis evaluates soil enzyme activity as a diagnostic component of soil biological functioning in irrigated Palmer mango orchards from the Brazilian semiarid region. The main indicators are β-glucosidase, arylsulfatase, and GMea.

In this project, GMea refers specifically to the geometric mean of β-glucosidase and arylsulfatase activities.

## Response variable

The main response variable is relative yield per plant. Relative yield per area is considered a secondary or sensitivity response when needed.

A relative yield cutoff of 70% is used to define the group of orchards with adequate relative performance for diagnostic threshold screening.

## Main diagnostic threshold method

Cate–Nelson analysis is used as the main method for estimating diagnostic enzyme thresholds.

This choice is based on its agronomic interpretation: it identifies a practical threshold that best separates low-yielding and adequate-yielding orchards according to the 70% relative yield criterion.

The diagnostic threshold should be interpreted as a minimum biological level associated with lower risk of poor relative yield, not as an optimum value.

## Sensitivity analyses

Shallow recursive partitioning is used as a complementary sensitivity method. It provides an alternative data-driven split in relative yield but is not treated as the primary diagnostic threshold method.

Linear, quadratic, linear-plateau, and quadratic-plateau models are compared to evaluate whether enzyme-yield relationships are adequately represented as purely linear responses.

Bootstrap resampling is used to assess the stability of correlations, diagnostic thresholds, high-yield reference ranges, and plateau breakpoints.

## Interpretation of plateau models

Plateau-type models are used to support the interpretation that enzyme activity may become less limiting beyond an intermediate level.

However, plateau breakpoints are not used as primary diagnostic thresholds unless they show acceptable stability across bootstrap resampling.

When plateau breakpoints are unstable, they are interpreted qualitatively rather than as definitive critical values.

## High-yield-associated reference range

The interquartile range of enzyme values observed among orchards with relative yield per plant ≥70% is used as a descriptive high-yield-associated reference range.

This range is not interpreted as a universal optimum range. It represents the central distribution of enzyme values observed in the adequate-yield group within this dataset.

## Biological interpretation

Low enzyme activity is interpreted as a potential indicator of biological limitation.

High enzyme activity does not necessarily imply high yield, because other soil, plant, or management factors may become limiting once minimum biological functioning is reached.

Therefore, the expected response is not necessarily monotonic. A threshold or plateau-like interpretation is considered more agronomically appropriate than assuming proportional yield increases across the full range of enzyme activity.

## Figure interpretation

In diagnostic figures:

- points represent individual observations;
- the horizontal line represents the 70% relative yield cutoff;
- the vertical line represents the Cate–Nelson diagnostic threshold;
- the shaded band represents the interquartile range of enzyme values among orchards with relative yield per plant ≥70%.

These graphical elements are intended to distinguish minimum diagnostic thresholds from descriptive high-yield-associated ranges.
