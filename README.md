# CAR T-cell Timing and Outcomes

This repository contains the analysis code for CAR T-cell timing and patient outcomes project.

## Repository Structure

The project consists of five main scripts:

### `01_tables.R`
- Prepares baseline tables and descriptive summaries for the study cohort.
- Includes summary statistics for demographics, clinical variables, and treatment characteristics.
- Includes outcomes tables.

### `02_logistic.R`
- Runs the main logistic regression models for primary and secondary binary outcomes.
- Includes model-based marginal effect estimates.
- Includes sensitivity analyses to test the robustness of findings.

### `03_survival.R`
- Performs KM analysis for overall and event-free survival.
- Calculates adjusted restricted mean survival time (RMST) at 2 years.

### `04_subgroup.R`
- Conducts subgroup analyses for logistic and survival outcomes.
- Stratifies by sex, cancer type, hospital, and other prespecified factors.
- Includes interaction terms and marginal effect plots for interpretation.

### `05_labs.R`
- Analyzes key laboratory values in relation to timing of CAR T-cell infusion.
- Compares distributions and trends across early and late infusion groups.

## Dependencies

The analysis was developed using R (≥ 4.0) and the following key packages:

- `tidytable`
- `collapse`
- `lme4`, `lmerTest`
- `broom.mixed`
- `survival`, `rms`
- `pseudo`
- `boot`
- `gtsummary`, `gt`
- `ggplot2`

