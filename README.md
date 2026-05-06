# High-Speed Rail Well-Being PLS-SEM Analysis in R

This repository presents a Partial Least Squares Structural Equation Modelling (PLS-SEM) analysis examining how high-speed rail (HSR) service quality influences passengers’ subjective well-being and customer loyalty

The study models the relationships between:

- Service Quality
- Travel Experience Affect
- Cognitive Evaluation
- Subjective Well-Being
- Customer Loyalty
- Willingness to Pay

The analysis was conducted in R using the `seminr` package and includes:

- PLS-SEM model estimation
- Reliability and validity assessment
- HTMT analysis
- Structural path analysis
- Bootstrapping
- Direct, indirect, and total effects
- R² and f² evaluation
- HTML reporting

## Repository Structure

```text
data/       -> Survey dataset
scripts/    -> R scripts and R Markdown report
outputs/    -> Exported CSV results
reports/    -> Generated HTML report
```

## Main Files

- `scripts/hsr_plssem_analysis.R`  
  Main analysis script

- `scripts/analysis_report.Rmd`  
  R Markdown report for reproducible analysis

- `reports/analysis_report.html`  
  Final rendered analysis report

## Methodology

The study applies PLS-SEM to analyse behavioural and experiential relationships in high-speed rail travel, focusing on the mediating role of subjective well-being between travel experience and behavioural outcomes.

## Packages Used

- seminr
- knitr
- kableExtra
- DiagrammeR

## Author

Alireza MOradpour 
MSc Transportation Engineering and Mobility
