# Cognitive Representations of Self-Harm as Predictors of Coping and Attitudes Towards Service Engagement

- This repository contains the archived data and sensitivity analysis code for Kennett, J., Friedrich, C., Rahman, R. (2026) Cognitive Representations of Self-Harm as Predictors of Coping and Attitudes Towards Service Engagement.
- For a detailed report of the methods and results, see the upcoming publication
- This project is licensed under the terms of the Creative Commons Attribution 4.0 International license (CC-BY 4.0) (<https://creativecommons.org/licenses/by/4.0/>).

## Contents:

- **data/**

  - s1_data.csv - raw and anonymised data for the study

  - data_dictionary.xlsx - information on variables in s1_data.csv

  - imputed_data/ - all 40 imputed datasets

- **util/**

  - imputation.R - code for the sensitivity analysis described in the publication. This can be executed when s1_data.csv is also downloaded.

- **imputation_diagnostics/**

  - convergence_plot.png - trace plots of multiple imputation algorithm convergence for the mean and standard deviation of items Q12_27 and Q12_28

  - distribution_Q12_27.csv - observed distribution and mean/min/max imputed distribution for item Q12_27

  - distribution_Q12_28.csv - observed distribution and mean/min/max imputed distribution for item Q12_28

  - Q12_27_imputed_values.csv - values for item Q12_27 for all imputed data sets

  - Q12_28_imputed_values.csv - values for item Q12_28 for all imputed data sets

  - strip_plot_Q12_27.png - scatter plot of real and imputed values for item Q12_27

  - strip_plot_Q12_28.png - scatter plot of real and imputed values for item Q12_28

- **imputation_results/**

  - cor_imputation_results.csv - observed correlation coefficient and mean, standard deviation, 2.5th and 97.5th percentiles across imputed data sets for every variable in the primary analysis

  - pooled_imputation_results.csv - pooled estimates, standard error, 95% confidence intervals, and p-values for each regression model using imputed data

Note: in imputation.R, to avoid constant re-generation of .csv and .png files in the diagnostics and results directories, code that creates these has been suppressed with a '\#'.
