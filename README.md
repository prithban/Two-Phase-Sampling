# Two-Phase-Sampling
Codes for Two-Phase Longitudinal Sampling Design

The codes for simulation studies, additional simulations, and the real data analysis is included in this repository. The description of the files are as follows:

1) simulation_study : Simulation codes for all 6 scenarios considered in the manuscript.

2) new_simulation_study : Additional simulations conducted as per reviewers request with high sensitivity and specificity.

3) data_analysis_hclustering: Simulated home healthcare study via hierarchical clustering.

4) data_analysis_mclustering: Simulated home healthcare study via model-based clustering.

5) Sensitivity_Specificity_check_hclustering: Sensitivity/Specificity analysis of simulated home healtcare study via hierarchical clustering.

6) Sensitivity_Specificity_check_mclustering: Sensitivity/Specificity analysis of simulated home healtcare study via model-based clustering.

P.S. The exact real data cannot be included in the github due to restriction that has been agreed upon with the PI of the original grant producing the data. Note, the Home Healthcare data is used in our paper for the demonstration purpose of the proposed method. So in lieu of that we took the real data and added small amount white noise keeping it's mean and s.d. unchanged for the continuous variables. The binary variables are sampled from a similar binomial random variable. This perturbations have very slight effect on the estimated incidence and prevalence rate and results are more or less reproducible, with difference appearing in the 2nd decimal place or after. 
