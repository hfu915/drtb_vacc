# drtb_vacc

This folder contains data files and codes for the manuscript entitled **"Modelling the global burden of drug-resistant tuberculosis averted by a post-exposure vaccine"**, submitted to _Nature Communications_. This study aims to evaluate the impact of a future vaccine on rifampicin-resistant tuberculosis (TB) incidence and mortality in 30 countries. Analyses were performed and codes were tested using MATLAB R2016b. Software information and installation guide can be found at https://uk.mathworks.com/products/matlab.html. The installation process takes less than 30 minutes on a normal desktop computer. There is no non-standard hardware required for running the analysis.

We categorised the 30 countries into four groups, representing different TB epidemiology and care cascade (see Table 1 in the manuscript and Supplementary information). We ran the analysis for each country and then aggregated the results to obtain regional and global estimates. Here we provide the codes and materials for analysis in India.


## Description of files
### script
* `Setup_model_ct1.m`: input data and set up model structure
* `Get_vacc_results_ct1.m`: obtain projections of population dynamics by different vaccination scenarios, after setting up the model (`Setup_model_ct1.m`) and loading the posterior samples (`MCMC_posterior_IND.mat`)
* `Get_summary_results.m`: obtain country-level estimates of vaccine impacts from model projections (`Vacc_Tx_IND_ve50.mat`)

### data
* `dataHBC_dmg.xlsx`: country-specific demographics
* `dataHBC_tbb.xlsx`: country-specific TB burden and care cascade 
* `dataHBC_tbb_yrs.xlsx`: additional country-specific TB burden (multiple time points)
* `dataHBC_hiv.xlsx`: country-specific human immunodeficiency virus care cascade 
* `dataHBC_ctt.xlsx`: country-specific mixing matrix (adpated from Prem et al., PLoS Comput Biol 2017)
* `dataHBC_antibiotic_DHS.xlsx`: country-specific proportions of antibiotic use
* `dataHBC_cough_IHME.xlsx`: country-specific burden of cough-related diseases
* `dataHBC_vacCov.xlsx`: country-specific vaccination coverages

### output
* `MCMC_posterior_IND.mat`: posterior samples and compartment sizes at 2019 acquried from calibration
* `Vacc_Tx_IND_ve50.mat`: expected results of model projections for vaccine impacts in India, by processing the script `Get_vacc_results_ct1.m`

### funtion
* `get_addresses.m`: naming structure for different states
* `goveqs_basis_ct1.m`: governing equations for ODEs and time-varying transitions  
* `goveqs_scaleup_ct1.m`: governing equations for linear scale-up, built on the function `goveqs_basis_ct1.m`  
* `goveqs_scaleup_int_ct1.m`: governing equations for linear scale-up during the improvement of durg-resistance management, built on the function `goveqs_basis_ct1.m`
* `make_model_ct1.m`: model structure of TB transmission and care cascade
 

## Demo: Estimating vaccine impacts on drug-resistant TB in India

By specifying the country index (`icty = 1`), we can run a country-specific analysis for India. The script `Get_vacc_results_ct1.m` produces projections of population dynamics in different vaccination scenarios and then saves these results in `Vacc_Tx_IND_ve50.mat`. The simulation time takes around 6 minutes on a normal desktop computer. These results can be further processed to calculate numbers and proportions of vaccine-averted drug-resistant TB incidence and mortality, by executing the script `Get_summary_results.m`.

