Data workflow
Project Author: Natalie Dupont
Code Authors: Natalie Dupont, Eric Pedersen, Ariella Fuzaylov
Code also sourced from Uri Simonsohn (https://osf.io/wdbmr)
Pedersen Quantitative Fisheries Ecology Lab
Supervisor: Dr. Eric Pedersen

This project is associated with "Applying A New Approach for Quantifying Peaks to the Intermediate Disturbance Hypothesis". It tests the performance of a shape-constrained additive model approach to detecting the presence of a peak in a dataset, and a generalized additive model first-derivative approach to estimating the location of a peak. These two methods are tested using simulation data and real IDH data from figures in published literature

This file describes data workflow through the scripts
Script authorship and descriptions are detailed within the respective scripts
-----------------------------------------------------------------------------------------------------
Simulation data

1. simulations.R: simulation data is generated and exported to ~data/simulation as "CLEAN_simulationData.csv" 
2. tests_sims.R: Simulation data ("CLEAN_simulationData.csv") is imported into to run peak detection and quantification tests. 
2. A. parametrized_peak_finder.R: contains functions for peak quantification 
2. B. testing_maxima.R: contains functions used by parametrized_peak_finder.R
Test results are exported to data/simulation
3. sim_results_analysis.R: Test results for detection tests, and quantification tests are analyzed and summarized to produce final statistics and figures

-----------------------------------------------------------------------------------------------------
Real data

1. data is procured from literature raw data table or digitized from a figure using web plot Digitizer (https://apps.automeris.io/wpd/)
2. data is saved as a csv in the ~data folder under the naming convention "AuthorLastNamesDate_DataSource.csv"
3. data_cleaning.R: Data is cleaned, transformed and cleaned data is output to the ~data/processed folder with the naming convention "CLEAN_AuthorLastNamesDate_DataSource.csv"
4. real_data_analysis.R: CLEAN csv files are imported and analyzed using functions from tests_sims.R and parametrized_peak_finder.R