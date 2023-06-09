# Transformed-linear Prediction for Extremes

# Author Contributions Checklist Form

## Data

### Abstract
This manuscript analyzes daily levels of $\text{NO}_2$ measured in parts per billion (ppb) by the US Environmental Protection Agency (EPA).
The air pollution dataset contains geographical information for gauging stations in Washington DC.
Additionally, we analyze the daily value-averaged returns of 30 US industry portfolios from 1950 to 2020.

### Availability
To access the air pollution data, users can visit the EPA website and download the daily data at (https://www.epa.gov/outdoor-air-quality-data/download-daily-data). Similarly, the industry portfolio data are publicly available at the Kenneth French Data Library (https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data\_library.html).

### Description
The preprocessing of both datasets is described in Section 6 of the manuscript, and the preprocessed data can be found in the 'Data' folder. The simulated data are also available in the 'Output' folder or can be generated using the script file "1_ScriptLinearPred" and its corresponding source code files. The EPA and the Kenneth French Data Library websites provide data dictionaries for both datasets.

## Code

### Abstract
The provided files offer guidance on how to apply transformed-linear prediction in both simulation studies and applications. By following the script files, users can generate optimized transformed-linear predictions, uncertainty quantification, and assess the coverage rate, etc.

### Description
The provided files contains comprehensive instructions and functions for applying transformed-linear prediction in simulation studies and applications, along with reproducing results in Sections 4, 5, and 6 of the manuscript:

### 1. Script and its source files for simulation studies.

* 1_ScriptLinearPred.R: A script file that outlines the steps to reproduce results in Section 4 and 5 of the manuscript. The script includes descriptions of the input arguments for each function used.
* 2_TransformedOperations.R: Functions that define transformed linear operations.
* 3_TailRatio.R: A function that illustrates a tail ratio in subset $V_+^q$ described in Section 4 of the manuscript.
* 4_PredictionError.R: A function that calculates the tail ratio of prediction error and reproduces Figure 3. 
* 5_TPDM_Ang.R: A function that simulates a vector X generated by a higher dimensional matrix A and returns a data matrix 'X_t'. 
* 6_TPDM_Est.R: A function that estimates the TPDM and obtains the 'prediction' inner product matrix (IPM).
* 7_Ang_CPD.R: A function that implements completely positive decomposition for a 2x2 prediction inner product matrix.
* 8_JointRegion.R: A function that creates the 95% 'joint polar region' by decomposing the prediction IPM and reproduces Figure 4.
    * Note: To create the 95% joint polar region, Xhats in the test set must be calculated first. See line 65-66 in the script file.
* 9_Ang_ked.R: A function that performs kernel density estimation for a bivariate angular density 'h'.
    * Requires the "VGAM" and "ks" packages to perform kernel density estimation.
* 10_CondDens.R: A function that provides an approximate conditional density given $\hat{X}$ is large.
* 11_CondInterval.R: A function that creates the 95% 'prediction interval' and returns 1 if its interval contains a predictand.
* 12_CondDenPlot.R: A function that plots an approximate conditional density with the 95% prediction interval given $\hat{X}$ is large.
* 13_AssessCoverage.R: A function that assesses a coverage rate.
* 14_CondIntervalPlot.R: A function that plots prediction intervals for the largest 5% given $\hat{X}$ is large.
* 16_Mix_ECDF_GPD.R: A function that returns a mixture distribution of an empirical CDF below the threshold and a generalized Pareto distribution above the threshold.

### 2. Script and its source files for applications.
* A1_DataCleaningNO2.R: A script that creates a new data frame with daily $\text{NO}_2$ data from five stations, after detrending.
    * Requires the "readr" and "dplyr" packages to load and process multiple CSV files, and the "ismev" package to estimate parameters.
* A2_Mov_Avg.R: A function that calculates the moving average and moving standard deviation of a time series.
* A3_ScriptNO2Mixture.R: A script that reproduces the results of the air pollution application in subsection 6.1 of the manuscript.
* A3_ScriptNO2MixtureMissing.R: A script that reproduces the results of the air pollution application when values at four stations are observed, but none are taken at Alexandria.
* A4_ScriptPortfolios.R: A script that reproduces the results of the industry portfolios application in subsection 6.2 of the manuscript.
    * Requires the "Hmisc" package to use the 'approxExtrap' function for linear extrapolation.

## Output files

### Rdata files for simulation study:

* AngCPD_Out_P.Rdata: contains the angles and angular masses obtained by decomposing the estimated prediction IPM using completely positive decomposition.
* PredIntTop5.Rdata: contains saved 95% prediction intervals for the largest 5% data determined by Xhat.

### Rdata files for applications:

####  Air pollution

* NewData.Rdata: contains preprocessed NO2 data for five stations.
* Train_NO2_Alx.Rdata: contains a training set with Pareto margins added.
* Test_NO2_Alx.Rdata: contains a test set with Pareto margins added.
* AngCPD_NO2_mix_Alx.Rdata: contains saved angles and angular masses by decomposition of a prediction IPM estimated in a train set.
* recentData.Rdata: contains a NO2 dataset for 2016-2020 when four stations are observed and no observation is taken at Alexandria.
* NewDat.Rdata: contains a dataset from 'NewData.Rdata' with Pareto margins added.
* AngCPD_NO2_mix_comb.Rdata: contains saved angles and angular masses when the NewData.Rdata is used as a whole.

####  Industry portfolios

* negPf_train.Rdata: contains a training set with Pareto margins added.
* negPf_test.Rdata: contains a test set with Pareto margins added.
* N_train.Rdata: contains a training set with Gaussian margins added.
* N_test.Rdata: contains a test set with Gaussian margins added.
* AngCPD_coal.Rdata, AngCPD_paper.Rdata, AngCPD_beer.Rdata: contains saved angles and angular masses by decomposition of a prediction IPM estimated in a training set for the corresponding company.

7_Ang_CPD.R uses completely positive decomposition to obtain a $2 \times q$ nonnegative matrix $B$. The computational time depends on the dimension of $q$ and the number of iterations in a for loop. By default, we choose $q$ to be 9 and set the number of iterations to 5000. This function takes approximately 2-3 minutes to obtain about 500 angles and angular masses. If users want to skip this function, they can directly use the corresponding Rdata files to reproduce the results.




