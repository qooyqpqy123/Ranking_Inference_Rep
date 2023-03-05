# Reproduction for the methodology of "Ranking Inferences Based on the Top Choice of Multiway Comparisons"

# 1. **Real Data Analysis** <br />

## Data  <br />

### 1.1. Data Description <br />
4.1 million continuous ratings of 100 jokes from 73421 users: collected between Apirl 1999-May 2003. <br />
### 1.2. Data Availability <br />
 The original Dataset and detailed description can be found in https://goldberg.berkeley.edu/jester-data/. <br />
 In our .zip file, the original Dataset is provided in the folder ``/Ranking_Inference_Reproduction/Real_Data/jester_1.csv" and `/Ranking_Inference_Reproduction/Real_Data/jester_2.csv" <br />
Description:
Data files contain anonymous ratings data from 73,421 users. <br />
Ratings are real values ranging from -10.00 to +10.00 (the value "99" corresponds to "null" = "not rated").<br />
One row per user <br />
The first column gives the number of jokes rated by that user. The next 100 columns give the ratings for jokes 01 - 100. <br />

## Codes Description and Implementation Details:
 There are two folders under the folder: ``/Ranking_Inference_Reproduction/Real_Data/"  <br /> 1. Table 4, <br />2. Table 5-6 <br /> 
### Folder `/Ranking_Inference_Reproduction/Real_Data/':
 It contains the codes for reproduce results in Table 4 and Table 5-6. The document `jester_1.csv and jester_2.csv' in that folder is downloaded from the website mentioned above. <br />
 
 The code **Table4.R** contains the codes for reproducing Table 4 of the real data analysis. The setting of the paper is given in line 96, one only need to run this code file directly. The outputs are given in line 133-134 and lines 265-271. We also mark their corresponding reproduced columns in Table 4.<br />
 
  **Run time** of **Table4.R** is less than 10 min.
 
  The code **Table5-6.R** contains the codes for reproducing Tables 5 and 6 of the real data analysis. The setting of the paper is given in line 101, one only need to set up these parameters according to the different settings of Table 5 and 6 in our paper and run this code file directly. The outputs are given in lines 269-279. We also mark their corresponding reproduced columns in Table 5-6.<br />
 
 **Run time** of **Table5-6.R** is less than 10 min.
 
 
# 2. **Simulation Studies** <br />

## Codes Description and Implementation Details:

### Folder `Code_Reproduction/simulation/Consitency_Figure1/':
There are three sub-folders inside this folder, namely, Gaussian, Uniform, and Heavy-tail. <br />

For folders Gaussian and Uniform, they contain codes **Gauss_consistency.R** and **Uniform_consistency.R**, respectively. By running these codes directly, we reproduce the methodology for generating results for Figure 1 (a) and (b). <br />

For folder Heavy-tail: it contains codes: **huber_estimation.R** and **T-consistency.R**. By running these codes directly, we reproduce the methodology for generating the red and blue lines of (c) of Figure 1.

**Run time** of **Gauss_consistency.R**, **Uniform_consistency.R**, **huber_estimation.R** and **T-consistency.R** are around 5 hrs, respectively.

### Folder `Code_Reproduction/simulation/Simulation_Table1/':
There is one file **pcr_adequate.R** under this folder. It reproduces the methodology of Table 1. <br /> 
There are several parameters in this file that need to be pre-specified. They correspond to different simulation settings of the paper. <br /> 

1. First, on line 7, if we let the parameter "choose_dim=1", it outputs the results with dimension p=200. If we let the parameter "choose_dim=2", it outputs the results with dimension p=500.  <br /> 

2. Second, on line 10, if we let the parameter "choose=1", it outputs the results under the setting where factors F and idiosyncratic components are generated via setting 1 in section 5.2 of this paper. If we let the parameter "choose=2",it outputs the results under the setting where factors F and idiosyncratic components are generated via setting 2 in section 5.2 of this paper. <br /> 

3. Third, on line 16, if we let the parameter "choose_noise=1", it outputs the results when the noise is generated based on Gaussian distribution. If we let the parameter "choose_noise=2", it outputs the results under the setting where the noise follows a uniform distribution.  <br /> 

**Run time** of **pcr_adequate.R** for any given setting described above is around 12 hrs.

### Folder `Code_Reproduction/simulation/Simulation_Table2/':

There is one file **inference.R** under this folder. It reproduces the methodology of Table 2. <br /> 
There are several parameters in this file that need to be pre-specified. They correspond to different simulation settings of the paper. <br /> 

1. First, on line 93, if we let the parameter "choose_p=1", it outputs the results with dimension p=250. If we let the parameter "choose_p=2", it outputs the results with dimension p=600.  <br /> 

2. Second, on line 94, if we let the parameter "choose_noise=1", it outputs the results where the noise is generated based on Gaussian distribution. If we let the parameter "choose_noise=2", it outputs the results under the setting where the noise follows a uniform distribution. <br /> 

3. Third, on line 95, if we let the parameter "choose_mix=1", it outputs the results under the setting where the covariates are generated, i.i.d. If we let the parameter "choose_mix=2",it outputs the results under the setting where the covariates are strong mixing. <br /> 

**Run time** of **inference.R** for any given setting described above is around 5 hrs.

### Folder `Code_Reproduction/simulation/Prediction_Section_B.1/':

There are four files inside this folder. **Prediction_Table1.R**, **Prediction_Table2.R**, **Prediction_Table3.R**, and **Prediction_Table4.R**   <br /> 
They reproduce Tables 1-4 in Appendix, respectively. 

**Run time** of **Prediction_Table1.R**, **Prediction_Table2.R**, **Prediction_Table3.R**, and **Prediction_Table4.R** for any given setting described above is less than 30 min.

### Folder `Code_Reproduction/simulation/Figure7_appendix/':

There is one file **qqplot_figure7.R** in the folder. It reproduces the methodology to generate Figure 7 in Appendix. The implementation is to run codes in 1-160 lines to generate "sure.csv". After that, one needs to read in this "sure.csv" to plot the Q-Q-plot (corresponds to lines 164-171). 


**Run time** of **qqplot_figure7.R** for any given setting described above is around 5 hrs.
