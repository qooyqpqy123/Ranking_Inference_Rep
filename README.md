# Reproduction for the methodology of "Ranking Inferences Based on the Top Choice of Multiway Comparisons"

# 1. **Real Data Analysis** <br />

## Data  <br />

### 1.1. Data Description <br />
FRED-MD is a large macroeconomic database designed for the empirical analysis of “big data.” The datasets of monthly and quarterly observations mimic the coverage of datasets already used in the literature, but they add three appealing features. They are updated in real-time through the FRED database. They are publicly accessible, facilitating the replication of empirical work.  <br />
### 1.2. Data Availability <br />
 The original Dataset can be found in https://research.stlouisfed.org/econ/mccracken/fred-databases/. <br />
 In our .zip file, the original Dataset is provided in the folder ``Code_Reproduction/Real_Data/Pre_process/current.csv." <br />
 The description of the data is provided in "Code_Reproduction/Real_Data/Pre_process/Description_of_Variables.pdf."  <br /> This description describes the meaning of all columns of this Dataset.

## Codes Description and Implementation Details:
 There are four folders under the folder: ``Code_Reproduction/Real_Data."  <br /> 1. Pre_process, <br />2. Prediction_Table_3 <br />  3. PCR_Adequate_Table4, <br /> 4. Sparse_Adequate_Table4. <br />
### Folder `Code_Reproduction/Real_Data/Pre_process':
 It contains the codes for real data pre-processing. The document `current.csv' in that folder is downloaded from the website mentioned above. <br />
 
 The code **realdata_read.R** contains the codes for pre-processing the data. The line 4 in "realdata_read.R", we have a variable called: "choose_time". If one lets choose_time=1 and run the codes, the output is the processed Dataset **realdata_115_126.csv**. If one lets choose_time=2, and run the codes, the output is the processed Dataset **realdata_189_126.csv**.<br />
 
 **Run time** of **realdata_read.R** is less than 3 min.
 
 Finally, for document "Description_of_Variables.pdf" inside that document, it contains a detailed description of all variables in the original dataset "current.csv". <br />
 

### Folder `Code_Reproduction/Real_Data/Prediction_Table3':
This folder contains the codes for reproducing the prediction results of Table 3 in our paper. <br />

*For code file **realdata_prediction_115.R**, it first reads the data "realdata_115_126.csv" in that folder. On line 54 of this file, there is a variable "j ."If we set j=49 and run the codes, we reproduce the line "HOUSTNE" line during "08.2010-02.2020" in Table 3. If we set j=81 and run the codes, we reproduce the line "GS5" line during "08.2010-02.2020" in Table 3. <br />

*For code file **realdata_prediction_189.R**, it first reads the data "realdata_189_126.csv" in that folder. On the line 54 of this file, there is a variable "j". If we set j=49, and run the codes, we reproduce the line "HOUSTNE" line during "02.1992-10.2007" in Table 3. If we set j=81 and run the codes, we reproduce the line "GS5" line during "02.1992-10.2007" in Table 3. <br />

For both .R files, the output is printed in lines 149-161 and 189 for different methods: Lasso (sparse linear regression), PCR (latent factor regression), Ridge (Ridge regression), El-Net (Elastic Net), RF (Random Forest), and FarmSelect (Factor adjusted Lasso), respectively.

 **Run time** of **realdata_prediction_115.R** and **realdata_prediction_189.R** are less than 10 min.

### Folder `Code_Reproduction/Real_Data/PCR_Adequate_Table4':
This folder contains the codes for reproducing the hypothesis testing results of column "LA_factor" of Table 4 in our paper. <br />

*For code file **Real_data_PCR_115.R**, it first reads the data "realdata_115_126.csv" in that folder. On line 17 of this file, there is a variable "j." If we set j=49 and run the codes, we reproduce the line "HOUSTNE", column "LA_factor," during "08.2010-02.2020" in Table 4. If we set j=81 and run the codes, we reproduce the line "GS5" line, column "LA_factor", during "08.2010-02.2020" in Table 4. <br />


*For code file **Real_data_PCR_189.R**, it first reads the data "realdata_189_126.csv" in that folder. Line 17 of this file has a variable "j". If we set j=49, and run the codes, we reproduce the line "HOUSTNE", column "LA_factor," during "02.1992-10.2007" in Table 4. If we set j=81 and run the codes, we reproduce the line "GS5" line, column "LA_factor", during "02.1992-10.2007" in Table 4. <br />

For both .R files, the outputs are p-values of the test on whether PCR is adequate.

**Run time** of **Real_data_PCR_115.R** and **Real_data_PCR_189.R** are less than 10 min.

### Folder `Code_Reproduction/Real_Data/Sparse_Adequate_Table4':
This folder contains the codes for reproducing the hypothesis testing results of column "SP_Linear" of Table 4 in our paper.

*For code file **Real_data_Sparse_115.R**, it first reads the data "realdata_115_126.csv" in that folder. On the line 94 of this file, there is a variable "j". If we set j=49, and run the codes, we reproduce the line "HOUSTNE", column "SP_Linear," during "08.2010-02.2020" in Table 4. If we set j=81 and run the codes, we reproduce the line "GS5" line, column "SP_Linear", during "08.2010-02.2020" in Table 4. <br />


*For code file **Real_data_Sparse_189.R**, it first reads the data "realdata_189_126.csv" in that folder. On line 94 of this file, there is a variable "j". If we set j=49, and run the codes, we reproduce the line "HOUSTNE", column "SP_Linear," during "02.1992-10.2007" in Table 4. If we set j=81 and run the codes, we reproduce the line "GS5" line, column "SP_Linear", during "02.1992-10.2007" in Table 4. <br />
 
For both .R files, the outputs are p-values of the test on whether Sparse Regression is adequate.

**Run time** of **Real_data_Sparse_115.R** and **Real_data_Sparse_189.R** are less than 10 min.

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
