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
 There are two files under the folder: ``/Ranking_Inference_Reproduction/Real_Data/"  <br /> 1. Table 4, <br />2. Table 5-6 <br /> 
### Folder `/Ranking_Inference_Reproduction/Real_Data/':
 It contains the codes for reproduce results in Table 4 and Table 5-6. The document `jester_1.csv and jester_2.csv' in that folder is downloaded from the website mentioned above. <br />
 
 The code file **Table4.R** contains the codes for reproducing Table 4 of the real data analysis. The setting is given in line 96 (corrsponds to setting of Table 4), one only needs to run this code file directly. The outputs are given in line 133-134 and lines 265-271. We also mark their corresponding reproduced columns in Table 4.<br />
 
  **Run time** of **Table4.R** is less than 10 min.
 
  The code file **Table5-6.R** contains the codes for reproducing Tables 5 and 6 of the real data analysis. The setting of is given in line 101, one only need to set up these parameters according to the different settings of Table 5 and 6 in our paper and run this code file directly. The outputs are given in lines 269-279. We also mark their corresponding reproduced columns in Table 5-6.<br />
 
 **Run time** of **Table5-6.R** is less than 10 min.
 
 
# 2. **Simulation Studies** <br />

## Codes Description and Implementation Details:

### File `/Ranking_Inference_Reproduction/Simulation/Figure1-2.R':
This is the code file Figure1-2.R contains codes for reproducing the methodology of Figures 1 and 2. <br />
The setting of the code is given in line 34. Here we let L be fixed and p vary. (Following the caption of Figure 2, we can also let p fixed and let L vary). The output of this simulation is given in lines 66-67, where the l_2 and l_{\infty} errors are recorded for 500 times. The mean value and standard errors can be computed via these outputs directly.<br />

**Run time** of **pcr_adequate.R** for any given setting described above is around xxx hrs.

### `/Ranking_Inference_Reproduction/Simulation/Figure_3.R':
This is the code file Figure_3.R contains codes for reproducing the methodology of right panel of Figure 3. <br />
The setting of the code is given in lines 91-92. Here we follow the setting of figure 3 in the codes (n=60,p=0.05,L=80). The output of this simulation is given in line 127, where the empirical p-values are recorded. With these data, one is able to reproduce the right_panel of figure 3.<br />

**Run time** of **pcr_adequate.R** for any given setting described above is around xxx hrs.

### `/Ranking_Inference_Reproduction/Simulation/Figure_7.R':

This is the code file Figure_7.R contains codes for reproducing the methodology of right panel of Figure 7 in appendix (also left panel of Figure 3). <br />
The setting of the code is given in lines 49 and 52. Here we follow the setting of figure 7 by letting L=10 and p=(0.008,0.015,0.03) respectively. One is also able to set $L$ by different numbers (e.g. 5, 10, 20 suggested in figure 7). The output of this simulation is given in line 79, where we record the distribution of the MLE. With these data, one is able to reproduce the right_panel of figure 7 and left panel of Figure 3.<br />

**Run time** of **inference.R** for any given setting described above is around xxx hrs.

### Folder `/Ranking_Inference_Reproduction/Simulation/Table1-3.R':

This is the code file Table1-3.R contains codes for reproducing the methodology of Tables 1-3. <br />
The setting of the code is given in lines 144. The output of this simulation is given in lines 435-481, where we specify the meanings of all outputs and their corresponding position in Tables 1-3.<br />

**Run time** of **Table1-3.R** for any given setting described above is xxx hrs.

