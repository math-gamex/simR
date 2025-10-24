#############################################################################
#                                                                           #  
# THE DYNAMICS OF SON PREFERENCE, TECHNOLOGY DIFFUSION, AND FERTILITY       #
# DECLINE UNDERLYING DISTORTED SEX RATIOS AT BIRTH: A SIMULATION APPROACH   #
#                                                                           #
#                                                                           #
# Ridhi KASHYAP - kashyap@demogr.mpg.de/ridhi.kashyap@nuffield.ox.ac.uk     #
# Francisco VILLAVICENCIO - villavicencio@imada.sdu.dk                      #                  
#                                                                April 2016 #
#############################################################################


INTRODUCTION
------------

This file offers a brief overview of the R code that executes the model. We have also commented the R-scripts as much as possible in order to help an external user better understand the code. Please contact us if something is not clear.

Files and folders available: 
  1) SexRatioBrith.R: Contains the code the user runs to execute the model.
  2) FunctionsSRB.R: Contains some of the functions used in the model. These functions are automatically uploaded when executing SexRatioBrith.R.
  3) Data: Folder with empirical data to be introduced as inputs into the model.
  4) Sensitivity: This folder contains several tables with different combinations of the model parameters to be used in the sensitivity analysis. The user can use those tables, or generate his/her own table with different parameter combinations to study model behavior.
  5) SimulationResults: This folder contains Rdata files with results needed to plot Figures 4, 6 and 7 in the Demography paper. 

The code is optimized to take advantage of parallel computing using the 'snowfall' R package. The 'parallel' package is also required in order to automatically detect the number of cores of local computer using the function detectCores(). If the number of cores is known, the user can replace detectCores() by that value in SexRatioBrith.R, and then the 'parallel' package is no longer necessary.


HOW TO RUN THE MODEL
--------------------

The file SexRatioBrith.R is divided into different sections:

1) PACKAGES, LIBRARIES AND WORKING DIRECTORY

In this section the user needs to specify the working directory (path) of the local machine in which the simulation code and input data are placed. The working directory should contain the two source code files, and the subfolders 'Data', 'Sensitivity', and 'SimulationResults'.

2) INITIAL PARAMETERS. In this section, the user needs to specify:

  - sensitivity = TRUE: This will carry out a sensitivity analysis. We suggest to use sensitivity = FALSE by default.
  - Which model is going to be tested: Empirical or Theoretical. We suggest to use Empirical by default. The empirical specification will take empirical fertility and mortality rates, and son preference values (Fig. 2 in the paper) as inputs. 
  - Several global parameters such as initial population size, SRB, or the age range. If initial year and / or number of iterations are changed, be aware to have input data corresponding to that period on the files of the 'Data' folder.
  - Model parameters: Stored in vectors;  all vector parameters should have the same length.
     * If n = 1, the same parameters are used in all parallel simulations.
     * If n > 1, n parallel simulations are run (assuming that the computer has at least n cores), using the first value of the vectors in the first core, the second value of the vectors in the second core, etc. Alternatively, the user may define a parameter vector like for example 'alphaPar <- c(0.075, .080)', which imply that the first core will use alpha = 0.075, and in the second core alpha = .080. 

3) DATA AND INITIAILIZATION OF THE MODEL

In this section, the user needs to specify which country is going to be tested (South Korea or India), and then all the necessary input data is uploaded.

4) RUN THE MODEL WITH PARALLEL COMPUTING

No modifications of the code need to be done by the user in this section. The number of parallel simulations (ncpus) is assinged depending on the model parameters, the scenario (Empirical or Theoretical) and if a sensitivity analysis is desired. Next, the code uses functions from the 'snowfall' R package to run the parallel simulations. 

5)  RESULTS - SAVE OUTPUT

The output of the model is stored in the list 'out', which collects the results for all the parallel simulations (e.g. out[[1]], out[[2]], etc.). Different indicators are provided, such as the proportion of boys or the TFR in each year of the simulation. Additionally, it is possible to save the output in an R.Data file. The user may change the name of that file which, by default, will be stored in the working directory previously specified.

6) PLOTTING

Code needed to reproduce Figures 4, 6 and 7 from the paper.


FUNCTIONS
---------

The functionsSRB.R file contains the following functions:
  - get.data(): This function uploads empirical data necessary to run the model (fertility rates, death rates, population structure). Data from two countries are available: South Korea and India (for calibration to Indian SRBs, see Kashyap and Villavicencio (2016)*. Most of the data are from the United Nations World Population Prospects 2012.
  - funcSonPref(): This function estimates the son preference by cohort and period of South Korea or India for the empirical scenarios by fitting logistic regressions to observed data points from survey data. It also provides a son preference matrix for the theoretical scenarios, such as those reported in Figure 6 of the paper.
  - initialization(): This function is used to obtain the initial population that will be used in the model. To approximate the initial parity distribution of women from 1980 onwards, we initialized the model 35 years earlier, in 1945, to allow all women in reproductive ages (15+) to complete their fertility careers by 1980 and have their children belong to the starting population of 1980.
  - run.simulations(): This one is the main function of the model, used to execute the simulations. When using parallel computing, this is the function that is sent to each core in order to run several parallel simulations.
  - findSonPref: This small function is just used to assign the son preference to each individual during the simulations.
  - There are a series of extract data functions at the end of the file. These functions extract data on TFR, SRBs, RMSEs, and other indicators from simulation RData files to study model behavior and for plotting results. 

---
*Kashyap, R. and F. Villavicencio (2016). An agent-based model of sex ratio at birth distortions. In A. Grow and J. Van Bavel (Eds.), Agent-Based Modelling in Population Studies: Concepts, Methods, and Applications, Chap. 12, Cham: Springer.
