# Project title

Program codes to generate results in "Mixture Prior for Sparse Signals with Dependent Covariance Structure".

## Overview

The folder “sim_study” has the codes for simulation study, and the folder “real_data_study” has the codes for real data study. The folder “matlab_functions” includes libraries and functions.

## How to run codes for simulation study
To run the simulation codes, go to “sim_study” folder and follow the steps below: 
1. Find R codes. “Simulation_Uninormal.R” and “Simulation_Binormal.R” are main R programs for simulating uninormal and binormal distribution type, respectively.  “GenSigmaMatrix.R”, “GenBeta.R” and “PFAdecomp.R” are functions that are called in the main R programs. 
   
   Before running the main R program for simulating uninormal distribution: “Simulation_Uninormal.R”, you need to do the following:
   1. Create output directories for each value of sparsity under each distribution type and each type of covariance matrix. For example

        ```c
        Your Hard Drive Name/sim_study/sim_dat_file/Uninormal/equalcorr4/0.6
        Your Hard Drive Name/sim_study/sim_dat_file/Uninormal/equalcorr4/0.7
        Your Hard Drive Name/sim_study/sim_dat_file/Uninormal/equalcorr4/0.8
        Your Hard Drive Name/sim_study/sim_dat_file/Uninormal/equalcorr4/0.9
        ………
        ………
        Your Hard Drive Name/sim_study/sim_dat_file/Uninormal/Threefactor/0.6
        Your Hard Drive Name/sim_study/sim_dat_file/Uninormal/Threefactor/0.7
        Your Hard Drive Name/sim_study/sim_dat_file/Uninormal/Threefactor/0.8
        Your Hard Drive Name/sim_study/sim_dat_file/Uninormal/Threefactor/0.9
        ```

   2. Modify the path_name in “Simulation_Uninormal.R” so that the output directory is correctly changed to “Your Hard Drive Name/sim_dat_file”
   3. Run “Simulation_Uninormal.R”

   Before running the main R program for simulating binormal distribution: “Simulation_Binormal.R”, do the same processes in items i, ii, iii described above.

2. Next, find Matlab programs “Gen_Sim_Uninormal.m” and “Gen_Sim_Binormal.m” for uninormal and binormal distribution type, respectively. “EstimationProc.m” is a Matlab function.

   Before running Matlab programs “Gen_Sim_Uninormal.m” and “Gen_Sim_Binormal.m”, you need to do the following:
   1. Go to the directory:
        ```c
        Your Hard Drive Name/matlab_functions/figtree-0.9.1/figtree-0.9.1
        ```
      You can find a “readme.txt” here. Follow the instructions in that “readme.txt” to compile the files in the folder:
        ```c
        Your Hard Drive Name/matlab_functions/figtree-0.9.1/figtree-0.9.1/matlab/
        ```
   2. Modify the addpath in “Gen_Sim_Uninormal.m” and “Gen_Sim_Binormal.m”, so that the addpath directory is correctly changed to “Your Hard Drive Name/matlab_functions/NonParametric” and “Your Hard Drive Name/matlab_functions/ figtree-0.9.1/figtree-0.9.1/matlab”.
   3. Modify the csvread(sprintf(‘….’)) in “Gen_Sim_Uninormal.m” and “Gen_Sim_Binormal.m”, so that the csvread(sprintf(‘….’))  directory is correctly changed to “Your Hard Drive Name/sim_dat_file/” 
   4.	Modify the save(sprintf(‘….’)) in “Gen_Sim_Uninormal.m” and “Gen_Sim_Binormal.m”, so that the save(sprintf(‘….’))  directory is correctly changed to “Your Hard Drive Name/sim_dat_file/”
   5. Run Matlab programs “Gen_Sim_Uninormal.m” and “Gen_Sim_Binormal.m” for uninormal and binormal distribution type, respectively.


## How to run codes for real data study
To run the real data analysis codes, go to “real_data_study” folder and follow the steps below:
1. Find R codes. “1_gen_cor_gene.R”, “2_imp_cor.R”, and “3_gen_coef.R” are the 3 main R programs. “PFAdecomp.R” is a function that are called by the main R program. 
   
   Before running the R programs, you need to do the following:
   1. Create a data directory to store input and output data files, For example:
        ```c
        Your Hard Drive Name/datafile/
        ```
      In this directory, store the gene expression data files.
   2. Open “1_gen_cor_gene.R”, modify the path name in as.matrix(read.delim(“”)) so that the data import directory is correctly changed to “Your Hard Drive Name/datafile”

      Modify the path name in write.table(“”) so that the data export directory is correctly changed to “Your Hard Drive Name/datafile”

      Run “1_gen_cor_gene.R”.
   3. Open “2_imp_cor.R”, modify the path name in as.matrix(read.delim(“”)) so that the data import directory is correctly changed to “Your Hard Drive Name/datafile” 

      Modify the path name in write.table(“”), write.csv(“”) so that the data export directory is correctly changed to “Your Hard Drive Name/datafile”

      Modify the path name in eps_graph_name, pdf_graph_name, png_graph_name so that the graph export directory is correctly changed to “Your Hard Drive Name/datafile”

      Run “2_imp_cor.R”.

   4. Open “3_gen_coef.R”, modify the path name in as.matrix(read.delim(“”)) so that the data import directory is correctly changed to “Your Hard Drive Name/datafile” 

      Modify the path_name so that the data export directory is correctly changed to “Your Hard Drive Name/datafile”

      Run “3_gen_coef.R”.

2. Next, find Matlab programs. “chrna_data.m” is the main Matlab program. “EstimationProc.m” is a function.
   1. If you have compiled the Matlab functions in the following directory “Your Hard Drive Name/matlab_functions/figtree-0.9.1/figtree-0.9.1/matlab”, then go to the next step; otherwise, you can find a “readme.txt” in the directory
        ```c
        Your Hard Drive Name/matlab_functions/figtree-0.9.1/figtree-0.9.1
        ```
      Follow the instructions in that “readme.txt” to compile the files in the folder:
        ```c
        Your Hard Drive Name/matlab_functions/figtree-0.9.1/figtree-0.9.1/matlab/
        ```
   2. Modify the addpath in “chrna_data.m” so that the addpath directory is correctly changed to “Your Hard Drive Name/matlab_functions/NonParametric” and “Your Hard Drive Name/matlab_functions/figtree-0.9.1/figtree-0.9.1/matlab”.
   3. Modify the readtable(‘….’), writetable(‘….’), csvread(‘….’)  in “chrna_data.m” so that the directory is correctly changed to “Your Hard Drive Name/datafile/”
   4.	Run Matlab program “chrna_data.m”
