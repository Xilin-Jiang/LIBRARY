# When should I organise code?
Organising code is a good habit -- however, due to the vast amount of scripts generated on a daily basis, it is impossible to keep track of all methods I used. In practice, whenever I need to **re-implement ** a previous pipeline, please organise all the code from that pipeline here. 

# What I should do?
First, keep a good habit when writing code! 

### ----- on the cluster: 
###### Step 1: Create the new folder under each project for the .sh files.
###### Step 2: Always contain a pipeline.txt file keep track of all the .sh files in each folder. 
###### Step 3: Use the same name of the .R file (which usually under the project folder) and .sh file (which is usually in the sub-folder).

### ----- on the local computer:
Organise the functions into a separate file, and source it; this will save so much efforts when packaging! Keep track of the code for the paper when actually started to wrap up the analysis and writing up. 

### ----- when organising code here:
When orgnising code here, try to be complete! Even you only use a small part of the pipeline, make sure you have small input example, .RData (for figures) and also keep track of other part of the pipeline as it is usually very little additional work. Follow the steps at the end of the README.md

# LIBRARY

Organising all the reusable code for future references. The code should be organised in the each folder with the name of the paper, which makes it easier to trace. 

# Structure
The scripts are structured at two level (1) top-level contains a couple of scripts that contains the code that are generated for projects that could not be  categories in the paper folder. For example, I wrote the code for a mixture of Dirichlet distribution which was not used in the ATM paper. (2) Paper folders, which contains the pipelines (organised) for this paper. It is strongly recommended to summarise all relevant code in the same pipeline. For example, using bolt_lmm to compute lmm and PRS should be summarised in the same scripts.   

### Top-level scripts
Top-level scripts should be carefully used. Only when it is really hard to organise the code in the paper folders, the relevant code could be put into the top level scripts. At the moment there are only three files: genetic_tools.R, inference.R, and plotting_functions.R.

### Paper folders
Each Relevant functions are organised into scripts. There are usually two scripts of the same name, one for .R code and one for .sh code; the namethe script: should describe the task using the published methods. (e.g. LDSC_seg, BOLT_lmm, SuSiE) There will also be a folder of the same name of the file: saving example Rdata and example figures for reproducibility.

Each section of the scripts represent a small task, and should be started with a one-line header describing the task. Functions should always be accompanied by a simple example. To increase readability example, save all examples objects in the folder. The naming ruld of the RData should be: SCRIPTNAME_task.RData  

Each paper folder should contain a script named Figures.R, which saves the implementation of the important figures. The figures should be organised as the order they presented in the paper and the discription should contain: Figure X (barplot/heatmap/line/point/boxplot). 

#### Summary

###### Step 1: Create the paper folder. (Name: Jiang_YEAR_JOURNAL)
###### Step 2: Create the Figures.R file.
###### Step 3: For each task XXX: create XXX.R, XXX.sh and XXX/
###### Step 4: Organise all the code; these step is usually done at the revision stage, when I need to rerun many analysis and redo the figures.  

# Note
No efforts is made to formally package the code, though relevant comments are kept to provide context. The goal of this LIBRARY is organising existing code (mostly from data analysis scripts) for future references. Important code should be packaged and should not be simply listed here.  
