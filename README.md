# Coevolution

These folders contain all scripts used in this project

Folder Experiment 1, Experiment 2, Experiment 3 are split up by specific questions asked 
in paper

Each of theses three folders have the experiments split up by different aspects applied 
to each experiment (i.e. genetice architecture combinations, geographical landscapes, 
simulation specific changes).

In each individual simulation experiment (i.e both, newt_change, no_hair ...) there will 
be a csv, msprime_out, slimtree_out, and snakemake folder. The folders csv, msprime_out, 
and slimtree_out collect data; log files from SLiM, msprime trees, and tree files 
produced by SLiM. The snakemake folders hold the snakefile which can be used to reproduce 
all of the simulations that were run in these experiments. The .yamal files hold the 
parameters used in each experiment.

In the Experiment1 folder there are three folders that represent different geographical 
landscapes (for cost gradient, flat map, interaction gradient). The cost_grade folder has 
folders depicting weither the map was applied to both species, one species 
(newt or snake), and for both species but without newt/snake interaction. The flat_map 
folder has the three experiments ran without a geographical landscape (normal, no 
heritability, and no interaction). The interaction_grade folder holds just one experiment. 

Experiment2 and Experiment3 each have four folders for different genetic 
architecture combinations (these were all ran without geographical landscape)

The folder MSprime_script contains the python script used to make msprime trees

The R_scripts folder contains the script used to gather the data into more condensed 
files and create figures for the paper

The SLiM_script folder contains all of the simulations used in this paper. There is an 
additional test_data folder containing a geographical landscape (gradent.png) and a 
msprime generated file so that each SLiM simulation can be run on a test basis,

Versions:

Python Script:
python 3.10.2 
pyslim 0.700
tskit 0.4.1
msprime 1.1.0
numpy 1.22.3

R Script:
R version 4.3.2
Rstudio 2023.09.1+494
cowplot_1.1.1
scales_1.2.1
RColorBrewer_1.1-3
patchwork_1.1.3   
ggpubr_0.6.0       
kableExtra_1.3.4   
knitr_1.45         
gridGraphics_0.5-1
plyr_1.8.9         
reshape2_1.4.4     
lubridate_1.9.3    
forcats_1.0.0     
stringr_1.5.1      
dplyr_1.1.4        
purrr_1.0.2        
readr_2.1.4       
tidyr_1.3.0        
tibble_3.2.1       
tidyverse_2.0.0    
gridtext_0.1.5    
gridExtra_2.3      
colorspace_2.1-0   
viridis_0.6.4      
viridisLite_0.4.2 
ggplot2_3.4.4 

SLiM version 3.7
 
