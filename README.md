Scripts used for Genetic architecture, spatial heterogeneity, and the arms race between 
newts and snakes: Exploring coevolution with simulations

These folders contain all scripts used in this project

********************************************************************************

This Study:

To better understand the dynamic interaction between spatial structure, 
genetic architecture, and coevolution, we conducted a simulation study, 
exploring a range of situations plausible for the Taricha newt Thamnophis 
garter snake system. To better understand the process of adaptive evolution 
across geographical space, this paper asks three main questions: 

(1) How does the genetic architecture of the traits in newts and snakes 
affect how they coevolve?
(2) Under what situations do we get spatial patterns of correlated traits 
as we see in the real world?
(3) How fast does coevolution increase resistance and toxicity in these 
organisms with different combinations of genetic architectures?

In particular, we compare different levels of mutational variance and polygenicity using 
individual-based simulations of continuous geographic space. The results complement field 
observations by describing situations that are consistent with empirical observations, 
and exploring other possible outcomes.

********************************************************************************

Caudill, V., & Ralph, P. L. (2023). Genetic architecture, spatial heterogeneity, and the 
coevolutionary arms race between newts and snakes. bioRxiv. 
https://www.biorxiv.org/content/10.1101/2023.12.07.570693v2

Authors and contact details
Victoria Caudill: vcaudill@uoregon.edu
Peter Ralph: plr@uoregon.edu

Victoria and Peter are responsible for writing the simulations. Victoria was responsible 
for writing the code to analize the data.

********************************************************************************

Folder Experiment 1, Experiment 2, Experiment 3 are split up by specific questions asked 
in paper. The letters next to mu's and v's correspond to table 1 in the paper. 

Each of theses three folders have the experiments split up by different aspects applied 
to each experiment (i.e. genetic architecture combinations, geographical landscapes, 
simulation specific changes).

In each individual simulation experiment there will 
be a csv, msprime_out, slimtree_out, and snakemake folder. The folders csv, msprime_out, 
and slimtree_out collect data; log files from SLiM, msprime trees, and tree files 
produced by SLiM. The snakemake folders hold the snakefile which can be used to reproduce 
all of the simulations that were run in these experiments. The .yamal files hold the 
parameters used in each experiment.

Each of the csv folders have three text files (all_cor, all_grid, and all_lit)
all_cor - contains spatial correlation information taken throughout the simulation
all_cor variables: 
generation - simulation time point,
mean_newt_pheno_By_mean_snake_pheno - correlation between average local newt toxicity and
average local snake resistance,
num_newts_By_num_snakes - correlation between local number of newts and local number of 
snakes, 
sum_newt_pheno_By_num_snake - correlation between average local newt toxicity and local 
number of snakes,
sum_snake_pheno_By_num_newt - correlation between average local snake resistance and 
local number of newts,
sum_newt_pheno_By_num_newt - correlation between average local newt toxicity and local 
number of newts,
sum_snake_pheno_By_snake_newt - correlation between average local snake resistance and 
local number of snakes,
rep - simulation trial,
snake_mu_rate - snake mutation rate,
snake_mu_effect_sd - snake mutation effect standard deviation,
newt_mu_rate - newt mutation rate,
newt_mu_effect_sd - newt mutation effect standard deviation

all_grid - contains information from smaller localized areas (grids) within the simulation 
(100 smaller squares within each simulation)
all_grid variables: 
generation  - simulation time point,
newt_num_ind_0 [newt_num_ind_1, newt_num_ind_2, ..., newt_num_ind_99] - number of newts,
snake_num_ind_0 [snake_num_ind_1, snake_num_ind_2, ..., snake_num_ind_99] - number of 
snakes,
newt_max_pheno_0 [newt_max_pheno_1, newt_max_pheno_2, ... newt_max_pheno_99] - max level 
of newt toxicity within the grid,
snake_max_pheno_0 [snake_max_pheno_1, snake_max_pheno_2, ... snake_max_pheno_99] - max
level of snake resistance within the grid,
newt_mean_pheno_0 [newt_mean_pheno_1, newt_mean_pheno_2, ... newt_mean_pheno_99] - average
level of newt toxicity within the grid,
snake_mean_pheno_0 [snake_mean_pheno_1, snake_mean_pheno_2, ... snake_mean_pheno_99] - 
average level of snake resistance within the grid,
newt_sd_pheno_0 [newt_sd_pheno_1, newt_sd_pheno_2, ... newt_sd_pheno_99] - standard 
deviation of newt toxicity within the grid,
snake_sd_pheno_0 [snake_sd_pheno_1, snake_sd_pheno_2, ... snake_sd_pheno_99] - standard 
deviation of snake resistance within the grid,
newt_min_pheno_0 [newt_min_pheno_1, newt_min_pheno_2, ... newt_min_pheno_99] - min level 
of newt toxicity within the grid,
snake_min_pheno_0 [snake_min_pheno_1, snake_min_pheno_2, ... snake_min_pheno_99] - min
level of snake resistance within the grid,
rep - simulation trial,
snake_mu_rate - snake mutation rate,
snake_mu_effect_sd - snake mutation effect standard deviation,
newt_mu_rate - newt mutation rate,
newt_mu_effect_sd - newt mutation effect standard deviation

all_lit -  contains averaged information from the entire simulation
all_lit variables:
generation  - simulation time point,
Newt_age - average age of newts (in generations),
Snake_age - average age of snakes (in generations),
Newt_density - newt density,
Snake_density - snake density,
Newt_min_Pheno - smallest level of toxicity,
Newt_max_Pheno - largest level of toxicity,,
Newt_mean_Pheno - average level of toxicity,,
Newt_sd_Pheno - standard deviation of toxicity,,
Snake_min_Pheno - smallest level of resistance,
Snake_max_Pheno - largest level of resistance,,
Snake_mean_Pheno - average level of resistance,,
Snake_sd_Pheno - standard deviation level of resistance,,
Newt_pop_size - number of newt individuals,
Snake_pop_size - number of snake individuals,
mean_newts_eaten - how many newts were killed by snakes,
sd_newts_eaten - standard deviation of how many newts were killed by snakes (number of 
newts killed/newts age),
mean_snakes_eaten - how many snakes were killed by newts,
sd_snakes_eaten - standard deviation of how many snakes were killed by newts (number of 
snakes killed/snakes age)
beta_n - correlation between level of newt toxicity and an individuals fitness,
beta_s - correlation between level of snake resistance and an individuals fitness,
snake_found_newt - number of newts snakes could find,
newt_found - number of interactions,
snake_deaths - number of snakes that failed to eat a newt,
newt_deaths - number of newts that failed to survive a snake encounter,
time - simulation output for how long it has been running (similar to seconds),
rep - simulation trial,
snake_mu_rate - snake mutation rate,
snake_mu_effect_sd - snake mutation effect standard deviation,
newt_mu_rate - newt mutation rate,
newt_mu_effect_sd - newt mutation effect standard deviation

In the Experiment1 folder there are three folders that represent different geographical 
landscapes (for cost gradient (cost_grade), flat map (flat_map), interaction gradient 
(interaction_grade)). The cost_grade folder has folders depicting whether the map was 
applied to both species (both), one species (newt_change or snake_change), and for 
both species but without newt/snake interaction (both_no_interaction). The flat_map 
folder has the three experiments ran without a geographical landscape (normal, no 
heritability (no_hair), and no interaction(no_interaction)). 
The interaction_grade folder holds just one experiment (simulation with an gradient 
for interaction rate). 


Experiment2 and Experiment3 each have four folders for different genetic 
architecture combinations (these were all ran without geographical landscape)

The folder MSprime_script contains the python script used to make msprime trees

The R_scripts folder contains the script used to gather the data into more condensed 
files and create figures for the paper

The SLiM_script folder contains all of the simulations used in this paper. There is an 
additional test_data folder containing a geographical landscape (gradent.png) and a 
msprime generated file so that each SLiM simulation can be run on a test basis,

Folder hierarchy/ description
Experiment1 (varies both mutation rates and mutation effect sizes.):
	cost_grade (simulations with a cost gradent)
		both (cost gradent for both newts and snakes)
			csv (all slim outputs)
				all_cor.txt (correlations)
				all_grid.txt (grid based)
				all_lit.txt (whole simulation values)
			msprime_out (trees from msprime)
			slimtree_out (trees from slim)
			snamemake (simulation set up)
		newt_change (cost gradient only for the newt)
			csv (all slim outputs)
				all_cor.txt (correlations)
				all_grid.txt (grid based)
				all_lit.txt (whole simulation values)
			msprime_out (trees from msprime)
			slimtree_out (trees from slim)
			snamemake (simulation set up)
		snake_change (cost gradient only for the snake)
			csv (all slim outputs)
				all_cor.txt (correlations)
				all_grid.txt (grid based)
				all_lit.txt (whole simulation values)
			msprime_out (trees from msprime)
			slimtree_out (trees from slim)
			snamemake (simulation set up)
		both_no_interaction (cost gradient for both newt and snake with no interaction)
			csv (all slim outputs)
				all_cor.txt (correlations)
				all_grid.txt (grid based)
				all_lit.txt (whole simulation values)
			msprime_out (trees from msprime)
			slimtree_out (trees from slim)
			snamemake (simulation set up)
	flat_map (simulations with no heterogeneous landscape)
		normal (regular simulation)
			csv (all slim outputs)
				all_cor.txt (correlations)
				all_grid.txt (grid based)
				all_lit.txt (whole simulation values)
			msprime_out (trees from msprime)
			slimtree_out (trees from slim)
			snamemake (simulation set up)
		no_hair (no heritability)
			csv (all slim outputs)
				all_cor.txt (correlations)
				all_grid.txt (grid based)
				all_lit.txt (whole simulation values)
			msprime_out (trees from msprime)
			slimtree_out (trees from slim)
			snamemake (simulation set up)
		no_interaction (no interaction)
			csv (all slim outputs)
				all_cor.txt (correlations)
				all_grid.txt (grid based)
				all_lit.txt (whole simulation values)
			msprime_out (trees from msprime)
			slimtree_out (trees from slim)
			snamemake (simulation set up)
	interaction_grade (simulation with heterogeneous interaction rate)
		csv (all slim outputs)
				all_cor.txt (correlations)
				all_grid.txt (grid based)
				all_lit.txt (whole simulation values)
			msprime_out (trees from msprime)
			slimtree_out (trees from slim)
			snamemake (simulation set up)
			
Experiment2 (mutation rate was fixed for both species in each simulation, but allows the 
species to have different mutation effect sizes (and hence mutational variance)):
	mu_a (genetic architecture b)
		csv (all slim outputs)
				all_cor.txt (correlations)
				all_grid.txt (grid based)
				all_lit.txt (whole simulation values)
			msprime_out (trees from msprime)
			slimtree_out (trees from slim)
			snamemake (simulation set up)
	mu_b (genetic architecture c)
		csv (all slim outputs)
				all_cor.txt (correlations)
				all_grid.txt (grid based)
				all_lit.txt (whole simulation values)
			msprime_out (trees from msprime)
			slimtree_out (trees from slim)
			snamemake (simulation set up)
	mu_c (genetic architecture d)
		csv (all slim outputs)
				all_cor.txt (correlations)
				all_grid.txt (grid based)
				all_lit.txt (whole simulation values)
			msprime_out (trees from msprime)
			slimtree_out (trees from slim)
			snamemake (simulation set up)
	mu_d (genetic architecture e)
		csv (all slim outputs)
				all_cor.txt (correlations)
				all_grid.txt (grid based)
				all_lit.txt (whole simulation values)
			msprime_out (trees from msprime)
			slimtree_out (trees from slim)
			snamemake (simulation set up)

Experiment 3 (mutational variance was the same for newts and snakes in each simulation,
although polygenicity could be different (by varying mutation rate and mutation effect 
size)):
	v_a (genetic architecture f)
		csv (all slim outputs)
				all_cor.txt (correlations)
				all_grid.txt (grid based)
				all_lit.txt (whole simulation values)
			msprime_out (trees from msprime)
			slimtree_out (trees from slim)
			snamemake (simulation set up)
	v_b (genetic architecture g)
		csv (all slim outputs)
				all_cor.txt (correlations)
				all_grid.txt (grid based)
				all_lit.txt (whole simulation values)
			msprime_out (trees from msprime)
			slimtree_out (trees from slim)
			snamemake (simulation set up)
	v_c (genetic architecture h)
		csv (all slim outputs)
				all_cor.txt (correlations)
				all_grid.txt (grid based)
				all_lit.txt (whole simulation values)
			msprime_out (trees from msprime)
			slimtree_out (trees from slim)
			snamemake (simulation set up)
	v_d (genetic architecture i)
		csv (all slim outputs)
				all_cor.txt (correlations)
				all_grid.txt (grid based)
				all_lit.txt (whole simulation values)
			msprime_out (trees from msprime)
			slimtree_out (trees from slim)
			snamemake (simulation set up)

MSprime_script (holds the python script used to make msprime trees)
	snake_newt_gv_sm.py (python script used to make msprime trees)
	
R_scripts
	coevo_figs.html (figures for the paper)
	coevo_figs.Rmd (figures for the paper)
	figures (figure folder for the paper)
	gather_data.html (script to put the data into one file)
	gather_data.Rmd (script to put the data into one file)
	
SLiM_script
	cost_1on1_sm_cost_grade_newt_change.slim (slim script for newt only cost gradient)
	cost_1on1_sm_cost_grade_no_interaction_cost_grade.slim (slim script for newt and 
	snake cost gradient with no interaction)
	cost_1on1_sm_cost_grade_snake_change.slim (slim script for snake only cost gradient)
	cost_1on1_sm_cost_grade.slim (slim script for newt and snake cost gradient)
	cost_1on1_sm_final.slim (normal slim simulation)
	cost_1on1_sm_interrate_grade.slim (slim simulation with interaction rate gradient)
	cost_1on1_sm_no_hair.slim (slim simulation with no heritability)
	cost_1on1_sm_no_interaction.slim (slim simulation with no interaction)
	test_data (to play around with this slim simulation)
		gradent.png (gradent used in the simulation)
		test_1000_su_1e-10_nu_1e-10_sue_0.1_nue_0.1_.init.trees (a test msprime tree)

********************************************************************************

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

********************************************************************************

MIT License

Copyright (c) 2024 Victoria Caudill Peter Ralph

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 
