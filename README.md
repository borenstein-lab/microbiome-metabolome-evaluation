# microbiome-metabolome-evaluation

This repository contains the analysis code for Noecker et al, 
"Defining and Evaluating Microbial Contributions to Metabolite Variation in Microbiome-Metabolome Association Studies" (submitted). 

**Contents:** 
- `FBA_functions.R`: A library of functions for processing and analyzing simulation data and performing analyses related to this work.
- `mainDatasetContributionsCorrelations.Rmd`: A notebook of analyses of the main simulation dataset described in the manuscript, including the code to generate Figures 1-4 and S1-S6.
- `variedSimulations_correlationResults.Rmd`: A notebook of analyses of the simulation datasets with varied nutrient inflow, including the code to generate Figures 5 and S9. 
- `LengthVaried_Simulations_correlationResults.Rmd`: Analysis of simulations of varying duration. Includes the code to generate Figure S7.
- `vmaxVaried_Simulations_correlationResults.Rmd`: Analysis of simulations with a varying maximum rate of reaction (vmax). Includes the code to generate Figure S8. 
- `ContributionDefinition_fluxInstantaneous.Rmd`: Analysis of the original ten-species dataset and simulations of longer duration using an alternative definition of contribution values.
- `HMPContributionsCorrelations.Rmd`: Analysis of simulations of 57 HMP-based communities. Includes the code to generate Figure 6 and Figure S10.
- `mimosaMainRun.Rmd`: A notebook documenting MIMOSA analysis of the simulation dataset, using both full model and the KEGG approximation. Includes the code to generate Figure 7. 

The input data files to re-run the analyses are available as Supplementary Data 1 and at http://borensteinlab.com/download.html. Save each sheet of the spreadsheet as a text file in the same directory. In order to run new simulations, you will need to separately download the pipeline code from http://borensteinlab.com/download.html. 
