# A case for multi-gear assessments: detection probabilities of nearshore fish with eDNA and seine nets vary by functional traits
## Overview
Repository contains data analyses for manuscript in prep using data collected in southern Southeast Alaska in 2021 and 2022. 
This repo uses seine and environmental DNA data. eDNA data were generated from kjledger's [nearshore_eDNA repo](https://github.com/kjledger-NOAA/nearshore_eDNA)
and then pulled into this repo for a combined analysis. 

## Scripts in repo

### 1. wholecomm_mv_analysis.rmd

Script that joins data types, cleans, and performs multivariate betadisper/permanova tests, generates PCoA, and identifies indicator species (ISA)
Creates table 1, table S4, figures 3, 4, S5, S6

### 2. spOccupancy_data_cleaning.R 

Cleans data in preparation for running a species occupancy model for spOccupancy

### 3. wc_spoccupancy.rmd

Runs joint species occupancy model and post-hoc analysis
Creates figure 5

### 4. Map_figure.rmd

Creates figure 2 

### ppc_integratedMsPGOcc.R

function sourced in wc_spoccupancy script to run the posterior predictive checks 
