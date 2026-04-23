# A case for multi-method assessments: detection probabilities of nearshore fish with eDNA and seine nets vary by functional traits
## Overview
Repository contains data analyses for manuscript accepted at CJFAS using data collected in southern Southeast Alaska in 2021 and 2022. 
This repo uses seine and environmental DNA data. eDNA data were generated from kjledger's [nearshore_eDNA repo](https://github.com/kjledger-NOAA/nearshore_eDNA)
and then pulled into this repo for a combined analysis. 

## Manuscript abstract
Ecological studies aim to understand species distributions, yet the sampling methods affect which species are detected and may be influenced by species traits. Detection of marine fishes based on species traits have been studied using traditional sampling methods, but such studies have generally not extended to environmental DNA (eDNA). Here, we investigated which functional species traits (scale type, schooling behavior, and position-in-water-column) are detected by eDNA metabarcoding and beach seines in nearshore eelgrass, mixed eelgrass, and understory kelp habitats. Using data from 35 sites across southeast Alaska, we applied occupancy modeling to estimate detection of species traits by each method. Detection probability with eDNA was 27 times greater for the species with deciduous scales (*Clupea pallasii*) compared to species with non-deciduous scales, and lower for species with plates (rather than scales). Conversely, species with plates showed greater odds of detection with beach seines. Given the novelty of eDNA sampling, quantifying interactions between detection and functional traits will be important to accurately characterize species distributions across marine habitats. 

## Scripts in repo

### 1. wholecomm_mv_analysis.rmd

Script that joins data types, cleans, and performs multivariate betadisper/permanova tests, generates PCoA, and identifies indicator species (ISA).
Creates table 1, table S4, figures 3, 4, S5, S6

### 2. spOccupancy_data_cleaning.R 

Cleans data in preparation for running a species occupancy model for spOccupancy

### 3. wc_spoccupancy.rmd

Runs joint species occupancy model and post-hoc analysis.
Creates figure 5

### 4. Map_figure.rmd

Creates figure 2 

### ppc_integratedMsPGOcc.R

function sourced in wc_spoccupancy script to run the posterior predictive checks 
