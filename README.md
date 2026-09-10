# Sterilization of *Amyelois transitella* by electron beam technology

<https://github.com/ChuckBV/y25_lars_laber_ebeam_stats>

This repository contains data sets and scripts for summary, analysis and 
presentation from an experiment conducted in 2024 by Lars Laber and Houston 
Wilson (University of California Riverside) and Chuck Burks (USDA-ARS, Parlier). 
The objective of these experiments was to find a sterilizing dose for 
*A. transitella* (navel orangeworm, NOW) using electron beam, and to determine 
the impact of substerile and sterilizing doses on NOW vigor as determined by 
flight cylinder and longevity assays.

This study involved a total of four cohorts that were taken from Parlier, CA 
approximately 150 miles to Fremont (near San Jose, CA) for irradiation at the 
Steri-tek facility there. Two of these cohorts were used for a flight cylinder 
and longevity experiment, and two were used in studies of impact on fertility. 
Details are given in a draft manuscript in the /doc sub directory.

## Subdirectories:

- /doc – Overview. Contains an early draft manuscript, and PDF slide deck
with figures and analysis of deviance tables for experiment 1 (flight cylinder
escapes and proportion of months surviving > 14 days) and experiment 2 
(proportion mated, eggs per mated pair, and count of F1 pupae  per mated pair.)
- /data – contains the data files used, as R data sets (*.Rds) and as 
comma-delimitted text (*.csv)
- /output – figures and summary csv files
- /scripts – R scripts used to access /data files and create /output file

## Scripts

### Initial
 - themes.R -- Tweaked themes for ggplot
 - script1_data_files.R -- provides a readRDS for each data set, and a 
short data description

### Experiment 1 -- Flight cylinders and longevity
 - script2\*.R -- Analyze flight cylinder data.
 - script3\*.R -- Analyze longevity data.

### Experiment 2 – Fertility assay
 - script4*.R -- Analyze mating success
 - script5*.R -- Analyzing total eggs per mated pair
 - script6*.R -- Analyzing F1 pupae per mated pair

### Finishing and supporting material
 - script7*.R -- All ggplots in presentation form
 - citations.R -- Bibliographic info for R packages, generated from citation()
 
 