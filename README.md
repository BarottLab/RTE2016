# RTE2016

This repository contains the data and code to accompany the manuscript: 

Katie L. Barott, Ariana S. Huffmyer, Jennifer M. Davidson, Elizabeth A. Lenz, Shayle B. Matsuda, Joshua R. Hancock, Teegan Innis, Crawford Drury, Hollie M. Putnam, Ruth D. Gates (2020). Coral bleaching response is unaltered following acclimatization to reefs with distinct environmental conditions 

**Repository contents:**

**Data/:**

PAR.csv: Summary Odyssey light logger data for both inner lagoon and outer lagoon

RTE_Frags.csv: Response variables for each fragment involved in the reciprocal transplant

RTE_Genotypes.csv: Averaged response variables for all fragments per genotype involved in the reciprocal transplant

RTE Metadata.csv: Column descriptions for "RTE_Frags.csv" and "RTE_Genotypes.csv"

EnvData.csv: SeapHOx environmental data for inner lagoon and outer lagoon (PR4, PR13)

SpawnDynamic.csv: Data for MCAP spawning activity

Stress_Master.csv: Response variables for each fragment involved in the thermal challenge

flow.csv: Flow regime data for inner lagoon and outer lagoon (PR4, PR13)

sedimentation.csv: Sedimentation rate data for inner lagoon and outer lagoon (PR4, PR13)

**Analysis/:**

MixedModelAnalysis.R: Linear mixed effect models for reciprocal transplant

Odyssey.R: Analysis of PAR at each lagoon site

ResponseSummaries.R: Code to generate summary figures

EnvSummary.R: Analysis of environmental regimes at each lagoon site

StressModelAnalysis.R: Linear mixed effect models for thermal challenge

StressTest.R: Code to generate thermal challenge figures

flow.R: Code to analyze flow dynamics at each lagoon site

sedimentation.R: Code to analyze sedimentation rates at each lagoon site

NMDS.R: Code for multivariate analyses
