# RTE2016

This repository contains the data and code to accompany the manuscript: 

Katie L. Barott, Ariana S. Huffmyer, Jennifer M. Davidson, Elizabeth A. Lenz, Shayle B. Matsuda, Joshua R. Hancock, Teegan Innis, Crawford Drury, Hollie M. Putnam, Ruth D. Gates (2020). Bleaching resistant corals retain heat tolerance following acclimatization to environmentally distinct reefs. bioRxiv doi: 10.1101/2020.09.25.314203

**Abstract:** Urgent action is needed to prevent the demise of coral reefs as the climate crisis leads to an increasingly warmer and more acidic ocean. Propagating climate change resistant corals to restore degraded reefs is one promising strategy; however, empirical evidence is needed to determine if resistance is retained following transplantation within or beyond a coral’s natal reef. Here we assessed the performance of bleaching-resistant individuals of two coral species following reciprocal transplantation between environmentally distinct reefs (low vs high diel variability) to determine if stress resistance is retained following transplantation. Critically, transplantation to either environment had no influence on coral bleaching resistance, indicating that this trait was relatively fixed and is thus a useful metric for selecting corals for reef restoration within their native range. In contrast, growth was highly plastic, and native performance was not predictive of performance in the novel environment. Coral metabolism was also plastic, with cross transplants of both species matching the performance of native corals at both reefs within three months. Coral physiology (autotrophy, heterotrophy, and metabolism) and overall fitness (survival, growth, and reproduction) were higher at the reef with higher flow and fluctuations in diel pH and dissolved oxygen, and did not differ between native corals and cross-transplants. Conversely, cross-transplants at the low-variability reef had higher fitness than native corals, thus increasing overall fitness of the recipient population. This experiment was conducted during a non-bleaching year, which suggests that introduction of these bleaching-resistant individuals will provide even greater fitness benefits to recipient populations during bleaching years. In summary, this study demonstrates that propagating and transplanting bleaching-resistant corals can elevate the resistance of coral populations to ocean warming while simultaneously maintaining reef function as the climate crisis worsens.

**Repository contents:**

**Data/:**

Coco.csv: CTD environmental data for inner lagoon (PR4)

Monty.csv: CTD environmental data for outer lagoon (PR13)

PAR.csv: Summary Odyssey light logger data for both inner lagoon and outer lagoon

RTE_Frags.csv: Response variables for each fragment involved in the reciprocal transplant

RTE_Genotypes.csv: Averaged response variables for all fragments per genotype involved in the reciprocal transplant

RTE Metadata.csv: Column descriptions for "RTE_Frags.csv" and "RTE_Genotypes.csv"

seaphox_final.csv: SeapHOx environmental data for inner lagoon and outer lagoon (PR4, PR13)

SpawnDynamic.csv: Data for MCAP spawning activity

Stress_Master.csv: Response variables for each fragment involved in the thermal challenge

**Analysis/:**

MixedModelAnalysis.R: Linear mixed effect models for reciprocal transplant

Odyssey.R: Analysis of PAR at each lagoon site

ResponseSummaries.R: Code to generate summary figures

SeapHOxAnalysis.R: Analysis of environmental regimes at each lagoon site

StressModelAnalysis.R: Linear mixed effect models for thermal challenge

StressTest.R: Code to generate thermal challenge figures
