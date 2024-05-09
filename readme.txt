# Title of the document
A meta-analysis on alternative mating tactics:  when the main and the alternative provide similar reproductive success

# Authors
da Silva, A.V. (mandyvieira15@gmail.com)
Almeida, J.G.L. (jgabriel.bio2@gmail.com)
Ventura, S.P. (stefania.ventura2@gmail.com)
Oliveira, R. (reislaxoliveira@gmail.com)
Peixoto, P.E.C. (pauloenrique@gmail.com) * corresponding author

# GENERAL INFORMATION #

### About the data
Because this review is a meta-analysis, 
it contains information from different original studies. Thus, the data were collected in different dates and different locations around the world

### Keywords
sexual selection
alternative mating tactics
alternative reproductive tactics
fixed tactics
flexible tactics
quantitative review

### Language
All data were written in English

### Funding
- CNPq
- Capes


# DATA AND FILE OVERVIEW #

### File: data_tactics.csv

#### General description
This file contains information about the reproductive success of species exhibiting alternative mating tactics and whether individuals of each species can change this tactics (flexible) or not (fixed) during their lifespan

#### Date of creation
September, 2019

#### Last update
September, 2022

### File: diag.xlsx

#### General description
This file contains the variance-covariance matrices used in meta-analytical models

#### Date of creation
February, 2024

#### Last update
May, 2024

File: analysis.R
General description
This file contains the script to generate all meta-analytical models

Date of creation
September, 2020

Last update
May, 2024

File: figure.R
General description
This file contains the script to generate the meta-analytical figures. It is a modification from orchard_plot function (from orchaRd v. 0.0.9 package) to remove the prediction interval

Date of creation
September, 2020

Last update
May, 2024

# METHODOLOGICAL INFORMATION #

### Brief description of methods
We conducted a search on Google Scholar, Scopus and Web of Science (core colletion) to find for studies published between 1945- September 2022. For Google Scholar, we recorded all studies up to the 100th page. We used a combination of the following keywords: (1) “alternative reproductive tactics AND mating success”; (2) “alternative reproductive strategies AND mating success”; (3) “alternative reproductive tactics AND reproductive success”; (4) “alternative reproductive strategies AND reproductive success”; (5) “alternative mating tactics AND mating success”; (6) alternative mating strategies AND mating success; (7) “alternative mating tactics AND reproductive success”; (8)” alternative mating strategies AND  reproductive success”.
We included int our dataset studies containing information about the reproductive success of species expressing alternative mating tactics. Whenever as possible, we investigated if individuals of species could alternate between tactics during its lifespan. Using measures of reproductive success provided in the studies, we calculated an overall and comparable metric of effect size (Hedges'g).


# DATA-SPECIFIC INFORMATION #

### Number of variables: 40

### Number of rows: 115 (including the header)

### Description of each variable included in the file
- n_obs: ID of each observation included
- id: study ID
- researcher: person responsible for collecting the data from that study
- author: author(s) of the study
- year: year of the publication of that study
- species: species used in the study (as described by the authors)
- species2: species used in Tree Web of Life to build the phylogenetic tree (some species has changed their name since the study was published or some species ware not found on Tree Web of Life. In both cases, we needed to change for a synonim or for the closest related species name)
- family: family of the species
- order: order of the species
- vertebrates: whether the species is a vertebrate or invertebrate
- experiment: if the study were conducted under controlled conditions (lab), or if it was observational (fieldwork) or if they conducted an experiment on the field (fieldwork experimental)
- type of tactic: if individuals of that species can change between tactics during its lifespan. Individuals adopting flexible tactics can change tactics while individuals adopting fixed tactics cannot change. For some species, we did not achieve the type of tactics adopted
- tactic_A: the main tactic
- tactic_B: the secondary tactic
- comparison: category to group similar tactics among species
- comparison2: we classified comparison into two different groups depending if the tactics was related to morphology or behavior
- difference: if there is any difference between individuals adopting main or secondary tactics
- type_difference: description of the difference when there is a difference between individuals adopting main or secondary tactics
-measure_1: how the study evaluated the reproductive success. We grouped in mating (e.g. duration or number of matings) or genetics (e.g. number of offspring)
- measure_2: how study evaluated the reproductive success as described by the study
- measure_timing: when the study evaluated the reproductive success (before, during or after mating)
- mean_sucA: mean of the reproductive success of males adopting the main tactic
- mean_sucB: mean of the reproductive success of males adopting the secondary tactic
- prop_sucA: proportion of males adopting the main tactic that mated
- prop_sucB: proportion of males adopting the secondary tactic that mated
- prop_failA: proportion of males adopting the main tactic that did not mate
- prop_failB: proportion of males adopting the secondary tactic that did not mate
- n_sucA: number of males adopting the main tactic that mated
- n_sucB: number of males adopting the secondary tactic that mated
- n_failA: number of males adopting the main tactic that did not mate
- n_failB: number of males adopting the secondary tactic that did not mate
- sd_A: standard deviation of the mean of the main tactic
- sd_B: standard deviation of the mean of the secondary tactic
- se_A: standard error of the mean of the main tactic
- se_B: standard error of the mean of the secondary tactic
- n_A: number of individuals adopting the main tactic
- n_B: number of individuals adopting the secondary tactic
- yi: Hedges'g effect size
- vi: variation of Hedges'g effect size
- measure: type of effect size calculated (SMD for Hedges'g)

# Missing data
Recorded as a blank cell
