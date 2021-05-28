# community-simulator

## Introduction
This repository contains the raw data files and code for the paper

Juan Diaz-Colunga, Nanxi Lu, Alicia Sanchez-Gorostiaga, Chang-Yu Chang, Helen S. Cai, Joshua E. Goldford,
Mikhail Tikhonov and Alvaro Sanchez (2021). Top-down and bottom-up cohesiveness in microbial community coalescence. *bioRxiv*:XXXXXX. 

The code includes R files to analyze the data and generate plots, as well as python scripts to run the simulations.
This repository also includes an updated version of [the Community Simulator package](https://github.com/Emergent-Behaviors-in-Biology/community-simulator)
with new functionalities described below.

## Data files and analysis
The ./data folder contains three files from the experiments:
- metaData.csv: includes the sample name, community numbers, supplied carbon source and other metadata for each experiment ID. The arrangement of the
96 well plates where coalescence experiments were carried out is illustrated in
the file Coalescence_setup.pptx.
- OTU_table.csv: contains the sequencing data for each experiment ID.
- taxonomy_silva.csv: contains the OTU-to-species (or families) map.

In addition, it contains a simul_data.txt with the results of the simulations (see below). The folder also contains a set of R scripts:
- analyzeData.R: this is the primary script for data analysis.
It loads the files mentioned above and runs the analyses described in the paper
(the R file is also intensively commented itself so it can be followed easily).
It generates a series of plots that are saved in the ./plots folder.
- myFunctions.R: contains some auxiliary functions called during the execution of
analyzeData.R, e.g. a function that returns the Bray-Curtis similarity of two vectors.
- isolates_abundance_from_ESVs.R: this file is an attempt at deconvoluting
Exact Sequence Variant abundances into actual spepcies abundances, given that
species can carry multiple variants of their 16S rRNA gene and that some overlap
across species is possible. Although the file is sourced in analyzeData.R,
the deconvolution was finally not done and the ESV-to-species mapping was done
one-to-one as described in the paper.

## Simulations
The file coalescence.py runs the simulations and generates the simul_data.txt file
for downstream analysis.

Simulations were carried out using
[the Community Simulator package](https://github.com/Emergent-Behaviors-in-Biology/community-simulator).
New features were added. First, we now give the option for individual species to
have unique metabolic matrices instead of a shared one. This behavior can be
controlled by passing a new entry to the dictionary of modeling assumptions
of the Community Simulator (see the package [README](./community-simulator/README.md) for details):
```
assumptions['metabolism'] = 'common'
```
or
```
assumptions['metabolism'] = 'specific'
```
The first choice enables the default behavior of the original Community Simulator.
The second makes it so the metabolic matrix is a list of matrices with as
many elements as species in the model. Each individual matrix is generated
randomly by calling the `MakeMatrices()` function of the original package.
This is equivalent to the metabolic matrix being three- instead of two-dimensional,
with the aditional dimension corresponding to the species index.

We also added a new control parameter (not used in the publication but
available in the updated package files) that can be passed when
the user chooses a species-specific metabolic matrix.
This parameter, that we denote *r<sub>s</sub>*, controls the correlation
between a species' resource preference and its secretion fluxes.
Setting
```
assumptions['rs'] = 0
```
will make it so species' secretions are independento of their resource
preferences. On the other hand, choosing
```
assumptions['rs'] = 1
```
will make it so species can only secrete byproducts that they can utilize themselves,
with secretion fluxes being biased towards higher values the larger the
species' preference for the resource (*c<sub>i&alpha;</sub>* for species _i_
preference for resource _&alpha;_).
Values between 0 and 1 are also accepted for intermediate behaviors.

The updated package files together with the original instructions to install
and operate it are in the ./community-simulator folder in this repository.

## Minimal model
The file coalescence_minimalModel.py runs the simulations corresponding to the
minimal model. It outputs a set of figures that are saved
in the ./minimal-model folder. Several parameter choices can be made to
control the behavior of the minimal model as described in the paper.