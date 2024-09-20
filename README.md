# Study of human cancer incidence as a function of age using mathematical models
This repository contains the source code of my final project for the MB151P133 Mathematical modelling in bioinformatics course at the Faculty of Science at Charles University in Prague.
It is based on an article <em>A Poisson distribution-based general model of cancer rates and a cancer risk-dependent theory of aging</em> available [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10522393/).

The main body of the project is available [here](https://github.com/bn-codes/mmb_project/blob/main/Projekt.pdf) (only in Czech).

## Structure of the Repository

* `Figures/` contains code that reproduces the graphs of Figure 3 from the original article
* `Leukemia/` contains code that tests the models on different data (here data on leukemia incidence in the population)
* `Models/` contains code that reproduces the models from the original article
* `Plots/` contains all the plots obtained from reproducing and testing the models
* `vysledky_final.csv` contains the results obtained from the Models/adaptive_model.py

## Abstract

The aim of this project is to reproduce the results of a model predicting age-related cancer incidence. The model uses a Poisson probability distribution to estimate the time 
of cancer onset and works with other parameters such as the number of cell turnovers per year and the rate of accumulation of deleterious mutations in the genetic information 
of cells. The replicated model is extended and analyzed to predict the onset of leukemia. Analysis of the models has provided information on age-dependent cancer prediction 
possibilities with the need for further steps for a better understanding of the problem.
