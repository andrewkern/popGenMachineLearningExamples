# popGenMachineLearningExamples

This repository is meant to house a series of jupyter notebooks that showcase some simple examples of using supervised machine learning for population genetics inference. 

The first notebook that we have added `demographicModelSelectionExample.ipynb` is meant as a companion to our recent review-- Schrider and Kern (2017) "Machine learning for population genetics: a new paradigm." That paper can be found on bioRxiv here: https://www.biorxiv.org/content/early/2017/10/20/206482

The subsequent notebooks really build on ideas presented in the first one. We aim to present a diversity of applications within population genetics that highlight a number of algorithms and ML practices. We recommend going through these in something like the following order:
1. `demographicModelSelectionExample.ipynb`- using a Random Forest classifier for doing demographic model selection
2. `sweepDetectionExample.ipynb` - using a Support Vector Machine classifier for detecting & categorizing sweeps
3. `ancPopSizeRegressionExample.ipynb` - using Random Forest regression to estimate past population size
