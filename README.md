# GPRTransform

Supplementary material for the GPR transformation method published in:

Machado, D., Herrgård, M. J., & Rocha, I. (2016). Stoichiometric Representation of Gene--Protein--Reaction Associations Leverages Constraint-based Analysis from Reaction to Gene-level Phenotype Prediction. (submitted)

It contains an implementation of the method that transforms an SBML model by integrating the GPR associations directly into the stoichiometric matrix. This enables gene-based analysis using several constraint-based methods.

This repository contains all the code, models and results used in the paper, including examples of application of the GPR transformation to:

+ Phenotype prediction
+ Random flux sampling
+ Gene essentiality analysis
+ Omics data integration
+ Elementary Mode Analysis
+ Rational strain design

It also contains implementations of new methods proposed in the paper such as:

+ gene-pFBA
+ gene-MOMA, gene-lMOMA
+ gene-GIMME
+ gene-EFlux

Dependencies:

+ framed v0.2: https://github.com/cdanielmachado/framed/releases/tag/v0.2
