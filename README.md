# GPRTransform

Supplementary material for the GPR transformation method published in:

Machado, D., Herrgård, M. J., & Rocha, I. (2016). Stoichiometric Representation of Gene--Protein--Reaction Associations Leverages Constraint-based Analysis from Reaction to Gene-level Phenotype Prediction. PLOS Computational Biology 12(10): e1005140.

http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005140

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

Note that for compatibility reasons you should use *framed 0.2*. For any other uses of *framed* checkout the latest release with many improved features and better documentation. 
