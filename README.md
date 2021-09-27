# Permuted_Ripley_K

This repositiory provides the R code for ***Statistical framework for studying the spatial architecture of the tumor immune microenvironment***. The code will illustrate the statistical methods outlined in Section 3. Briefly, 
  1. Observed and theoretical Ripley's K.
  2. Permute type labels across all cell locations (1000 times) and compute Ripley's K. Then average these values to use as a baseline for CSR (Permuted CSR).
  3. Cross validation for the cut points. Each patient is split into 10 folds and cut points are derived by different quantiles of the cell abundance and degree of spatial clustering (Observed Ripley's K - Permuted CSR) leading to five groups None (image with no positive cells), High Abundance & High clustering (HH), High Abundance & Low clustering (LH), Low Abundance & High clustering (LH), Low Abundance & Low clustering (LL). Cox proportional hazard (CPH) models are constructed using the cluster option for the repeated samples per subject for every combination of cell abundance and degree of spatial clustering and the effect of the five groups was measured by comparing the log liklihood of a model with onlu clinical features and a model with clinical and spatial features. The final model is built using the cutpoint combination that lead to the model where the spatial groups explained the most additional variation.   

Steps 1 and 2 were complete using the spatialTIME package which is avialible on both [CRAN](https://cran.r-project.org/web/packages/spatialTIME/index.html) and [Github](https://github.com/FridleyLab/spatialTIME). 
