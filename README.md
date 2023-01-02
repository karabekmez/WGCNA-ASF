Since the WGCNA R-package is based only on the R2 value, matrix files  can be created in MATLAB to prepare the datasets for the hard-thresholding and ASF transformation.
INPUT: Expression Data
Calculates lowest hard-thresholding value (tau) satisfying scale-free criteria
Calculates lowest alpha value for transformation of correlation matrix by asymmetric sigmoid function (ASF) satidfying scale-free criteria (setting mu=tau)
OUTPUT: ".csv" files containing transformed similarity matrix by hard-thresholding and by ASF.Then, with using these matrix files as an input file, WGCNA tutorial can be applied by using ASF. 
