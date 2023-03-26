Our manuscript results can be easily reproduced using script files. The workflow is divided into two parts: the first half focuses on obtaining a point estimate, while the second half deals with uncertainty quantification.

To achieve this, we define transformed-linear operations in a source file "2_TransformedOperations.R" and use the source code file "5_TPDM_Ang.R" to generate a data matrix 'X_t' through matrix multiplication. A source file "5_TPDM_Ang.R" returns the optimized vector 'b' and $\hat{X}$.

Subsequently, we use "6_TPDM_Est.R" to estimate the TPDM and "7_Ang_CPD.R" to decompose a prediction inner product matrix to obtain angles and angular masses. We then create the 95% joint polar region using "8_JointRegion.R".

To calculate an approximate conditional density given that $\hat{X}$ is large, we use "9_Ang_kde.R" to obtain an estimated kernel density, and "10_CondDens.R" to produce the density calculation. Finally, "11_CondInterval.R" is used to produce the 95% prediction interval.

Using these script files, users can easily reproduce our results and conduct their own analysis.
