This is preliminary software implementing stratified LD fourth moments regression (S-LD4M). In order to run this software, you will need to compute LD fourth moments and LD second moments (LD scores), and you will need unsigned summary statistics (chi^2 statistics, or per-normalized-genotype effect size estimates; estimates of polygenicity do not depend on the scaling of your summary statistics). We plan to make LD fourth moments available for download; alternatively, you may compute them yourself from a sequencing reference panel.

Package contents:

SLD4M.m

INPUT ARGUMENTS: chisq,  Mx1 matrix of chi^2 statistics; ell2, MxP matrix of LD scores (LD second moments), where P is #annotations; ell4, MxP vector of LD 4th moments; idxldsc, Mx1 vector of indices corresponding  to the reference SNPs that will be used in the regression; annot, MtotxP matrix of annotation values, where Mtot is the total number of reference SNPs (may be greater than M); weights_ldsc, Mx1 vector of weights for LDSC; weights_ld4m, Mx1 vector of weights for LD4M; no_blocks, number of jackknife blocks; ref_annot <= P', index of the annotation that will be used in denominator to compute enrichments; report_annot_mat, a PxP' matrix whose columns correspond to linear combinations of input annotations to output estimates for. For example, if P=2 with annotations for common and LF SNPs, then set report_annot_mat to [[1; 1], [1; 0], [0; 1]] and SLD4M will output Ma estimates for 3 annotations: all SNPs, common SNPs, LF SNPs. 

To download LD scores and LD 4th moments, go to this link: https://www.dropbox.com/sh/iiyftw01gdpt6un/AACU7AmWK45RxTmDJvRkdKhIa?dl=0

OUTPUT ARGUMENTS: maest, Estimated Ma for each category;  maerr, Ma standard error; h2est, heritability estimate; h2err,  heritability standard error; maenrichest, estimated Ma enrichment; maenricherr, Ma enrichment standard error; h2enrichest, estimated  heritability enrichment; h2enricherr, heritability enrichment standard  error; ma_jk, leave-one-out jackknife estimates of Ma; h2_jk, leave-one-out jackknife estimates of heritability; kurtexpest, kurtosis explained by functional annotations; kurtexperr, standard error; propkurtest, proportion of kurtosis explained by functional annotations (log scale); propkurterr, standard error.

NOTES: This implementation is parallelized across jackknife blocks, which is a large speed-up (approximately equal to the number of blocks), but it makes the code difficult to read and modify. Please contact me if you would like non-parallelized code.
