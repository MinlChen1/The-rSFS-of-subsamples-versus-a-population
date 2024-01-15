#run R script
NITER=100
PATHR=`realpath`
R --vanilla --args 1 100 $NITER $PATHR < ./Sweep_SFSline_diversity_sims.R #calculation for haplotypes
R --vanilla --args 2 100 $NITER $PATHR < ./Sweep_SFSline_diversity_sims.R #calculation for diploid individuals
