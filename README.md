# OikScrambling

This package contains a collection of vignettes to perform the analysis of
pairwise genome alignments between _Oikopleura_ genomes.

The core functions used here are maintained in our _GenomicBreaks_ R package,
which is fully documented at: <https://oist.github.io/GenomicBreaks>.

## How to build this site

See the main vignette for details

    /bucket/LuscombeU/common/Singularity/GenomicBreaks-0.10.0.sif Rscript -e "remotes::install_github('oist/GenomicBreaks', repos=BiocManager::repositories())"
    /bucket/LuscombeU/common/Singularity/GenomicBreaks-0.10.0.sif Rscript -e "install.packages('BreakpointsData', repos='https://oist.github.io/plessy_oikgenomes_drat/')"
    # And so on for the BSgenome packages.
    /bucket/LuscombeU/common/Singularity/GenomicBreaks-0.10.0.sif R CMD INSTALL .
    srun -pcompute -c8 --mem 300G --pty /bucket/LuscombeU/common/Singularity/GenomicBreaks-0.10.0.sif Rscript -e "pkgdown::build_article('LoadGenomicBreaks')"
    srun -pcompute -c8 --mem 300G --pty /bucket/LuscombeU/common/Singularity/GenomicBreaks-0.10.0.sif Rscript -e "pkgdown::build_site(devel=TRUE, lazy=TRUE)"
    srun -pcompute -c8 --mem 300G --pty /bucket/LuscombeU/common/Singularity/GenomicBreaks-0.10.0.sif Rscript -e "pkgdown::build_site()"
