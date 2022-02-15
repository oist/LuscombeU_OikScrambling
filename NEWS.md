# OikScrambling 1.0.0

* Major change in the BSgenome packages, needing the removal of the old ones
  and the recomputation of the data.  The new packages allow to much easier
  link the GBreaks objects to their genome sequence or annotations.

# OikScrambling 0.8.3

* Move some regions of interest in a separate vignette.

# OikScrambling 0.8.2

* Added a new `SegmentMovements` vignette.

# OikScrambling 0.8.1

* Added a new `wgo` whole-genome object indicating if a base is part of a
  syntenic block, or just aligned, or unaligned.

# OikScrambling 0.8.0

* Added `nonCoa` flags to the `gbs`, `unal` and `coa` objects.

# OikScrambling 0.7.1

* Created compound plots of width distribution for Figure 2.

# OikScrambling 0.7.0

* Overlap with repeats on both strands.  The difference is not so massive,
  considering that the unaligned and unmapped regions are already strandless
  anyway.

# OikScrambling 0.6.6

* New functions `OikScrambling:::compDistClass()` and
  `OikScrambling:::compDistClass()`, and new option `short` for
  `OikScrambling:::compDistance()`, to produce more compact summary graphs.

# OikScrambling 0.6.5

* Identify the flipped repeat as MITE and search for other copies.

# OikScrambling 0.6.4

* Add regions of interest with repeat flipping and insertion.

# OikScrambling 0.6.3

* Add Hox regions of interest for Figure 1

# OikScrambling 0.6.2

* Add 10-Mbp dot plots for Figure 1

# OikScrambling 0.6.1

* Add Mike's vignette on coalescing cutoffs.

# OikScrambling 0.6.0

* Adding _Drosophila_ data for a plot on CP's MBSJ2021 poster.

# OikScrambling 0.5.0

* `load_one_to_ones()` function moved from _GenomicBreaks_ to this package.
* Internal functions updated to load _C. intestinalis_ annotations.

# OikScrambling 0.4.0

* From now on the `main` branch must remain buildable.
* Added a `NEWS.md` file to track changes to the package.

# OikScrambling before that

* Vignettes to build with `pkgdown`
* Plenty of other things

# Before OikScrambling

* OikScrambling stems from the old version of the vignettes
  of the `GenomicBreaks` package.
