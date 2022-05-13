# OikScrambling 4.0.0

* Reformed the `wgo` object to match the 4 classes: _isolated alignment_,
  _breakpoint region_, _collinear alignment_ and _bridge region_.

# OikScrambling 3.0.6

* Alignment features plotted on operon boundaries for Figure 3.
* New chromosome plots for Figure 4.
* New panel showing that short bridge regions are mostly made of short introns.

# OikScrambling 3.0.5

* Explore strand proportion indices.
* New region width plots for Figure 2.

# OikScrambling 3.0.4

* Draft chromosome plots for Figure 4 (ex-3, ex-2b).
* Natural scale width plots for Figure 2
* 5-panel feature plots for Figure 2

# OikScrambling 3.0.3

* Mid-centered feature plots for Figure 2.

# OikScrambling 3.0.2

* Plot genomic features on white background for Figure 2.
* Compute strand randomisation indexes.
* Provide alignment information for manual brush-up of Hox panels in Figure 1.
* Add Osaka–Aomori Oxford plots.

# OikScrambling 3.0.1

* Plot operons density near alignment boundaries.
* Compute confidence intervals on width distributions.
* Plot width distributions in PDF format.

# OikScrambling 3.0.0

* Corrected arm names on query genomes.
* New functions `isSyntenic()`, `isSynbrachial()`,
  `nameSyntenic()` and `nameSynbrachial()`.
* Renamed _mapped unaligned_ regions to _bridge_ regions and renamed `unalMap`
  object to `bri`.
* Added an 'isolated' category in the `RegionWidths` vignette.

# OikScrambling 2.2.2

* Parallel plot of chr2
* Export aggregated numbers about segment widths.

# OikScrambling 2.2.1

* New Oki-Osa Oxford plots for Figure 1.

# OikScrambling 2.2.0

* Strand-colored version of the 10-Mb window Oxford plots.  Needs
  `GenomicBreaks` version `0.13.1` or superior.
* Human-mouse 10-Mb window Oxford plot.  Needs `BreakpointsData`
  version `3.9.0` or superior.

# OikScrambling 2.1.0

* Add a vignette about the movement of the PAC3 gene.

# OikScrambling 2.0.0

* Adopt semantic versionning: any backwards-incompatible change increases
  the major version number, no matter how (un)important it is.

* Use the new `GenomicBreaks::bridgeRegions()` function that produces
  zero-width ranges instead of the local function that artificially
  added 1nt to the ends.

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
