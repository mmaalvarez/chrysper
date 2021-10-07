# chrysper

Generalized Linear Model (Negative Binomial or Poisson)

Regression: 
`geneKO read counts ~ time + offset`

Literature: Tzelepis et al., 2016 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5081405/ (mainly page 1196, and Figure 2D)


MAIN INTEREST
-------------
Develop a method to detect non-linear time trends for gene essentiality, i.e. identify whether the inactivation of a gene has a non-linear effect on cell fitness (convex: V or concave: /\ ) <-- BETTER HOCKEY STICK \\_  \_/  ¯\  /¯

For example, the loss of a gene that is needed for proliferation could cause an early cell count drop: \\_ , while in the case of a gene that is needed for things other than proliferation it could cause a late drop: ¯\ 


BACKGROUND
----------
"Gene essentiality": The effect that inactivating a single gene has on cell fitness. It can be inferred by the increase or decrease of cell counts after some days of selection, relative to the cell culture pool. Cell counts are proxied by sgRNA counts, which have to be normalized across samples

The DepMap database contains information on gene essentiality for many different tumor cell lines (from the Achilles, Score, and Demeter projects), but without specifying whether e.g. an essential gene causes an "early" or "late" cell count drop

This "early" or "late" cell count drop is explored in Tzelepis et al., 2016 for the HT29 cell line, but the idea here is to explore more cell lines and develop a method that efficiently tests for and detects these "hockey stick shapes"


METHODS
-------
This uses own and public longitudinal (time) CRISPR screening data, namely the untreated samples (so that the gene essentiality is non-conditional to very stressful conditions) from each cell line separately, and finds gene-KOs with significant non-linearity for `counts ~ time + offset`. Currently 'time' is categorized into t0, t-mid, and t-late, with the default orthogonal polynomial contrasts (linear and quadratic) <-- CHANGE CONTRASTS

- FILTER1: FDR<0.25 in time.Q (quadratic) for a cell line and a trend shape <-- CHANGE (see TODO)

- FILTER2: hit overlap across cell lines (mean FDR <0.25) for that trend shape <-- MAYBE NOT NECESSARY, hits found in at least 1 cancer cell line is fine

Draw predicted `counts ~ time` curves of gene hits

GO set enrichment of gene hits can provide more insights


NEXT STEP
---------
gene\*TP53\*time interactions

`counts ~ time*TP53 + offset`
