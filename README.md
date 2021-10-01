# chrysper

Generalized Linear Model (Negative Binomial or Poisson)

Regression: 
`geneKO read counts ~ time + offset`


MAIN INTEREST
-------------

Detect non-linear time trends for gene essentiality: this differs from classic gene essentiality analyses in that it could identify genes whose effect on fitness is non-monotonic (convex: V or concave: /\ ) <-- BETTER HOCKEY STICK \\_ /¯

This uses untreated samples from each cell line separately, and finds gene-KOs with significant non-linearity for `counts ~ time + offset`. Currently 'time' is categorized into t0, t-mid, and t-late, with the default orthogonal polynomial contrasts (linear and quadratic) <-- CHANGE CONTRASTS

- FILTER1: FDR<0.25 in time.Q (quadratic) for a cell line and a trend shape

- FILTER2: hit overlap across cell lines (mean FDR <0.25) for that trend shape <-- MAYBE NOT NECESSARY, hits found in at least 1 cancer cell line is fine

Draw predicted `counts ~ time` curves of gene hits

GO set enrichment of gene hits can provide more insights


NOTES
-----

"Gene essentiality": The effect that inactivating a single gene has on cell fitness. It can be inferred by the increase or decrease of cell counts after some days of selection, relative to the cell culture pool. Cell counts are proxied by sgRNA counts, which have to be normalized across samples.

Using only untreated (no A3A, no ATRi, no TMZ, normal O2%, no olaparib), so that the gene essentiality is non-conditional to very stressful conditions
