# chrysper

Generalized Linear Model (family Negative Binomial)

Regression: 
`geneKO read counts ~ time + offset`


MAIN INTEREST
-------------

Detect non-linear time trends for gene essentiality: this differs from classic gene essentiality analyses in that it could identify genes whose effect on fitness is non-monotonic (/\ or \/)

This uses untreated samples from each cell line separately, and finds gene-KOs with significant non-linearity for `counts ~ time + offset`. Currently 'time' is categorized into t0, t-mid, and t-late, with the default orthogonal polynomial contrasts (linear and quadratic)

	FILTER1: FDR<0.25 in time_cat.Q (quadratic)

	FILTER2: hit overlap across control cell lines (mean FDR <0.25)

	FILTER3: draw predicted counts ~ time curves of gene hits to check that they match across cell lines

GO set enrichment of gene hits can provide more insights


NOTES
-----

"Gene essentiality": The effect that inactivating a single gene has on cell fitness. It can be inferred by the increase or decrease of cell counts after some days of selection, relative to the cell culture pool. Cell counts are proxied by sgRNA counts, which have to be normalized across samples.


We want to detect PANCANCER non-linear gene essentiality:

	hits found across different cancer lines, with some of them even having altered genotypes, thus they overlap different cancer genetic backgrounds
	
	using only untreated (no A3A, no ATRi, no TMZ, normal O2%, no olaparib), so that the gene essentiality is non-conditional to very stressful conditions
