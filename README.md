# chrysper

Generalized Linear Model (family Negative Binomial)

Regression: 
`geneKO read counts ~ treatment*time + control variables + offset + ε`


MAIN INTEREST
-------------

Detect non-linear time trends for the effect that inactivating a single gene has on cell fitness (i.e. increase or decrease of cell -sgRNA- counts)

	a) For controls and treated samples separately, find geneKOs with significant non-linearity for `counts ~ time + offset` (>2 time categories with orthogonal polynomial contrasts)

		- E.g. time.Q (quadratic) or time.C (cubic) significant

	b) For all samples together, and the previous gene hits, check `counts ~ treatment * time + offset` with backward contrasts to see whether the non-linearity above depends on presence/absence of treatment, at least at 1 time point (i.e. gene×time×treatment epistasis)

		- E.g. time_mid-vs-0:treatment.L and/or time_late-vs-mid:treatment.L significant


GO set enrichment of gene hits can provide more insight



Tumor treatment improvement
---------------------------

Find genes whose KO makes the treated samples counts to decrease more gradually, so the treatment may not be too aggressive for the patient

	- ideally, for the controls the line should be flat (βs = 0), so the loss of these genes per se does NOT affect fitness

	- approaches to do so:

		i) in an absolute manner, i.e. defining β and FDR thresholds

		ii) relative to the "default" shape, i.e. more "straight" and gradual than the default effect of drug without gene function loss

			- 1st check how tumor cell count decreases just due to treatment (e.g. A3A presence), without the influence of the loss of a gene (expectedly sudden/fast). This could be proxied by the shape of *treated* sample counts across time when using

				a) 350 non-essential genes only

				b) all genes together?

	- validations that a hit is real:

		i) it replicates across cell lines (not NECESSARY though, it could really have this effect only in one cell line, due to genetic background, but if it's replicated across cell lines, it is proof enough that it is a real hit)

		ii) there are pathways/GO terms enriched for these hits -- ideally, the mechanism should make sense (e.g. related to A3A - deamination - BER - DSBs..., or others)
