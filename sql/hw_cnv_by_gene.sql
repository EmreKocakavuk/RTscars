/*
---
Created: April 30, 2020
Updated: May 7, 2020
Author: FP Barthel
---
*** MODIFIED FROM `cnv_by_gene_gatk.sql` from GLASS (original by FP Barthel) ***
This pull copy number from Hartwig
Updated to ensure homdels of CDKN2A
---
*/
WITH
selected_genes AS
(
	SELECT dr.gene_symbol,chrom,pos
	FROM ref.driver_genes dr
	LEFT JOIN ref.genes ge ON ge.gene_symbol = dr.gene_symbol
	WHERE has_cnv IS TRUE
),
gene_seg_intersect AS
(
    SELECT sampleid, gene_symbol, gs.chrom, (upper(t0.pos * gs.pos) - lower(t0.pos * gs.pos) -1) AS w, major_allele_ploidy
    FROM hwc gs
    INNER JOIN selected_genes t0 ON t0.chrom = gs.chrom AND t0.pos && gs.pos
),
gene_sample_call AS
(
    SELECT sampleid, gene_symbol, 
		sum(w * major_allele_ploidy::numeric) / sum(w) AS wploidy
    FROM gene_seg_intersect
    GROUP BY sampleid, gene_symbol
)
SELECT
	gc.sampleid,
	gc.gene_symbol,
	(CASE
	 WHEN gc.wploidy >= 0.5 AND gc.wploidy <= 1.5 THEN 0
	 WHEN gc.wploidy < 0.5 THEN -1
	 WHEN gc.wploidy > 1.5 THEN 1
	 ELSE NULL
	 END) ploidy_based_cn,
	gc.wploidy
FROM gene_sample_call gc