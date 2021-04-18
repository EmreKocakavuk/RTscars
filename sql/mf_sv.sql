/*
---
Date: April 21, 2020
Author: FP Barthel
---
Determine the mutation burden by variant type and fraction for each case in the RT-scars cohort
- Comparable analysis to `mf_snv_indel.sql`
- Limited to SVs called by LUMPY
- Coverage adjustment is used and a threshold of 10x coverage is used
- Genotype of 0/1 or 1/1 is required for variant calling
- Genotypes that failed to call (./.) were excluded
---
*/
WITH
svtypes (svtype) AS
(
	VALUES ('BND'),('DEL'),('DUP'),('INV')
),
paired_geno_sv AS
(
	SELECT DISTINCT sp.*, ta.eventid, ta.svtype, ta.genotype geno_a, tb.genotype geno_b--, (SELECT bool_or(ti.chrom = cap.chrom AND cap.pos && ti.pos) FROM ref.intersect_exome_hg19 cap)
	FROM analysis.scars_pairs sp
	INNER JOIN svg ta ON ta.aliquot_barcode = sp.tumor_barcode_a
	INNER JOIN svg tb ON tb.aliquot_barcode = sp.tumor_barcode_b AND tb.eventid = ta.eventid
	INNER JOIN svi ti ON ti.case_barcode = sp.case_barcode AND ti.eventid = tb.eventid --AND ti.eventid = ta.eventid
	--INNER JOIN ref.intersect_exome_hg19 cap ON cap.chrom = ti.chrom AND cap.pos && ti.pos
	WHERE
		is_mate IS NOT TRUE AND -- svi contains multiple rows for each BND, these need to be removed
		(ta.genotype != './.' AND tb.genotype != './.') AND -- discard varians w/o genotype
		ta.read_depth >= 10 AND -- filter variants with insufficient read depth
		tb.read_depth >= 10 AND -- filter variants with insufficient read depth
		ta.alt_split_count + tb.alt_split_count > 0 AND -- require multiple forms of read evidence
		(ta.alt_clipped_count + tb.alt_clipped_count > 0 OR ta.alt_pe_count + tb.alt_pe_count > 0)
		--imprecise IS NOT TRUE -- discard variants with imprecise breakpoints
),
pairs_squared AS -- all possible tumor/tumor + svtype combinations
(
	SELECT *
	FROM analysis.scars_pairs sp
	CROSS JOIN svtypes
),
called_pairs_squared AS -- all tumor/tumor + svtype combinations with succesfull SV calls (some LUMPY runs failed)
(
	SELECT DISTINCT case_barcode, tumor_barcode_a, tumor_barcode_b, svt.svtype
	FROM paired_geno_sv
	CROSS JOIN svtypes svt
),
paired_geno_sv_counts_by_case_type_altcount AS
(
	SELECT cps.case_barcode, cps.tumor_barcode_a, cps.tumor_barcode_b, cps.svtype,
		COUNT(geno_a)::integer AS n,
		COUNT(CASE WHEN geno_a != '0/0' AND geno_b != '0/0' THEN 1 END)::integer AS s,
		COUNT(CASE WHEN geno_a != '0/0' AND geno_b = '0/0' THEN 1 END)::integer AS i,
		COUNT(CASE WHEN geno_a = '0/0' AND geno_b != '0/0' THEN 1 END)::integer AS r,
		COUNT(CASE WHEN geno_a = '0/0' AND geno_b = '0/0' THEN 1 END)::integer AS e
	FROM paired_geno_sv ps
	RIGHT OUTER JOIN called_pairs_squared cps ON cps.tumor_barcode_a = ps.tumor_barcode_a AND cps.tumor_barcode_b = ps.tumor_barcode_b AND cps.svtype = ps.svtype
	GROUP BY 1,2,3,4
),
coverage_pairs AS
(
	SELECT case_barcode, tumor_barcode_a, tumor_barcode_b, ca.cumulative_coverage::numeric / 10e6 AS cov_a, cb.cumulative_coverage::numeric / 10e6 AS cov_b
	FROM analysis.scars_pairs sp
	INNER JOIN analysis.coverage ca ON ca.aliquot_barcode = sp.tumor_barcode_a
	INNER JOIN analysis.coverage cb ON cb.aliquot_barcode = sp.tumor_barcode_b
	WHERE ca.coverage = 10 AND cb.coverage = 10
),
coverage_adj_rates AS
(
	SELECT
		ps.case_barcode, ps.tumor_barcode_a, ps.tumor_barcode_b, ps.svtype,
		n, s, r, i, e,
		ROUND(s/((cov_a+cov_b)/2),4) AS mf_s,
		ROUND(r/cov_b,4) AS mf_r,
		ROUND(i/cov_a,4) AS mf_i,
		ROUND((i+s)/cov_a,4) AS mf_is,
		ROUND((r+s)/cov_b,4) AS mf_rs
	FROM pairs_squared ps
	LEFT JOIN paired_geno_sv_counts_by_case_type_altcount dat ON dat.tumor_barcode_a = ps.tumor_barcode_a AND dat.tumor_barcode_b = ps.tumor_barcode_b AND dat.svtype = ps.svtype
	INNER JOIN coverage_pairs cp ON cp.tumor_barcode_a = ps.tumor_barcode_a AND cp.tumor_barcode_b = ps.tumor_barcode_b 
)
SELECT *
FROM coverage_adj_rates
ORDER BY 1,4