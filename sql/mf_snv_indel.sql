/*
---
Date: April 20, 2020
Author: FP Barthel
Update: December 22, 2020 (By Emre Kocakavuk)
---
Determine the mutation burden by variant type and fraction for each case in the RT-scars cohort
- Limited to SNVs and indels called by M2 in multi-sample and single-sample mode
- Coverage adjustment is used and a threshold of 10x coverage is used
- Both queries that use alt_counts > 0 and ssm2_pass_call calling thresholds are included
- AF > 0.10 in at least one sample for each tumor pair is required
- Integrate repeatmasker info
---
*/
WITH
vtypes (variant_type) AS
(
	VALUES ('DEL'),('INS'),('SNP'),('DNP'),('TNP')
),
pairs_squared AS -- all possible tumor/tumor + variant type combinations
(
	SELECT *
	FROM analysis.scars_pairs sp
	CROSS JOIN vtypes
),
paired_geno AS -- filter and reformat paired genotypes (SNVs + Indels)
(
	SELECT
		tumor_pair_barcode,
		pg.case_barcode,
		pg.tumor_barcode_a,
		pg.tumor_barcode_b,
		pg.variant_id,
		pg.variant_type,
		rm.reptype,
		ref_count_a,
		alt_count_a,
		af_a,
		mutect2_call_a AS ssm2_pass_call_a,
		ref_count_b,
		alt_count_b,
		af_b,
		mutect2_call_b AS ssm2_pass_call_b
	FROM variants.pgeno pg -- 
	INNER JOIN analysis.scars_pairs sp ON sp.tumor_barcode_a = pg.tumor_barcode_a AND sp.tumor_barcode_b = pg.tumor_barcode_b
	LEFT JOIN variants.passanno pa ON pg.variant_id = pa.variant_id
	LEFT JOIN ref.repeatmasker rm ON pg.chrom = rm.chrom AND pg.pos && rm.pos
	WHERE
		(alt_count_a + ref_count_a >= 10) AND
		(alt_count_b + ref_count_b >= 10) AND -- coverage threshold because cannot call mutations on low coverage sites
		(af_a > 0.10 OR af_b > 0.10) -- filter variants to limit those where AF > 0.10 in at least one sample in tumor pair
),
called_pairs_squared AS -- all tumor/tumor + svtype combinations with succesfull SV calls (some LUMPY runs failed)
(
	SELECT DISTINCT case_barcode, tumor_barcode_a, tumor_barcode_b, svt.variant_type
	FROM paired_geno
	CROSS JOIN vtypes svt
	--CROSS JOIN ref.repeatmasker_reptypes rt
),
paired_geno_counts_by_case_type_altcount AS
(
	SELECT cps.case_barcode, cps.tumor_barcode_a, cps.tumor_barcode_b, cps.variant_type, ps.reptype,
		COUNT(alt_count_a)::integer AS n,
		COUNT(CASE WHEN alt_count_a > 0 AND alt_count_b > 0 THEN 1 END)::integer AS s,
		COUNT(CASE WHEN alt_count_a = 0 AND alt_count_b > 0 THEN 1 END)::integer AS r,
		COUNT(CASE WHEN alt_count_a > 0 AND alt_count_b = 0 THEN 1 END)::integer AS i,
		COUNT(CASE WHEN alt_count_a = 0 AND alt_count_b = 0 THEN 1 END)::integer AS e
	FROM paired_geno ps
	RIGHT JOIN called_pairs_squared cps ON cps.tumor_barcode_a = ps.tumor_barcode_a AND cps.tumor_barcode_b = ps.tumor_barcode_b AND cps.variant_type = ps.variant_type
	GROUP BY 1,2,3,4,5
),
paired_geno_counts_by_case_type_ssm2 AS
(
	SELECT case_barcode, tumor_barcode_a, tumor_barcode_b, variant_type, reptype,
		COUNT(*)::integer AS n,
		COUNT(CASE WHEN ssm2_pass_call_a IS TRUE AND ssm2_pass_call_b IS TRUE THEN 1 END)::integer AS s,
		COUNT(CASE WHEN ssm2_pass_call_a IS NOT TRUE AND ssm2_pass_call_b IS TRUE THEN 1 END)::integer AS r,
		COUNT(CASE WHEN ssm2_pass_call_a IS TRUE AND ssm2_pass_call_b IS NOT TRUE THEN 1 END)::integer AS i,
		COUNT(CASE WHEN ssm2_pass_call_a IS NOT TRUE AND ssm2_pass_call_b IS NOT TRUE THEN 1 END)::integer AS e
	FROM paired_geno
	GROUP BY 1,2,3,4,5
),
coverage_pairs AS
(
	SELECT case_barcode, tumor_barcode_a, tumor_barcode_b, ca.cumulative_coverage::numeric / 10e6 AS cov_a, cb.cumulative_coverage::numeric / 10e6 AS cov_b
	FROM analysis.scars_pairs sp
	INNER JOIN analysis.coverage ca ON ca.aliquot_barcode = sp.tumor_barcode_a
	INNER JOIN analysis.coverage cb ON cb.aliquot_barcode = sp.tumor_barcode_b
	WHERE ca.coverage = 10 AND cb.coverage = 10
),
coverage_adj_rates AS -- substitute LEFT JOIN for `ssm2` table above if wanting to use that as calling criterium
(
	SELECT
		ps.case_barcode, ps.tumor_barcode_a, ps.tumor_barcode_b, ps.variant_type, reptype,
		n, s, r, i, e,
		ROUND(s/((cov_a+cov_b)/2),4) AS mf_s,
		ROUND(r/cov_b,4) AS mf_r,
		ROUND(i/cov_a,4) AS mf_i,
		ROUND((i+s)/cov_a,4) AS mf_is,
		ROUND((r+s)/cov_b,4) AS mf_rs
	FROM pairs_squared ps
	LEFT JOIN paired_geno_counts_by_case_type_altcount t1 ON t1.tumor_barcode_a = ps.tumor_barcode_a AND t1.tumor_barcode_b = ps.tumor_barcode_b AND t1.variant_type = ps.variant_type
	INNER JOIN coverage_pairs cp ON cp.tumor_barcode_a = t1.tumor_barcode_a AND cp.tumor_barcode_b = t1.tumor_barcode_b 
)
SELECT *
FROM coverage_adj_rates
ORDER BY 1,4