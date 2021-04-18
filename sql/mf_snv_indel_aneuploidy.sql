/*
---
Date: April 21, 2020
Author: FP Barthel
---
Determine the number of SNPs/indels per sample and arm
- Filter variants within the generic GLASS capture region
- Normalize using the number of base pairs on each arm that lie within the captured region
---
*/
WITH
fractions (vfraction) AS
(
	VALUES ('S'),('I'),('R')
),
vtypes (variant_type) AS
(
	VALUES ('DEL'),('INS'),('SNP'),('DNP'),('TNP')
),
selected_arms AS 
(
  SELECT
		arm,
		chrom,
		pos
	FROM ref.chr_arms ca
),
coverage_by_arm AS -- calculate number of base pairs for each arm that lies within common captured region
(
	SELECT
		arm,
		sum(upper(ca.pos * cap.pos) - lower(ca.pos * cap.pos)) AS covered_arm_size
	FROM ref.chr_arms ca
	INNER JOIN ref.intersect_exome_hg19 cap ON cap.chrom = ca.chrom AND cap.pos && ca.pos
	GROUP BY 1
),
paired_cnv AS
(
	SELECT
		case_barcode,
		tumor_barcode_a,
		tumor_barcode_b,
		c1.arm,
		ca.covered_arm_size,
		c1.arm_call arm_call_a,
		c2.arm_call arm_call_b,
		CASE
			WHEN c2.arm_call = -1 AND c1.arm_call = 1 THEN 'Unstable'
			WHEN c2.arm_call = 1 AND c1.arm_call = -1 THEN 'Unstable'
			WHEN c2.arm_call = 0 AND c1.arm_call != 0 THEN 'Acquired euploid'
			WHEN c2.arm_call = 1 AND c1.arm_call != 1 THEN 'Acquired gain'
			WHEN c2.arm_call = -1 AND c1.arm_call != -1 THEN 'Acquired loss'
			WHEN c2.arm_call = 0 OR c1.arm_call = 0 THEN 'Stable euploid'
			WHEN c2.arm_call = 1 OR c1.arm_call = 1 THEN 'Stable gain'
			WHEN c2.arm_call = -1 OR c1.arm_call = -1 THEN 'Stable loss'
			ELSE 'Stable euploid'
		END arm_change_relaxed,
		CASE
			WHEN c2.arm_call IS NULL OR c1.arm_call IS NULL THEN NULL
			WHEN c2.arm_call = -1 AND c1.arm_call = 1 THEN 'Unstable'
			WHEN c2.arm_call = 1 AND c1.arm_call = -1 THEN 'Unstable'
			WHEN c2.arm_call = 0 AND c1.arm_call != 0 THEN 'Acquired euploid'
			WHEN c2.arm_call = 1 AND c1.arm_call != 1 THEN 'Acquired gain'
			WHEN c2.arm_call = -1 AND c1.arm_call != -1 THEN 'Acquired loss'
			WHEN c2.arm_call = 0 OR c1.arm_call = 0 THEN 'Stable euploid'
			WHEN c2.arm_call = 1 OR c1.arm_call = 1 THEN 'Stable gain'
			WHEN c2.arm_call = -1 OR c1.arm_call = -1 THEN 'Stable loss'
			ELSE 'Stable euploid'
		END arm_change_strict
	FROM analysis.scars_pairs sp
	INNER JOIN analysis.gatk_cnv_by_arm c1 ON c1.aliquot_barcode = sp.tumor_barcode_a
	INNER JOIN analysis.gatk_cnv_by_arm c2 ON c2.aliquot_barcode = sp.tumor_barcode_b AND c1.arm = c2.arm
	INNER JOIN coverage_by_arm ca ON ca.arm = c2.arm
),
paired_cnv_vfraction_vtype AS
(
	SELECT *
	FROM paired_cnv
	CROSS JOIN fractions, vtypes
),
selected_variants AS
(	
	SELECT
		tumor_pair_barcode,
		pg.case_barcode,
		pg.tumor_barcode_a,
		pg.tumor_barcode_b,
		arm,
		variant_id,
		variant_type,
		ref_count_a,
		alt_count_a,
		af_a,
		mutect2_call_a AS ssm2_pass_call_a,
		ref_count_b,
		alt_count_b,
		af_b,
		mutect2_call_b AS ssm2_pass_call_b,
		CASE
			WHEN alt_count_a > 0 AND alt_count_b > 0 THEN 'S'
			WHEN alt_count_a > 0 AND alt_count_b = 0 THEN 'I'
			WHEN alt_count_a = 0 AND alt_count_b > 0 THEN 'R'
		END vfraction
	FROM variants.pgeno pg -- 
	INNER JOIN analysis.scars_pairs sp ON sp.tumor_barcode_a = pg.tumor_barcode_a AND sp.tumor_barcode_b = pg.tumor_barcode_b
	INNER JOIN ref.intersect_exome_hg19 cap ON cap.chrom = pg.chrom AND cap.pos && pg.pos -- filter those variants within captured region
	INNER JOIN selected_arms sa ON pg.chrom = sa.chrom AND pg.pos && sa.pos -- annotate arm
	WHERE
		(alt_count_a + ref_count_a >= 10) AND
		(alt_count_b + ref_count_b >= 10) AND -- coverage threshold because cannot call mutations on low coverage sites
		(af_a > 0.10 OR af_b > 0.10) -- filter variants to limit those where AF > 0.10 in at least one sample in tumor pair
),
arm_counts AS
(
	SELECT
		pc.tumor_barcode_a,
		pc.tumor_barcode_b,
		pc.case_barcode,
		pc.arm,
		pc.arm_call_a,
		pc.arm_call_b,
		covered_arm_size,
		arm_change_strict,
		arm_change_relaxed,
		pc.variant_type,
		pc.vfraction,
		COUNT(sv.vfraction) AS abs_arm_count,
		ROUND((COUNT(sv.vfraction)/covered_arm_size::numeric) * 1e6,4) AS rel_arm_freq
	FROM paired_cnv_vfraction_vtype pc
	LEFT JOIN selected_variants sv ON sv.tumor_barcode_a = pc.tumor_barcode_a AND sv.tumor_barcode_b = pc.tumor_barcode_b AND sv.arm = pc.arm AND sv.vfraction = pc.vfraction AND sv.variant_type = pc.variant_type
	GROUP BY 1,2,3,4,5,6,7,8,9,10,11
)
SELECT * FROM arm_counts