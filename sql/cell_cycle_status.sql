/*
---
Date: April 28, 2020
Author: FP Barthel
---
Annotate cell cycle alteration status for the scars set
- List of alterations included:
	* RB1 mutation, RB1 deletion
	* TP53 mutation, TP53 deletion
	* CDKN2A deletion
	* CDK4 amplification
	* CDK6 amplification
	* CCND2 amplification
	* MDM2 amplification
	* MDM4 amplification
--
*/
--SELECT * FROM analysis.scars_pairs

WITH
muts (gene_symbol) AS 
(
	VALUES ('TP53'), ('RB1')
),
cnvs (gene_symbol, direction) AS
(
	VALUES ('TP53', -2), ('RB1', -2), ('CDKN2A', -2), ('CDK4', 2), ('CDK6', 2), ('CCND2', 2), ('MDM2', 2), ('MDM4', 2)
),
called_muts AS
(
	SELECT
		tumor_pair_barcode,
		pg.case_barcode,
		pg.tumor_barcode_a,
		pg.tumor_barcode_b,
		pg.gene_symbol,
		variant_id,
		variant_type,
		pg.variant_classification,
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
	INNER JOIN muts ON muts.gene_symbol = pg.gene_symbol
	INNER JOIN variants.variant_classifications vs ON vs.variant_classification = pg.variant_classification
	WHERE
		(af_a > 0.10 OR af_b > 0.10) AND -- filter variants to limit those where AF > 0.10 in at least one sample in tumor pair
		vs.variant_classification_priority > 0
),
cnvs_squared AS
(
	SELECT
		sp.case_barcode,
		sp.tumor_barcode_a,
		sp.tumor_barcode_b,
		cnvs.gene_symbol || (CASE WHEN cnvs.direction = -2 THEN ' del' WHEN cnvs.direction = 2 THEN ' amp' ELSE ' NA' END) AS variant,
		CASE
			WHEN cnv_b.hlvl_call = cnvs.direction  AND cnv_a.hlvl_call = cnvs.direction THEN 'S'
			WHEN cnv_b.hlvl_call != cnvs.direction AND cnv_a.hlvl_call = cnvs.direction THEN 'I'
			WHEN cnv_b.hlvl_call = cnvs.direction AND cnv_a.hlvl_call != cnvs.direction THEN 'R'
			WHEN cnv_b.hlvl_call != cnvs.direction AND cnv_a.hlvl_call != cnvs.direction THEN 'wt'
		ELSE NULL END status
	FROM analysis.scars_pairs sp 
	CROSS JOIN cnvs
	INNER JOIN analysis.gatk_cnv_by_gene cnv_a ON sp.tumor_barcode_a = cnv_a.aliquot_barcode AND cnv_a.gene_symbol = cnvs.gene_symbol
	INNER JOIN analysis.gatk_cnv_by_gene cnv_b ON sp.tumor_barcode_b = cnv_b.aliquot_barcode AND cnv_b.gene_symbol = cnvs.gene_symbol
),
muts_squared AS
(
	SELECT 
		sp.case_barcode,
		sp.tumor_barcode_a,
		sp.tumor_barcode_b,
		muts.gene_symbol || ' mut' AS variant,
		CASE
			WHEN bool_or(alt_count_a > 0 AND alt_count_b > 0) THEN 'S'
			WHEN bool_or(alt_count_a > 0 AND alt_count_b = 0) THEN 'I'
			WHEN bool_or(alt_count_a = 0 AND alt_count_b > 0) THEN 'R'
			WHEN bool_or(alt_count_a = 0 AND alt_count_b = 0) THEN 'wt'
		ELSE 'wt' END status
	FROM analysis.scars_pairs sp
	CROSS JOIN muts
	LEFT JOIN called_muts cm ON sp.tumor_barcode_a = cm.tumor_barcode_a AND sp.tumor_barcode_b = cm.tumor_barcode_b
	GROUP BY 1,2,3,4
),
muts_cnvs_merge AS
(
	SELECT * FROM muts_squared
	UNION
	SELECT * FROM cnvs_squared
)
SELECT * FROM muts_cnvs_merge