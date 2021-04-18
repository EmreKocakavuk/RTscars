/*
---
Date: April 21, 2020
Author: FP Barthel
---
Use Taylor aneuploidy metric to quantify chromosome missegregation events
- Determine for each chromosome (autosomes only):
	+ when both arms are zero (euploid) or both arms are null (indeterminate): euploid
	+ when either arm is not euploid:
		* Simple missegregation: both arms aneuploid in same direction
		* Complex missegregation: aneuploidy varies per arm
- Aggregate counts per aliquot
	+ euploid + simple + complex = 22
	+ simple = simple_gain + simple_loss
	+ complex = complex_gain_neut + complex_loss_neut + complex_loss_gain
*/
WITH 
missegregation_counts AS
(
	SELECT
		aliquot_barcode,
		chrom,
		CASE WHEN bool_and(arm_call = 0) OR bool_and (arm_call IS NULL)THEN TRUE ELSE FALSE END AS euploid,
		CASE WHEN bool_or(arm_call != 0) THEN min(arm_call) = max(arm_call) END AS simple_missegregation,
		CASE WHEN bool_or(arm_call != 0) THEN min(arm_call) != max(arm_call) END AS complex_missegregation,
		CASE WHEN min(arm_call) = max(arm_call) THEN min(arm_call) END AS chr_call,
		COUNT(CASE WHEN arm_call = 0 THEN 1 END) n_neut,
		COUNT(CASE WHEN arm_call = -1 THEN 1 END) n_loss,
		COUNT(CASE WHEN arm_call = 1 THEN 1 END) n_gain
	FROM analysis.gatk_cnv_by_arm
	GROUP BY 1,2
	ORDER BY 1,2
),
missegregation_aggregates AS
(
	SELECT
		aliquot_barcode,
		COUNT(CASE WHEN euploid IS TRUE THEN 1 END):: integer n_euploid,
		COUNT(CASE WHEN simple_missegregation IS TRUE THEN 1 END)::integer n_simple_missegregation,
		COUNT(CASE WHEN complex_missegregation IS TRUE THEN 1 END)::integer n_complex_missegregation,
		COUNT(CASE WHEN simple_missegregation IS TRUE AND n_gain > 0 AND n_loss = 0 AND n_neut = 0 THEN 1 END)::integer n_simple_missegregation_gain,
		COUNT(CASE WHEN simple_missegregation IS TRUE AND n_gain = 0 AND n_loss > 0 AND n_neut = 0 THEN 1 END)::integer n_simple_missegregation_loss,
		COUNT(CASE WHEN complex_missegregation IS TRUE AND n_gain > 0 AND n_loss = 0 AND n_neut > 0 THEN 1 END)::integer n_complex_missegregation_gain_neut,
		COUNT(CASE WHEN complex_missegregation IS TRUE AND n_gain = 0 AND n_loss > 0 AND n_neut > 0 THEN 1 END)::integer n_complex_missegregation_loss_neut,
		COUNT(CASE WHEN complex_missegregation IS TRUE AND n_gain > 0 AND n_loss > 0 AND n_neut = 0 THEN 1 END)::integer n_complex_missegregation_gain_loss
	FROM missegregation_counts
	GROUP BY 1
	ORDER BY 1
),
missegregation_pairs AS
(
	SELECT
		case_barcode, tumor_barcode_a, tumor_barcode_b,
		ta.n_euploid AS n_euploid_a, 
		ta.n_simple_missegregation AS n_simple_missegregation_a,
		ta.n_complex_missegregation AS n_complex_missegregation_a,
		ta.n_simple_missegregation_gain AS n_simple_missegregation_gain_a,
		ta.n_simple_missegregation_loss AS n_simple_missegregation_loss_a,
		ta.n_complex_missegregation_gain_neut AS n_complex_missegregation_gain_neut_a,
		ta.n_complex_missegregation_loss_neut AS n_complex_missegregation_loss_neut_a,
		ta.n_complex_missegregation_gain_loss AS n_complex_missegregation_gain_loss_a,
		tb.n_euploid AS n_euploid_b, 
		tb.n_simple_missegregation AS n_simple_missegregation_b,
		tb.n_complex_missegregation AS n_complex_missegregation_b,
		tb.n_simple_missegregation_gain AS n_simple_missegregation_gain_b,
		tb.n_simple_missegregation_loss AS n_simple_missegregation_loss_b,
		tb.n_complex_missegregation_gain_neut AS n_complex_missegregation_gain_neut_b,
		tb.n_complex_missegregation_loss_neut AS n_complex_missegregation_loss_neut_b,
		tb.n_complex_missegregation_gain_loss AS n_complex_missegregation_gain_loss_b
	FROM analysis.scars_pairs sp
	INNER JOIN missegregation_aggregates ta ON ta.aliquot_barcode = sp.tumor_barcode_a
	INNER JOIN missegregation_aggregates tb ON tb.aliquot_barcode = sp.tumor_barcode_b
)
SELECT * FROM missegregation_pairs