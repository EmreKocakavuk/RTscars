/*
---
Date: April 30, 2020
Author: FP Barthel
---
Modified from `cnv_missegregation.sql` for Hartwig analysis
--
*/
WITH 
missegregation_counts AS
(
	SELECT
		sampleid,
		chrom,
		CASE WHEN bool_and(arm_call = 0) OR bool_and (arm_call IS NULL)THEN TRUE ELSE FALSE END AS euploid,
		CASE WHEN bool_or(arm_call != 0) THEN min(arm_call) = max(arm_call) END AS simple_missegregation,
		CASE WHEN bool_or(arm_call != 0) THEN min(arm_call) != max(arm_call) END AS complex_missegregation,
		CASE WHEN min(arm_call) = max(arm_call) THEN min(arm_call) END AS chr_call,
		COUNT(CASE WHEN arm_call = 0 THEN 1 END) n_neut,
		COUNT(CASE WHEN arm_call = -1 THEN 1 END) n_loss,
		COUNT(CASE WHEN arm_call = 1 THEN 1 END) n_gain
	FROM hwa
	GROUP BY 1,2
	ORDER BY 1,2
),
missegregation_aggregates AS
(
	SELECT
		sampleid,
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
)
SELECT * FROM missegregation_aggregates