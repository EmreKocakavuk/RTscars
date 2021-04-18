/*
---
Created: April 30, 2020
Updated: May 7, 2020
Author: FP Barthel
---
** MODIFIED FROM `taylor_aneuploidy.sql` from GLASS for HWMF copy number data analysis
---
This specific script takes arm level aneuploidy and computes the Taylor aneuploidy score
---
*/
WITH
taylor_aneuploidy AS
(
	SELECT
		sampleid,
		sum(CASE WHEN arm_call <> 0 THEN 1 ELSE 0 END) OVER w::integer AS aneuploidy_score,
		sum(CASE WHEN arm_call = 1 THEN 1 ELSE 0 END) OVER w::integer AS aneuploidy_amp_score,
		sum(CASE WHEN arm_call = -1 THEN 1 ELSE 0 END) OVER w::integer AS aneuploidy_del_score,
		first_value(arm_num_seg) OVER w2::integer AS max_loss_arm_n,
		first_value(arm_ploidy_wmean) OVER w2 AS max_loss_arm_wploidy,
		first_value(arm_ploidy_wsd) OVER w2 AS max_loss_arm_wsd,
		first_value(arm_num_seg) OVER w3::integer AS max_gain_arm_n,
		first_value(arm_ploidy_wmean) OVER w3 AS max_gain_arm_wploidy,
		first_value(arm_ploidy_wsd) OVER w3 AS max_gain_arm_wsd,
		row_number() OVER w AS rn
	FROM hwa
	WINDOW
		w  AS (PARTITION BY sampleid),
		w2 AS (PARTITION BY sampleid ORDER BY arm_ploidy_wmean ASC NULLS LAST),
		w3 AS (PARTITION BY sampleid ORDER BY arm_ploidy_wmean DESC NULLS LAST)
),
ploidy AS
(
	SELECT
		sampleid,
		COUNT(*) AS num_seg,
		sum((upper(pos)::numeric - lower(pos)) * copy_number) / sum(upper(pos)::numeric - lower(pos)) AS ploidy,
		sum(upper(pos)::numeric - lower(pos)) AS seg_size
	FROM hwc
	WHERE chrom < 23
	GROUP BY 1
)
SELECT
	ta.sampleid,
	ploidy AS sample_ploidy,
	aneuploidy_score,
	aneuploidy_amp_score,
	aneuploidy_del_score,
	max_loss_arm_n,
	max_loss_arm_wploidy,
	max_loss_arm_wsd,
	max_gain_arm_n,
	max_gain_arm_wploidy,
	max_gain_arm_wsd
FROM taylor_aneuploidy ta
LEFT JOIN ploidy ON ploidy.sampleid = ta.sampleid
WHERE rn = 1