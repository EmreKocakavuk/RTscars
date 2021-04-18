/*
---
Created: April 30, 2020
Updated: May 7, 2020
Author: FP Barthel
---
** MODIFIED FROM `taylor_aneuploidy.sql` from GLASS for HWMF copy number data analysis
+ HWMF copy number data differs because it is not a conventional segmentation file with log2 copy ratios
+ Instead, it provides ploidy per segment
---
This particular script computes arm level copy number, see the original script by FP Barthel in the GLASS project for documentation
---
*/
WITH
hwc_calls AS
(
	SELECT
		sampleid,
		chrom,
		pos,
		copy_number AS ploidy,
		minor_allele_ploidy,
		major_allele_ploidy,
		(CASE 
		 WHEN minor_allele_ploidy = major_allele_ploidy AND minor_allele_ploidy != 0 THEN 0 -- WGD should be considered as "euploid"
		 WHEN minor_allele_ploidy = 0 OR major_allele_ploidy = 0 THEN -1 
		 WHEN major_allele_ploidy > 1 THEN 1
		 ELSE NULL END) AS cnv -- a specific modification made for the HW dataset is to manually set thresholds of copy number
	FROM hwc
),
arrange_segments_by_arm AS
(
    SELECT
		sampleid,
		gs.chrom,
		ca.arm,
		ca.pos * gs.pos AS pos,
		ploidy,
		cnv,
		(upper(ca.pos * gs.pos)::decimal - lower(ca.pos * gs.pos)) as seg_size,
		(sum((upper(ca.pos * gs.pos)::decimal - lower(ca.pos * gs.pos))) OVER w) as arm_size,
		row_number() OVER w2 - row_number() OVER w3 AS grp
    FROM hwc_calls gs
	INNER JOIN ref.chr_arms ca ON ca.chrom = gs.chrom AND ca.pos && gs.pos
	WHERE 
		gs.chrom NOT IN (23,24) AND 
		ca.arm <> '21p'
	WINDOW 
		w AS (PARTITION BY gs.sampleid, ca.arm),
		w2 AS (PARTITION BY gs.sampleid, ca.arm ORDER BY ca.pos * gs.pos),
		w3 AS (PARTITION BY gs.sampleid, ca.arm, cnv ORDER BY ca.pos * gs.pos)
),
unfiltered_join_adjacent_segments AS
(
	SELECT
		sampleid,
		chrom,
		arm,
		grp,
		int4range(min(lower(pos)),max(upper(pos))) AS pos,
		cnv,
		COUNT(*) AS num_seg,
		sum(seg_size * ploidy) / sum(seg_size) AS wploidy,
		(CASE 
			WHEN COUNT(*) > 1 AND sum(t1.seg_size * (t1.ploidy ^ 2::numeric)) > ((sum(t1.seg_size * t1.ploidy) ^ 2::numeric) / sum(t1.seg_size)) THEN sqrt( (sum(seg_size * ploidy^2) - (sum(seg_size * ploidy)^2)/sum(seg_size)) / (sum(seg_size)-1) ) 
			ELSE 0 END) AS wsd,
		sum(seg_size) AS seg_size,
		min(arm_size) AS arm_size
	FROM arrange_segments_by_arm t1
	--WHERE ploidy >= 0 AND seg_size > 0
	GROUP BY 1,2,3,4,6
	ORDER BY 1,2,3,5
),
filtered_join_adjacent_segments AS
(
	SELECT
		t1.sampleid,
		t1.chrom,
		t1.arm,
		t1.grp,
		int4range(min(lower(t1.pos)),max(upper(t1.pos))) AS pos,
		t1.cnv,
		COUNT(*) AS fnum_seg,
		sum(t1.seg_size * ploidy) / sum(t1.seg_size) AS fwploidy,
		(CASE 
			WHEN COUNT(*) > 1 AND sum(t1.seg_size * (t1.ploidy ^ 2::numeric)) > ((sum(t1.seg_size * t1.ploidy) ^ 2::numeric) / sum(t1.seg_size)) THEN sqrt( (sum(t1.seg_size * ploidy^2) - (sum(t1.seg_size * ploidy)^2)/sum(t1.seg_size)) / (sum(t1.seg_size)-1) ) 
			ELSE 0 END) AS fwsd,
		sum(t1.seg_size) AS seg_size,
		min(t1.arm_size) AS arm_size
	FROM arrange_segments_by_arm t1
	INNER JOIN unfiltered_join_adjacent_segments us ON us.sampleid = t1.sampleid AND us.chrom = t1.chrom AND us.arm = t1.arm AND us.grp = t1.grp AND us.cnv = t1.cnv
	WHERE
		--ploidy >= 0 AND
		--t1.seg_size > 0 AND
		wsd = 0 OR ((ploidy - wploidy) > -2.0 * wsd AND (ploidy - wploidy) < 2.0 * wsd)
	GROUP BY 1,2,3,4,6
	ORDER BY 1,2,3,5
),
sort_by_longest_segment AS
(
	SELECT
		t2.sampleid,
		t2.chrom,
		t2.arm,
		t2.grp,
		t2.cnv,
		row_number() OVER w AS partion_id,
		t2.seg_size/t2.arm_size AS prop_arm,
		num_seg,
		wploidy,
		wsd,
		fss.fnum_seg,
		fss.fwploidy,
		fss.fwsd
	FROM unfiltered_join_adjacent_segments t2
	LEFT JOIN filtered_join_adjacent_segments fss ON fss.sampleid = t2.sampleid AND fss.chrom = t2.chrom AND fss.arm = t2.arm AND fss.grp = t2.grp AND fss.cnv = t2.cnv
	WINDOW w AS (PARTITION BY t2.sampleid, t2.chrom, t2.arm ORDER BY t2.seg_size/t2.arm_size DESC)
),
cnv_by_arm AS
(
	SELECT
		sampleid,
		chrom,
		arm,
		(CASE
		 WHEN prop_arm > 0.80 AND cnv = 1 THEN 1
		 WHEN prop_arm > 0.80 AND cnv = -1 THEN -1
		 WHEN prop_arm > 0.80 AND cnv = 0 THEN 0
		 ELSE NULL END) AS arm_call,
		fnum_seg AS arm_num_seg,
		fwploidy AS arm_ploidy_wmean,
		fwsd AS arm_ploidy_wsd
	FROM sort_by_longest_segment t3
	WHERE partion_id = 1
)
SELECT * FROM cnv_by_arm