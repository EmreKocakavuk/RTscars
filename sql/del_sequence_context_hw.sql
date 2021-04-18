/*
---
Created: May 4, 2020
Author: FP Barthel
---
Get sequence context analysis for HW data
---
*/
WITH
del_context AS
(
	SELECT
		hwg.sampleid,
		eventid,
		chrom,
		pos,
		ref,
		alt,
		variant_type,
		genotype,
		depth,
		ref_count,
		alt_count,
		af,
		seq,
		substring(seq,1,20) AS seq_before,
		substring(seq,21,length(ref)-1) AS seq_del,
		substring(seq,20+length(ref),20) AS seq_after
	FROM hwg
	INNER JOIN hmf_meta hm ON hm.sampleid = hwg.sampleid
	WHERE 
		variant_type = 'DEL' AND 
		length(ref) < 22 AND -- maximum deletion size = 20nt, the `ref` includes a single (the first) nt that is not actually deleted
		genotype != '0/0'
)
SELECT
	dc.*,
	ROUND(length(regexp_replace(seq_before, '[AT]','', 'g')) / 20::numeric,4) AS gc_before,
	ROUND(length(regexp_replace(seq_del, '[AT]','', 'g')) / length(seq_del)::numeric,4) AS gc_del,
	ROUND(length(regexp_replace(seq_after, '[AT]','', 'g')) / 20::numeric,4) AS gc_after,
	
	infer_homology(seq_before, seq_del, seq_after) AS homology,
	length(infer_homology(seq_before, seq_del, seq_after)) AS homology_length,
	
	t3.repeat_sequence AS mh_repeat_sequence,
    t3.num_repeats AS mh_num_repeats,
    t3.gap AS mh_gap,
    t3.mismatches AS mh_mismatches
FROM
	del_context dc,
	walk_repeats(infer_homology(seq_before, seq_del, seq_after)) AS t3