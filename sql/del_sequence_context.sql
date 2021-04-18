/*
---
Date: April 22, 2020
Author: FP Barthel, Javad Noorbakhsh (helped with repeat and homology function)
---
Details the sequence context of small deletions in the RT scars cohort
- The function `infer_homology` takes the deleted sequence and flanking 
  sequences and determines whether there is evidence for microhomology.
  It outputs the homologous bases, if any.
- The function `walk_repeats` takes a single sequence and analyses it for
  repeated sequences. It returns four outputs:
  	* repeat_sequence: the repeat unit (the sequence that is repeated)
  	* num_repeats: the number of repeats
  	* gap: number of nucleotides that partially match the repeat unit
  	* mismatches: number of nucleotides that do not match repeat sequence (max 1 allowed)
  Examples:
  	+ TATAG is a TAx2 repeat with a gap and mismatch
	+ TATAT is TA repeat with a gap and no mismatches
	+ TAGTAGTG is a TAG repeat (x2) with a gap and a single mismatch
	+ TATACATATAG is not a repeat (there is a mismatch in the gap region, G, and a C  mismatch, so 2 mismatches --> no repeat)
---
*/
WITH 
del_context AS
(
	SELECT
		rm.reptype, 
		pg.case_barcode,
		pg.tumor_barcode_a,
		pg.tumor_barcode_b,
		pa.variant_type,
		pa.chrom,
		pa.gene_symbol,
		lower(pa.pos) AS pos,
		pa.ref,
		pa.alt,
		af_a,
		af_b,
		CASE
			WHEN pg.alt_count_a > 0 AND pg.alt_count_b > 0 THEN 'S'
			WHEN pg.alt_count_a > 0 AND NOT pg.alt_count_b > 0 THEN 'P'
			WHEN pg.alt_count_b > 0 AND NOT pg.alt_count_a > 0 THEN 'R'
			ELSE NULL::text
		END AS fraction,
		cov1.cumulative_coverage AS cumulative_coverage_a,
		cov2.cumulative_coverage AS cumulative_coverage_b,
		length(pa.ref) AS length_ref,
		length(pa.alt) AS length_alt,    
		length(reference_context) AS length_ref_contxt,
		reference_context,
		--substring(reference_context,10-(length(pa.ref)-1),length(pa.ref)) AS ref_before_del,
		--substring(reference_context,10+1,length(pa.ref)-1) AS actual_deletion
		substring(reference_context,1,10) AS ref_before_del,
		substring(reference_context,10+length(pa.ref),10) AS ref_after_del,
		substring(reference_context,10+1,length(pa.ref)-1) AS actual_deletion
	FROM variants.pgeno pg
	INNER JOIN variants.passanno pa ON pa.variant_id = pg.variant_id
	INNER JOIN analysis.scars_pairs sp ON sp.tumor_barcode_a = pg.tumor_barcode_a AND sp.tumor_barcode_b = pg.tumor_barcode_b
	INNER JOIN analysis.coverage cov1 ON cov1.aliquot_barcode = pg.tumor_barcode_a
	INNER JOIN analysis.coverage cov2 ON cov2.aliquot_barcode = pg.tumor_barcode_b
	LEFT JOIN ref.repeatmasker rm ON pg.chrom = rm.chrom AND pg.pos && rm.pos
	WHERE
		(pg.alt_count_a > 0 OR pg.alt_count_b > 0) AND
		pa.variant_type = 'DEL' AND
		cov1.coverage = 10 AND
		cov2.coverage = 10 AND
		(pg.ref_count_a + pg.alt_count_a >= 10) AND
		(pg.ref_count_b + pg.alt_count_b >= 10)
),
del_context_homology AS
(
    SELECT 
        *,
        infer_homology(ref_before_del, actual_deletion, ref_after_del) AS homology,
        length(infer_homology(ref_before_del, actual_deletion, ref_after_del)) AS homology_length
    FROM del_context
)
SELECT
    t2.*,

    ROUND(length(regexp_replace(ref_before_del, '[AT]','', 'g')) / length(ref_before_del)::numeric,4) AS gc_before,
	ROUND(length(regexp_replace(actual_deletion, '[AT]','', 'g')) / length(actual_deletion)::numeric,4) AS gc_del,
	ROUND(length(regexp_replace(ref_after_del, '[AT]','', 'g')) / length(ref_after_del)::numeric,4) AS gc_after,
    
    t3.repeat_sequence AS mh_repeat_sequence,
    t3.num_repeats AS mh_num_repeats,
    t3.gap AS mh_gap,
    t3.mismatches AS mh_mismatches,
    
    t4.repeat_sequence AS del_repeat_sequence,
    t4.num_repeats AS del_num_repeats,
    t4.gap AS del_gap,
    t4.mismatches AS del_mismatches,
    
    t5.repeat_sequence AS before_repeat_sequence,
    t5.num_repeats AS before_num_repeats,
    t5.gap AS before_gap,
    t5.mismatches AS before_mismatches,
    
    t6.repeat_sequence AS after_repeat_sequence,
    t6.num_repeats AS after_num_repeats,
    t6.gap AS after_gap,
    t6.mismatches AS after_mismatches
FROM 
    del_context_homology AS t2,
    walk_repeats(homology) AS t3,
    walk_repeats(actual_deletion) AS t4,
    walk_repeats(ref_before_del) AS t5,
    walk_repeats(ref_after_del) AS t6