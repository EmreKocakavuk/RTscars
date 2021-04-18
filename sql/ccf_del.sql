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
paired_geno AS
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
		mutect2_call_b AS ssm2_pass_call_b, 
		p1.pyclone_ccf AS ccf_a, 
		p2.pyclone_ccf AS ccf_b, 
	(CASE
	 WHEN alt_count_a > 0 AND alt_count_b > 0 			THEN 'S'
	 WHEN alt_count_a > 0 AND alt_count_b = 0 			THEN 'P'
	 WHEN alt_count_a = 0 AND alt_count_b > 0 			THEN 'R'
	 ELSE NULL END) AS fraction
	FROM variants.pgeno pg -- 
	INNER JOIN analysis.scars_pairs sp ON sp.tumor_barcode_a = pg.tumor_barcode_a AND sp.tumor_barcode_b = pg.tumor_barcode_b
	LEFT JOIN variants.passanno pa ON pg.variant_id = pa.variant_id
	LEFT JOIN ref.repeatmasker rm ON pg.chrom = rm.chrom AND pg.pos && rm.pos
	LEFT JOIN variants.passgeno p1 ON p1.variant_id = pg.variant_id AND p1.aliquot_barcode = pg.tumor_barcode_a
	LEFT JOIN variants.passgeno p2 ON p2.variant_id = pg.variant_id AND p2.aliquot_barcode = pg.tumor_barcode_b
	WHERE
		(alt_count_a + ref_count_a >= 10) AND
		(alt_count_b + ref_count_b >= 10) AND -- coverage threshold because cannot call mutations on low coverage sites
		(af_a > 0.10 OR af_b > 0.10) -- filter variants to limit those where AF > 0.10 in at least one sample in tumor pair
), 
paired_geno_ccf AS
(
SELECT *, 
(CASE
	 WHEN fraction = 'S' THEN (ccf_a+ccf_b)/2
	 WHEN fraction = 'P' THEN ccf_a
	 WHEN fraction = 'R' THEN ccf_b
	 ELSE NULL END) AS ccf
FROM paired_geno 
WHERE reptype IS NULL and variant_type = 'DEL'
)
SELECT case_barcode, variant_type, fraction, ccf, ccf_a, ccf_b
FROM paired_geno_ccf
WHERE ccf IS NOT NULL