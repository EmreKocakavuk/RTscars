SELECT
        pa.variant_id,
        pg.case_barcode,
        pg.tumor_barcode_a,
        pg.tumor_barcode_b,
        pa.variant_type,
		rm.reptype,
        pa.chrom,
        lower(pa.pos) AS pos,
        pa.ref,
        pa.alt,
        feature_type, 
        dist,
        (CASE
         WHEN pg.alt_count_a > 0 AND pg.alt_count_b > 0 THEN 'S'
         WHEN pg.alt_count_a > 0 AND NOT pg.alt_count_b > 0 THEN 'P'
         WHEN pg.alt_count_b > 0 AND NOT pg.alt_count_a > 0 THEN 'R'
         ELSE NULL::text
         END) AS fraction
    FROM variants.pgeno pg
    INNER JOIN variants.passanno pa ON pa.variant_id = pg.variant_id
     INNER JOIN analysis.variants_nonbdb vn ON vn.variant_id = pg.variant_id
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