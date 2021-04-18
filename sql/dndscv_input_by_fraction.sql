SELECT
    pg.case_barcode,
    pa.variant_type,
	(CASE WHEN pg.chrom = 23 THEN 'X' ELSE pg.chrom::varchar(2) END) AS chrom,
 	 lower(pg.pos) AS pos,
 	 pg.ref,
 	 pg.alt AS mut,
    (CASE
     WHEN pg.alt_count_a > 0 AND pg.alt_count_b > 0 THEN 'S'
     WHEN pg.alt_count_a > 0 AND NOT pg.alt_count_b > 0 THEN 'P'
     WHEN pg.alt_count_b > 0 AND NOT pg.alt_count_a > 0 THEN 'R'
     ELSE NULL::text
     END) AS fraction,
	 st.idh_codel_subtype AS subtype, 
	 pg.gene_symbol
FROM variants.pgeno pg
INNER JOIN analysis.scars_pairs sp ON sp.tumor_barcode_a = pg.tumor_barcode_a AND sp.tumor_barcode_b = pg.tumor_barcode_b
LEFT JOIN clinical.subtypes st ON st.case_barcode = pg.case_barcode
INNER JOIN variants.passanno pa ON pa.variant_id = pg.variant_id
WHERE
  (pg.alt_count_a > 0 OR pg.alt_count_b > 0) AND
  (pg.ref_count_a + pg.alt_count_a >= 10) AND
  (pg.ref_count_b + pg.alt_count_b >= 10) AND
  (af_a > 0.10 OR af_b > 0.10)