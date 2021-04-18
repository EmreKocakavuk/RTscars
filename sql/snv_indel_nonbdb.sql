/*
---
Date: April 22, 2020
Author: FP Barthel
---
Map variants to non-B DNA feature distances
---
*/

SELECT
    pa.variant_id,
    pg.case_barcode,
    pg.tumor_barcode_a,
    pg.tumor_barcode_b,
    pa.variant_type,
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
INNER JOIN analysis.scars_pairs sp ON sp.tumor_barcode_a = pg.tumor_barcode_a AND sp.tumor_barcode_b = pg.tumor_barcode_b
INNER JOIN variants.passanno pa ON pa.variant_id = pg.variant_id
INNER JOIN analysis.variants_nonbdb vn ON vn.variant_id = pg.variant_id
WHERE
  (pg.alt_count_a > 0 OR pg.alt_count_b > 0) AND
  (pg.ref_count_a + pg.alt_count_a >= 10) AND
  (pg.ref_count_b + pg.alt_count_b >= 10) AND
  (af_a > 0.10 OR af_b > 0.10)