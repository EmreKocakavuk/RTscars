SELECT
        sampleid,
        variant_type,
        chrom,
        lower(pos) AS start_pos,
        upper(pos)-1 AS end_pos,
        ref,
        alt
FROM hwg
WHERE alt_count > 0 AND ref_count >= 10