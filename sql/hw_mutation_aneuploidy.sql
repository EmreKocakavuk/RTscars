/*
---
Created: April 21, 2020
Updated: May 7, 2020
Author: FP Barthel
---
Determine the number of SNPs/indels per sample and arm
- Filter variants within the generic GLASS capture region
- Normalize using the number of base pairs on each arm that lie within the captured region
---
*/
WITH
vtypes (variant_type) AS
(
	VALUES ('DEL'),('INS'),('SNP'),('DNP'),('TNP')
),
coverage_by_arm AS -- calculate number of base pairs for each arm that lies within common captured region
(
	SELECT
		chrom,
		arm,
		pos,
		upper(ca.pos) - lower(ca.pos) AS covered_arm_size
	FROM ref.chr_arms ca
),
sample_arms AS
(
	SELECT
		sampleid,
		hwa.arm,
		arm_call,
		arm_ploidy_wmean AS arm_ploidy,
		covered_arm_size,
		variant_type
	FROM hwa
	INNER JOIN coverage_by_arm ca ON hwa.arm = ca.arm
	CROSS JOIN vtypes
	--WHERE sampleid = 'CPCT02020314T'
),
selected_variants AS
(	
	SELECT
		sampleid,
		eventid,
		hwg.chrom,
		hwg.pos,
		arm,
		variant_type,
		genotype
	FROM hwg
	INNER JOIN coverage_by_arm sa ON hwg.chrom = sa.chrom AND hwg.pos && sa.pos -- annotate arm
	WHERE
		genotype != '0/0' --AND sampleid = 'CPCT02020314T'
),
arm_counts AS
(
	SELECT
		sa.sampleid,
		sa.arm,
		sa.variant_type,
		sa.arm_call,
		sa.arm_ploidy,
		sa.covered_arm_size,
		COUNT(sv.genotype) AS abs_arm_count,
		ROUND((COUNT(sv.genotype)/covered_arm_size::numeric) * 1e6,4) AS rel_arm_freq
	FROM sample_arms sa
	LEFT JOIN selected_variants sv ON sv.sampleid = sa.sampleid AND sv.arm = sa.arm AND sv.variant_type = sa.variant_type
	GROUP BY 1,2,3,4,5,6
)
SELECT *
FROM arm_counts
ORDER BY 1,2,3
-- END --