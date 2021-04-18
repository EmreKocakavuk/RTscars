/*
---
Date: April 24, 2020
Author: FP Barthel
Update: December 26, 2020 (by Emre Kocakavuk)
---
Count somatic variants (indels + SNVs) from the Hartwig dataset
- Aggregated counts by sample
- Age and treatment are annotated
- Repeat region info added
---
*/
SELECT
	hwg.sampleid,
	hwm.patientid,
	variant_type,
	--reptype
	COUNT(*)::integer nmut
FROM hwg
INNER JOIN hmf_meta hwm ON hwm.sampleid = hwg.sampleid
--LEFT JOIN ref.repeatmasker rm ON hwg.chrom = rm.chrom AND hwg.pos && rm.pos
WHERE variant_type IN ('SNP','INS','DEL') --= 'DEL' AND rm.reptype IS NULL
GROUP BY 1,2,3