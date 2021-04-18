/*
---
Date: May 1, 2020
Author: E. Kocakavuk
---
---
Query deletion size from HMF data
---
*/
SELECT sampleid, UPPER(pos) - LOWER (pos) -1 AS length, ref, alt, variant_type 
FROM hwg
WHERE variant_type = 'DEL' 