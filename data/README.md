Datafiles:

* all_dnsv.csv: variant info for the *de novo* SVs based on the GRCh38 reference genome.
 <details>
 * Chrom: chromosome of the SV
 * Pos: start position of the SV
 * End: end position of the SV
 * Project: CEPH or SFARI
 * Sample: sample ID (de-identified)
 * SVLen: length of the SV
 * SVType: type of the SV
 * Call_pipeline: derivation of the original call (Smoove or GATK-SV)
 * Father: ID of the father
 * Father_age_birth_years: age of the father at birth of offspring
 * Generation: only for CEPH (SV from the F1 or F2 generation)
 * Mosaic: True if SV labeled mosaic
 * Mother: ID of the mother
 * Mother_age_birth_years: age of the mother at birth of offspring
 * Parent_of_origin: NA if unidentified
 * Pred_mechanism: predicted causal mechanism of the SV
 * Role: sibling or proband
 * Sex: male or female
 * DN_SNVs: count of *de novo* SNVs in sample with the *de novo* SV
 </details>

* all_dnsv_metadata.csv: some of the above fields (dn_snvs, parent info, sex) but for all samples rather than only those with dnSVs
* dnsv_dataframe.csv: dnSV data summarized for creation of figure plots in R script. Count of dnSVs by type in all samples, count of dnSNVs where known samples, parent information
* hg38.blacklist.bed: regions not included (centromeres, telomeres, chrM, PARs, unmappable regions, HLAs, etc.)
* size_freqs.csv: brief summary of size frequencies for R analysis of length enrichment in probands
