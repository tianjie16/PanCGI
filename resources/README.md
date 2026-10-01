# HPRC2 PanCGI resource

This resource comprises 47,039 CpG island (CGI) alleles at 33,463 loci, 12,639,213 assembly-specific members, representative sequences, CGI state matrices and evolutionary trajectories.

## Files

| Directory | Contents |
| --- | --- |
| `catalogue/` | Loci, alleles, members and allele-haplotype membership |
| `genotypes/` | Allele and locus CGI states for 462 HPRC haplotypes and CHM13/GRCh38 |
| `sequences/` | One representative sequence per allele |
| `evolution/` | Trajectory classifications for 2,279 loci |
| `calls/masked/` | Masked CGI calls for 462 HPRC haplotypes, including HG002 |
| `calls/unmasked/` | Unmasked CGI calls for the same haplotypes |
| `metadata/` | Assembly sources, contig lengths, field definitions and BED file inventory |

## Catalogue

The allele catalogue corresponds to Supplementary Table 9 (`HPRC2_catalogue=1`). Join alleles to loci by `locus_id`, members to alleles by `allele_id`, and haplotypes to `metadata/assemblies.tsv` by `haplotype_id`. Each member is an assembly-specific CGI call; an allele groups similar member sequences. FASTA record names are allele IDs, also used as representative member IDs. Member coordinates and assembly sources identify individual member sequences.

Representative and member assembly intervals are 0-based, half-open. CHM13 locus coordinates and associated-SV positions use the 1-based convention of Supplementary Table 9. Trajectory coordinates follow Supplementary Table 10.

FASTA sequences are in the forward orientation of the source assembly. The `representative_sequence_*` fields describe these bases, including ambiguous bases and assembly gaps. CGI-call statistics describe the masked sequence used for calling; FASTA statistics describe the underlying assembly sequence.

## States and frequencies

State matrices use `1` for CGI-positive/allele-present, `0` for CGI-negative/allele-absent and `.` for a missing state. Other tables use empty fields for missing or inapplicable values.

`fraction_of_462_HPRC_haplotypes` uses all 462 HPRC haplotypes. The locus field `frequency_among_called_HPRC_haplotypes` uses called haplotypes. Reference assemblies are excluded from these denominators.

`CHM13_CGI` and `GRCh38_CGI` indicate CGI presence at a locus; `*_CGI_in_released_alleles` indicate reference membership within the allele catalogue. Reference CGI alleles can lie outside this catalogue. Reference members belonging to catalogue alleles are included.

The trajectory table contains 2,034 resolved classifications and 245 loci with an unresolved root state or unavailable genealogy. Empty trajectory labels in the catalogue indicate loci outside this analysis.

## CGI calls

Both call sets use UCSC `cpg_lh`. Masked calls use RepeatMasker and TRF period-1/2 masking. The BED sets contain all calls for the 462 HPRC haplotypes; catalogue membership is given in `catalogue/`.

BED files use 0-based, half-open coordinates and the 10-column UCSC cpgIslandExt layout (BED4+6): `chrom, start0, end0, name, length, cpgNum, gcNum, perCpg, perGc, obsExp`. Column 5 is interval length; column 6 is CpG count. Member IDs use `chrom:start0-end0`; `CpG: N` labels describe CpG counts.

## Download

The six archives extract into `PanCGI_resource_v1.0.0/`. Verify downloaded files using the accompanying `MD5SUMS.txt`; after extraction, run `md5sum -c MD5SUMS.txt` within the resource directory to verify its contents.

`SCHEMA.json` and `metadata/FIELDS.tsv` define the fields. `metadata/assemblies.tsv` provides assembly names, accessions and source URLs. `FILES.tsv` lists file sizes and MD5 checksums. Resource version: 1.0.0.

## Sources

- [HPRC assembly index](https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/41aa47dd3430fbb250cdb6a78efde43313d35557/data_tables/assemblies_pre_release_v0.6.1.index.csv)
- [UCSC cpgIslandExt schema](https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=regulation&hgta_track=cpgIslandExt&hgta_table=cpgIslandExt&hgta_doSchema=describe+table+schema)
- [BED coordinates](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
