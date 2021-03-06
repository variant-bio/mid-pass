# Methods for mid-pass whole genome sequencing

This document provides an overview of how the sequencing data in  “Mid-pass whole genome sequencing enables biomedical genetic studies of diverse populations” were processed. 

This code is freely available for academic and non-commercial research purposes [LICENSE](LICENSE.md).

## Processing workflow to generate fully phased and imputed VCF

<img src="workflow.png" alt="Workflow overview" width="800"/>
 
 
## Preprocessing
Sample level processing down to GVCFs followed the GATK Best Practices guidelines as detailed [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-). Specifically, after adapter trimming with cutadapt v2.10, reads were mapped against GRCh38 with ALT contigs using BWA-mem v0.7.15. Subsequent steps for duplicate marking, sorting, base call quality recalibration and creating GVCFs with HaplotypeCaller were done with versions v2.14 of Picard, v1.10 of samtools, and v4.1.4.0 of GATK. Most steps were run with default parameter settings with few exceptions, e.g. to remove extremely short alignments. 

## Joint variant calling and site level filtering
Mid-pass (1-7x) GVCFs and 30x GVCFs were jointly genotyped (using GATK genomicsDBImport and GenotypeGVCFs) to produce one multi-sample VCF that was subsequently run through VQSR (SNP sensitivity=99.8%, indel sensitivity=99%) for site-level filtering. Only PASS sites were used subsequently.

## GQ filtering and imputation
The PASS-only variant calls were further filtered to remove low quality genotype calls, to exclude those from potentially misguiding imputation.  By measuring imputation performance over varying GQ thresholds, we determined that keeping only genotypes with GQ>17 was optimal. All GT fields with GQ<=17 were hence set to "./." (this can be done for example with vcftools --minGQ (see [here](https://vcftools.github.io/man_latest.html)).
The resulting VCF file was then used as input to Beagle 5.1 using the “gt=” input parameter without reference panel (i.e. to impute within the provided cohort).

## Call flagging
We finally compared the imputed genotypes with the sequencing-based calls before GQ filtering, in order to identify and flag calls where sequencing data (even if at low confidence) was in agreement or disagreement with the imputed call. The script for merging unfiltered and imputed VCF files is provided here ([flag_calls.py](flag_calls.py)). It annotates each genotype call with an IM value, where a lower values represents higher confidence in the call. Depending on the downstream application, these flags can be used to treat certain flagged variants with caution or filter them entirely.


<img src="IMflagging.png" alt="IM flagging overview" width="550"/>

A call is flagged with IM=0 if sequencing-based genotype and imputed genotype agree fully. IM=1 means the imputed call does not disagree with the sequencing-based call (either because it was missing or we may have only observed one of two alleles in sequencing). IM=2 and IM=3 flag sites with disagreement between sequencing-based and imputed calls, where IM=2 calls were heterozygous in sequencing (potentially due to sequencing errors or contamination) and IM=3 calls were homozygous. Especially IM=3 calls where sequencing and imputation disagree most strongly may better be excluded or further refined for downstream analysis. 


## Performance assessment

We assessed recall and precision for called variant sites and non-reference concordance rate (NCR = 1-NDR) of genotype calls for each test sample. The following code snippets detail the commands used to subset test and "truth" VCF files using bcftools (v1.10). This was done for each sample, variant type (snps|indels) and within HC regions or otherwise genomewide.

```# extract sample “truth” calls
bcftools view [-R $regions] -s $sample -v $type joint.truth.vcf.gz | bcftools view -i ‘FMT/GQ>20’ | bgzip -c > $sample.truth.vcf.gz
bcftools view -c 1 $sample.truth.vcf.gz | bgzip -c > $sample.truth.sites.vcf.gz

## extract sample test calls - different cases for different experiments
# 1. from imputed VCF
bcftools view [-R $regions] -s $sample -v $type joint.imputed.vcf.gz | bgzip -c > $sample.vcf.gz

# 2. from jointly genotyped PASS filter VCF
bcftools view [-R $regions] -s $sample -v $type joint.vcf.gz | bcftools view -i ‘FMT/GQ>17’ | bgzip -c > $sample.vcf.gz

# 3. from single sample VCF
bcftools view [-R $regions] -v $type $sample.single.vcf.gz | bcftools view -i ‘FMT/GQ>17’ | bgzip -c > $sample.vcf.gz
```

The different test VCFs were then compared to the corresponding sample truth VCFs to examine site recall and precision as well as genotype accuracy, using vcf-compare (v0.1.14-12-gcdb80b8).

```
# get only those sites where sample has a non-ref allele
bcftools view -c 1 $sample.vcf.gz > $sample.sites.vcf.gz

# site comparison for recall/precision 
vcf-compare $sample.truth.sites.vcf.gz $sample.sites.vcf.gz | grep ^VN

# genotype comparison for NCR = 1-NDR
vcf-compare -g $sample.truth.vcf.gz $sample.vcf.gz | grep -A 1 ^GS
```
## Publicly available mid-pass data for benchmarking

Sequencing data from 120 libraries generated for the 12 HapMap individuals has been deposited to the Sequence Read Archive under acessession [PRJNA697982](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA697982/). Briefly, 4 replicate libraries were generated using the [plexWell LPS384](https://seqwell.com/products/plexwell-lp-384/) library preparation kit (designed for low-pass applications), and 2 replicate libraries were generated using the [plexWell WGS24](https://seqwell.com/products/plexwell-wgs-24/) library preparation kit (designed for high-pass applications). LPS384 libraries were sequenced at target coverages of 1x and 4x, and WGS24 libraries were sequenced at a target coverage of 4x.
