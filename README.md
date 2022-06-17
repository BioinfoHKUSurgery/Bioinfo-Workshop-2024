# Pipeline for evaluating actionable findings from VCF file

## Prerequisites
### _Software needed_
- GATK v4
- bcftools

### _Data_
- AnVIL_CCDG_Broad_NP_Epilepsy_HKOSB_GRU_WES_Year4.vcf.gz.*

### _Resources_
- GATK GRCh38 resource files (via FTP Server Access; see details [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle))

## Create a pedigree file
- Create a pedigree file with six columns
1. Family ID (same as individual ID if unrelated)
2. Individual ID
3. Paternal ID (0 if founder)
4. Maternal ID (0 if founder)
5. Sex (2=Female; 1=Male)
6. Phenotype (2=Case; 1=Control)

```
zcat AnVIL_CCDG_Broad_NP_Epilepsy_HKOSB_GRU_WES_Year4.vcf.gz | head -5000 | egrep CHROM | \
  awk '{ for (i=10;i<=NF;++i) print $i,$i,0,0,2,1}' \
  > AnVIL_CCDG_Broad_NP_Epilepsy_HKOSB_GRU_WES_Year4.vcf.ped
```

## Prepare a dbSNP file by excluding variants after version 129
```bash
# specify the path to your gatk resource folder 
gatk_resources=<path>

refseq=$gatk_resources/Homo_sapiens_assembly38.fasta
dbsnp=$gatk_resources/dbsnp_146.hg38.vcf.gz

# remove variants after version 129 to evaluate TiTv (>3.0 for exome and >2.0 for genome)
bcftools filter -e "dbSNPBuildID>129" /psychipc01/disk2/references/GATK_hg38_bundle/dbsnp_146.hg38.vcf.gz | bgzip > $gatk_resources/dbsnp_146.hg38.excluding_sites_after_129.vcf.gz
tabix -p vcf -f $gatk_resources/dbsnp_146.hg38.excluding_sites_after_129.vcf.gz
dbsnp129=$gatk_resources/dbsnp_146.hg38.excluding_sites_after_129.vcf.gz
```

## Variant Evaluation using GATK to determine the TiTv
```bash
./gatk --java-options "-Xms10g -Xmx10g" VariantEval \
   -R $refseq \
   --dbsnp $dbsnp129 \
   --eval set1:AnVIL_CCDG_Broad_NP_Epilepsy_HKOSB_GRU_WES_Year4.vcf.gz \
   -ped AnVIL_CCDG_Broad_NP_Epilepsy_HKOSB_GRU_WES_Year4.ped \
   -noST -ST Sample -ST Novelty -ST Filter \
   -noEV -EV CountVariants \
   -noEV -EV TiTvVariantEvaluator \
   -S LENIENT \
   -O AnVIL_CCDG_Broad_NP_Epilepsy_HKOSB_GRU_WES_Year4.VariantEval.gatk-report
```
