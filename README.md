# Pipeline for evaluating actionable findings from VCF file

## Prerequisites
### _Software needed_
- GATK v4
- bcftools
- PLINK2 (https://www.cog-genomics.org/plink/2.0/)

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
   --eval AnVIL_CCDG_Broad_NP_Epilepsy_HKOSB_GRU_WES_Year4.vcf.gz \
   -ped AnVIL_CCDG_Broad_NP_Epilepsy_HKOSB_GRU_WES_Year4.ped \
   --no-st -ST Sample -ST Novelty -ST Filter \
   --no-ev -EV CountVariants -EV TiTvVariantEvaluator \
   --lenient \
   -O AnVIL_CCDG_Broad_NP_Epilepsy_HKOSB_GRU_WES_Year4.VariantEval.gatk-report
```

## Individual-based and variant-based quality control using PLINK
```bash
# PLINK2
plink2 --vcf /lustre1/u/u3579068/project/AnVIL_CCDG_Broad_NP_Epilepsy_HKOSB_GRU_WES_Year4_AnVIL_CCDG_Broad_NP_Epilepsy_HKOSB_GRU_WES_Year4.vcf.gz \
  --double-id \
  --vcf-min-gq 20 \
  --vcf-min-dp 8 \
  --out AnVIL_CCDG_Broad_NP_Epilepsy_HKOSB_GRU_WES_Year4.QCpos

plink2 --pfile AnVIL_CCDG_Broad_NP_Epilepsy_HKOSB_GRU_WES_Year4.QCpos \
  --maf 0.05 \
  --geno 0.05 \
  --max-alleles 2 \
  --make-bed \
  --out AnVIL_CCDG_Broad_NP_Epilepsy_HKOSB_GRU_WES_Year4.QCedpos.biallelic.maf05

# PLINK 1.9
plink --bfile AnVIL_CCDG_Broad_NP_Epilepsy_HKOSB_GRU_WES_Year4.QCedpos.biallelic.maf05 \
  --indep-pairwise 200 50 0.2 \
  --out AnVIL_CCDG_Broad_NP_Epilepsy_HKOSB_GRU_WES_Year4.QCedpos.biallelic.maf05.pruned 

plink --bfile AnVIL_CCDG_Broad_NP_Epilepsy_HKOSB_GRU_WES_Year4.QCedpos.biallelic.maf05 \
  --extract AnVIL_CCDG_Broad_NP_Epilepsy_HKOSB_GRU_WES_Year4.QCedpos.biallelic.maf05.pruned.prune.in \
  --genome \
  --out AnVIL_CCDG_Broad_NP_Epilepsy_HKOSB_GRU_WES_Year4.QCedpos.biallelic.maf05.pruned 
```
# KGGseq annotation
- Allele frequency
- In silico damaging scores
```bash
java -Xmx37g -jar kggseq.jar \
        --resource resources/ \
        --no-lib-check --buildver hg19 \
        --vcf-file /lustre1/u/u3579068/project/AnVIL_CCDG_Broad_NP_Epilepsy_HKOSB_GRU_WES_Year4_AnVIL_CCDG_Broad_NP_Epilepsy_HKOSB_GRU_WES_Year4.vcf.gz \
        --seq-qual 30 --seq-mq 20 --gty-qual 20 --gty-dp 8 \
        --db-gene refgene --gene-feature-in 0,1,2,3,4,5,6,7,8,9,10 \
        --db-filter 1kg201204,1kgafr201204,1kgeur201204,1kgasn201204,1kgeur201305,1kgeas201305,1kgafr201305,ESP6500AA,ESP6500EA,exac,ehr,gadexome,gadgenome \
        --rare-allele-freq 1.1 \
        --db-score dbnsfp \
        --mendel-causing-predict best \
        --out AnVIL_CCDG_Broad_NP_Epilepsy_HKOSB_GRU_WES_Year4.impact \
        --nt 10       # number of nodes
```
