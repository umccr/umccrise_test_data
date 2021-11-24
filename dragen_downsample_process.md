# Downsampling DRAGEN data
## TODO
* check requirement for other features
    - e.g. hla, viral integration detection
* get read mates
* spike in viral reads
    - requires remapping
    - alternatively, use a sample that has existing viral reads or calls
* remove unused sites/variants e.g. peddy
* reference data subset/downsampling if testing needs to be expanded
* may need to downsample tumour SV VCFs depending on variant count

## Setup
Directories
```bash
mkdir -p \
  data/{reference,sample}/ \
  output/SEQC50/germline/2021-11-22--3.9.3/SEQC-II_dragen_germline/ \
  output/SEQC50/somatic/2021-11-22--3.9.3/SEQC-II_dragen_somatic/ \
  working/
```

Data
```bash
PREDISPOSE_GENES_URL='https://raw.githubusercontent.com/umccr/NGS_Utils/master/ngs_utils/reference_data/key_genes/predispose_genes.hg38.coding.bed'
PREDISPOSE_GENES_FP=data/reference/"${PREDISPOSE_GENES_URL##*/}"
curl "${PREDISPOSE_GENES_URL}" > "${PREDISPOSE_GENES_FP}"

KEY_GENES_URL='https://raw.githubusercontent.com/umccr/NGS_Utils/master/ngs_utils/reference_data/key_genes/umccr_cancer_genes.hg38.coding.bed'
KEY_GENES_FP=data/reference/"${KEY_GENES_URL##*/}"
curl "${KEY_GENES_URL}" > "${KEY_GENES_FP}"

PEDDY_SITES_URL='https://raw.githubusercontent.com/brentp/peddy/master/peddy/GRCH38.sites'
PEDDY_SITES_FP=data/reference/peddy_"${PEDDY_SITES_URL##*/}"
curl "${PEDDY_SITES_URL}" > "${PEDDY_SITES_FP}"

CONPAIR_SITES_URL='https://raw.githubusercontent.com/nygenome/Conpair/master/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed'
CONPAIR_SITES_FP=data/reference/conpair_sites_hg38.bed
curl "${CONPAIR_SITES_URL}" > "${CONPAIR_SITES_FP}"

# umccrise reference genome index
GENOME_FA_FAI_URI='s3://umccr-refdata-dev/gpl-nf/genome/umccrise_hg38/hg38.fa.fai'
GENOME_FA_FAI_FP=data/reference/hg38.fa.fai
aws s3 cp "${GENOME_FA_FAI_URI}" "${GENOME_FA_FAI_FP}"
```

Config
```bash
SLOP=200
```

> DRAGEN data downloaded from GDS: `gds://development/validation_data/wgs/SEQC50/analysis/`


## Germline VCF
Process:
* combine sites used by normal workflow
    * peddy sites
    * conpair sites
    * predispose gene sites
* add slop at boundaries
* extract these regions of interest and index new VCF
* create new BED containing regions of only extracted elements

```bash
NORMAL_VCF_FP=data/sample/SEQC50/germline/2021-11-22--3.9.3/SEQC-II_dragen_germline/SEQC-II.hard-filtered.vcf.gz
NORMAL_VCF_DS_FP=output/"${NORMAL_VCF_FP##data/sample/}"

CONPAIR_SITES_MIN_FP=working/"${CONPAIR_SITES_FP##*/}"
PEDDY_SITES_BED_FP=working/"${PEDDY_SITES_FP##*/}.bed"
PREDISPOSE_GENES_MIN_FP=working/"${PREDISPOSE_GENES_FP##*/}"
NORMAL_BED_FP=working/normal.bed
NORMAL_VARIANTS_BED_FP=working/normal_variants.bed

shuf -n500 --random-source=/dev/zero "${PEDDY_SITES_FP}" | \
  awk -F: '{print "chr" $1 "\t" $2-1 "\t" $2}' > "${PEDDY_SITES_BED_FP}";
shuf -n500 --random-source=/dev/zero "${CONPAIR_SITES_FP}" > "${CONPAIR_SITES_MIN_FP}"
cut -f1-3 -d$'\t' "${PREDISPOSE_GENES_FP}" > "${PREDISPOSE_GENES_MIN_FP}";

cat "${CONPAIR_SITES_MIN_FP}" "${PEDDY_SITES_BED_FP}" "${PREDISPOSE_GENES_MIN_FP}" | \
  bedtools sort | \
  bedtools merge | \
  bedtools slop -b "${SLOP}" -g "${GENOME_FA_FAI_FP}" > "${NORMAL_BED_FP}"

bcftools view -Oz -R "${NORMAL_BED_FP}" "${NORMAL_VCF_FP}" > "${NORMAL_VCF_DS_FP}"

tabix "${NORMAL_VCF_DS_FP}"
bcftools query -f '%CHROM\t%POS\n' "${NORMAL_VCF_DS_FP}" | \
  awk '{print $1 "\t" $2-1 "\t" $2}' | \
  bedtools slop -b "${SLOP}" -g "${GENOME_FA_FAI_FP}" > "${NORMAL_VARIANTS_BED_FP}"
```

## Somatic VCF
Process:
* combined sites in key genes and selected conpair sites
* add slop at gene bounds
* extract these regions of interest and index new VCF
* create new BED containing regions of only extracted elements

```bash
TUMOUR_VCF_FP=data/sample/SEQC50/somatic/2021-11-22--3.9.3/SEQC-II_dragen_somatic/SEQC-II.hard-filtered.vcf.gz
TUMOUR_VCF_DS_FP=output/"${TUMOUR_VCF_FP##data/sample/}"
KEY_GENES_MIN_FP=working/"${KEY_GENES_FP##*/}"

TUMOUR_BED_FP=working/tumour.bed
TUMOUR_VARIANTS_BED_FP=working/tumour_variants.bed

cut -f1-3 -d$'\t' "${KEY_GENES_FP}" > "${KEY_GENES_MIN_FP}";

cat "${CONPAIR_SITES_MIN_FP}" "${KEY_GENES_MIN_FP}" | \
  bedtools sort | \
  bedtools merge | \
  bedtools slop -b "${SLOP}" -g "${GENOME_FA_FAI_FP}" > "${TUMOUR_BED_FP}"
bcftools view -Oz -R "${TUMOUR_BED_FP}" "${TUMOUR_VCF_FP}" > "${TUMOUR_VCF_DS_FP}"

tabix "${TUMOUR_VCF_DS_FP}"
bcftools query -f '%CHROM\t%POS\n' "${TUMOUR_VCF_DS_FP}" | \
  awk '{print $1 "\t" $2-1 "\t" $2}' | \
  bedtools slop -b "${SLOP}" -g "${GENOME_FA_FAI_FP}" > "${TUMOUR_VARIANTS_BED_FP}"
```

## BAMs
Process
* take all selected sites (peddy, conpair, etc) and sites of downsampled variants
* subset BAMs these sites
* index BAMs
```bash
TUMOUR_BAM_FP=data/sample/SEQC50/somatic/2021-11-22--3.9.3/SEQC-II_dragen_somatic/SEQC-II_tumor.bam
NORMAL_BAM_FP=data/sample/SEQC50/germline/2021-11-22--3.9.3/SEQC-II_dragen_germline/SEQC-II.bam
TUMOUR_BAM_DS_FP=output/"${TUMOUR_BAM_FP##data/sample/}"
NORMAL_BAM_DS_FP=output/"${NORMAL_BAM_FP##data/sample/}"
BAM_BED_FP=working/bam.bed

cat "${TUMOUR_VARIANTS_BED_FP}" "${NORMAL_VARIANTS_BED_FP}" "${CONPAIR_SITES_MIN_FP}" "${PEDDY_SITES_BED_FP}" | \
  bedtools sort | \
  bedtools merge > "${BAM_BED_FP}"

(sambamba view -f bam -t 2 -L "${BAM_BED_FP}" "${TUMOUR_BAM_FP}" > "${TUMOUR_BAM_DS_FP}") &
(sambamba view -f bam -t 2 -L "${BAM_BED_FP}" "${NORMAL_BAM_FP}" > "${NORMAL_BAM_DS_FP}") &
wait
samtools index "${TUMOUR_BAM_DS_FP}" &
samtools index "${NORMAL_BAM_DS_FP}" &
wait
```

## Others (not downsampled)
```bash
TUMOUR_DIR=output/SEQC50/somatic/2021-11-22--3.9.3/SEQC-II_dragen_somatic/
NORMAL_DIR=output/SEQC50/germline/2021-11-22--3.9.3/SEQC-II_dragen_germline/

cp data/sample/SEQC50/somatic/2021-11-22--3.9.3/SEQC-II_dragen_somatic/SEQC-II.sv.vcf.gz "${TUMOUR_DIR}"
cp data/sample/SEQC50/somatic/2021-11-22--3.9.3/SEQC-II_dragen_somatic/SEQC-II-replay.json "${TUMOUR_DIR}"
cp data/sample/SEQC50/germline/2021-11-22--3.9.3/SEQC-II_dragen_germline/SEQC-II-replay.json "${NORMAL_DIR}"

samtools view -H data/sample/SEQC50/somatic/2021-11-22--3.9.3/SEQC-II_dragen_somatic/SEQC-II_normal_normal.bam > "${TUMOUR_DIR}/SEQC-II_normal_normal.bam"

qc_files=".vc_metrics.csv .ploidy_estimation_metrics.csv .mapping_metrics.csv .fragment_length_hist.csv .fastqc_metrics.csv"
qc_files_tumour=".wgs_fine_hist_tumor.csv .wgs_coverage_metrics_tumor.csv .wgs_contig_mean_cov_tumor.csv"
qc_files_normal=".wgs_fine_hist.csv .wgs_coverage_metrics.csv .wgs_contig_mean_cov.csv"

for file_suffix in ${qc_files} ${qc_files_tumour}; do
  file=$(find data/sample/SEQC50/somatic/2021-11-22--3.9.3/SEQC-II_dragen_somatic/ -name "*${file_suffix}")
  file_n=$(echo "${file}" | wc -l)
  if [[ -z ${file} ]]; then
    echo "no matches found for ${file_suffix}"
  elif [[ ${file_n} -ne 1 ]]; then
    echo "got unexpected number of files (${file_n}) for ${file_suffix}: ${file}"
  else
    cp "${file}" "${TUMOUR_DIR}"
  fi
done

for file_suffix in ${qc_files} ${qc_files_normal}; do
  file=$(find data/sample/SEQC50/germline/2021-11-22--3.9.3/SEQC-II_dragen_germline/ -name "*${file_suffix}")
  file_n=$(echo "${file}" | wc -l)
  if [[ -z ${file} ]]; then
    echo "no matches found for ${file_suffix}"
  elif [[ ${file_n} -ne 1 ]]; then
    echo "got unexpected number of files (${file_n}) for ${file_suffix}: ${file}"
  else
    cp "${file}" "${NORMAL_DIR}"
  fi
done
```

## Organise
```bash
mkdir -p final/
rsync -aP output/SEQC50/somatic/2021-11-22--3.9.3/SEQC-II_dragen_somatic final/
rsync -aP output/SEQC50/germline/2021-11-22--3.9.3/SEQC-II_dragen_germline final/
```
