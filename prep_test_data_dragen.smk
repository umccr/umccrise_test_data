import os, subprocess
from ngs_utils.dragen import DragenProject
from os.path import join, abspath, relpath, basename, isfile, dirname
from ngs_utils.reference_data import get_key_genes_bed
from ngs_utils import logger as log
from reference_data import api as refdata


shell.executable(os.environ.get('SHELL', 'bash'))
log.init(config.get('debug', False))

GENOME = 'hg38'
#DRAGEN_DIR = config.get('dragen_dir', '/g/data/gx8/extras/umccrise_017_2020_Jan/umccrise_test_data/running_dragen/T_SRR7890902_20pc')
DRAGEN_DIR = config.get('dragen_dir', '/g/data/gx8/extras/vlad/synced/umccr/umccrise_test_data/running_dragen/P025')
run = DragenProject(DRAGEN_DIR)
included_names = [s.name for s in run.samples]
batch_by_name = {b.tumor.name: b for b in run.batch_by_name.values() if not b.is_germline()}
assert included_names

PROJECT_NEW_PATH = config.get('out', 'data/dragen_test_project')

ONCOVIRAL_READS_DIR = join('data', 'viral_reads_from_neverresponder')
Out_PON_PATH = f'data/genomes/{GENOME}/panel_of_normals'
SORT_SEED = "--random-source <(echo AAAAAAAAAAAAAAAAAAA)"


def extract_genes_from_vcf(inp='{input}', out='{output}'):
    key_genes_bed = get_key_genes_bed(GENOME)
    return f'bcftools view -R {key_genes_bed} {inp} > {out} && grep -v ^# {out} | wc'

def vcf_to_bed():
    return "bcftools view -H {input} | py -x \"'\\t'.join([x.split()[0], str(int(x.split()[1])-1), x.split()[1]])\"" \
           " > {output}"

def downsample_vcf(num, inp='{input}', out='{output}'):
    return 'bcftools view -h ' + inp + ' > ' + out + ' && ' \
           'bcftools view -H ' + inp + ' -f.,PASS | sort -R ' + SORT_SEED + ' | head -n' + str(num) + '' \
          ' | sort -k1,1V -k2,2n >> ' + out


rule all:
    input:
        expand(join(PROJECT_NEW_PATH, '{batch}.bam'), batch=batch_by_name.keys()),
        expand(join(PROJECT_NEW_PATH, '{batch}.bam.bai'), batch=batch_by_name.keys()),
        expand(join(PROJECT_NEW_PATH, '{batch}_tumor.bam'), batch=batch_by_name.keys()),
        expand(join(PROJECT_NEW_PATH, '{batch}_tumor.bam.bai'), batch=batch_by_name.keys()),
        expand(join(PROJECT_NEW_PATH, '.populated_downsampled_{batch}.done'), batch=batch_by_name.keys()),
        expand(join(PROJECT_NEW_PATH, '.populated_other_files_{batch}.done')   , batch=batch_by_name.keys()),
        'work_snake/roi.bed',
        # reference data:
        # pon_snps_vcf    = join(Out_PON_PATH, 'panel_of_normals.snps.vcf.gz'),
        # pon_indels_vcf  = join(Out_PON_PATH, 'panel_of_normals.indels.vcf.gz'),
        # gnomad          = f'data/genomes/{GENOME}/gnomad_genome.vcf.gz',
        # purple_gc       = f'data/genomes/{GENOME}/hmf/GC_profile.1000bp.cnp',
        # purple_het      = f'data/genomes/{GENOME}/hmf/germline_het_pon.vcf.gz',
        # hmf_hotspot     = f'data/genomes/{GENOME}/hmf/KnownHotspots.tsv.gz',
        # hmf_giab_conf   = f'data/genomes/{GENOME}/hmf/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed.gz'


######################################
#### SOMATIC ####
rule downsample_somatic:
    input:
        lambda wc: batch_by_name[wc.batch].somatic_vcf
    output:
        'work_snake/somatic/{batch}.vcf'
    shell:
        extract_genes_from_vcf()

rule somatic_roi:
    input:
        rules.downsample_somatic.output[0]
    output:
        bed = 'work_snake/somatic/{batch}.bed'
    shell:
        vcf_to_bed()

# ######################################
# #### GERMLINE ####
# rule downsample_germline:
#     input:
#         lambda wc: join(run.date_dir, f'{batch_by_name[wc.batch].normal.name}-germline-ensemble-annotated.vcf.gz')
#     output:
#         'work_snake/germline/{batch}-germline-ensemble.vcf'
#     shell:
#         extract_genes_from_vcf()
#
# rule downsample_germline_random100:
#     input:
#         rules.downsample_germline.output[0]
#     output:
#         'work_snake/germline/{batch}-germline-ensemble-100.vcf'
#     shell:
#         downsample_vcf(100)
#
# rule germline_roi:
#     input:
#         rules.downsample_germline_random100.output[0]
#     output:
#         'work_snake/germline/{batch}-germline-ensemble.bed'
#     shell:
#         vcf_to_bed()


######################################
#### SV VCF ####
rule downsample_sv:
    input:
        lambda wc: batch_by_name[wc.batch].sv_vcf
    output:
        'work_snake/sv/{batch}.vcf'
    shell:
        'bcftools view -f.,PASS {input} -o {output}'

rule sv_roi:
    input:
        rules.downsample_sv.output[0]
    output:
        bed = 'work_snake/sv/{batch}.bed'
    shell:
        vcf_to_bed()

######################################
#### BAMS ####
rule batch_roi:
    input:
        bed = rules.somatic_roi.output.bed,
        # rules.germline_roi.output[0]
    output:
        bed = 'work_snake/{batch}-roi.bed'
    shell:
        'cat {input.bed} | bedtools sort -i stdin | bedtools merge -i stdin > {output}'

# adding few chrY-unique locations to make goleft happy
# rule sex_bed:
#     output:
#         'work_snake/sex.bed'
#     run:
#         sex_bed = 'Y\t2655029\t2655030\n' + \
#                   'Y\t2713681\t2713682\n' + \
#                   'Y\t6939595\t6939596\n'
#         with open(output[0], 'w') as out:
#             out.write(sex_bed)

rule conpair_roi:
    input:
        bed = 'data/conpair_markers/GRCh37.bed' if GENOME == 'GRCh37' else \
              'data/conpair_markers/hg38.liftover.bed',
        fai = refdata.get_ref_file(GENOME, key='fa') + '.fai'
    output:
        bed = 'work_snake/conpair_roi.bed'
    shell:
        'head -n100 {input.bed} '
        ' | bedtools slop -b 1000 -i stdin -g {input.fai}'
        '> {output}'

rule project_roi:
    input:
        beds = expand(rules.batch_roi.output.bed, batch=batch_by_name.keys()) +
               [rules.conpair_roi.output.bed],
        fai = refdata.get_ref_file(GENOME, key='fa') + '.fai'
    output:
        'work_snake/roi.bed'
    shell:
        'cat {input.beds}'
        ' | bedtools sort -i stdin'
        ' | bedtools slop -b 10 -i stdin -g {input.fai}'
        ' | bedtools merge -i stdin'
        ' > {output}'

rule project_somatic_roi:
    input:
        beds = expand(rules.somatic_roi.output.bed, batch=batch_by_name.keys()),
        fai = refdata.get_ref_file(GENOME, key='fa') + '.fai',
    output:
        bed = 'work_snake/somatic_roi.bed'
    shell:
        'cat {input.beds}'
        ' | bedtools sort -i stdin'
        ' | bedtools slop -b 10 -i stdin -g {input.fai}'
        ' | bedtools merge -i stdin'
        ' > {output}'

rule subset_bam:
    input:
        bam = lambda wc: getattr(batch_by_name[wc.batch], wc.phenotype).bam,
        roi_bed = rules.project_roi.output[0]
    output:
        bam_namesorted = 'work_snake/subset_bam/{batch}_{phenotype}.bam'
    shell:
        'sambamba slice {input.bam} -L {input.roi_bed}'
        ' | samtools sort -n -Obam -o {output}'

rule extract_reads:
    input:
        bam_namesorted = rules.subset_bam.output.bam_namesorted,
    output:
        fq1 = 'work_snake/bam_subset_fq/{batch}_{phenotype}.R1.fq',
        fq2 = 'work_snake/bam_subset_fq/{batch}_{phenotype}.R2.fq',
    shell:
        'samtools fastq {input.bam_namesorted} -1 {output.fq1} -2 {output.fq2} -s /dev/null'

rule add_viral_reads:
    input:
        fq1 = rules.extract_reads.output.fq1,
        fq2 = rules.extract_reads.output.fq2,
        viral_fq1 = join(ONCOVIRAL_READS_DIR, 'step6_HPV18_bridging.R1.fq'),
        viral_fq2 = join(ONCOVIRAL_READS_DIR, 'step6_HPV18_bridging.R2.fq'),
    output:
        fq1 = 'work_snake/bam_subset_plusviral_fq/{batch}_{phenotype}.plusviral.R1.fq',
        fq2 = 'work_snake/bam_subset_plusviral_fq/{batch}_{phenotype}.plusviral.R2.fq',
    shell:
        'cat {input.fq1} {input.viral_fq1} > {output.fq1} &&'
        'cat {input.fq2} {input.viral_fq2} > {output.fq2}'

rule remap_reads:
    input:
        fq1 = rules.add_viral_reads.output.fq1,
        fq2 = rules.add_viral_reads.output.fq2,
    output:
        bam = 'work_snake/bam_remap/{batch}_{phenotype}.bam',
    params:
        bwt_index = refdata.get_ref_file(GENOME, 'bwa', must_exist=False),
        sample = lambda wc: getattr(batch_by_name[wc.batch], wc.phenotype).rgid,
    shell:
        "test -e {params.bwt_index}.bwt &&"
        " bwa mem -R '@RG\\tID:{params.sample}\\tSM:{params.sample}' {params.bwt_index} {input.fq1} {input.fq2} "
        " | samtools sort -Obam -o {output.bam}"

rule extract_hla_contigs:
    input:
        bam = lambda wc: getattr(batch_by_name[wc.batch], wc.phenotype).bam,
    output:
        'work_snake/hla_bam/{batch}_{phenotype}/hla_contigs.txt'
    shell:
        "samtools view -H {input.bam} | "
        "grep 'SN:HLA' | cut -f2 | sed 's/SN://' > {output}"

rule extract_hla_bam:
    input:
        bam = lambda wc: getattr(batch_by_name[wc.batch], wc.phenotype).bam,
        hla_contigs = rules.extract_hla_contigs.output[0],
    output:
        bam = 'work_snake/hla_bam/{batch}_{phenotype}/hla.bam'
    shell:
        'samtools view {input.bam} $(cat {input.hla_contigs}) -Obam -o {output.bam}'

rule subsample_hla_bam:
    input:
        bam = rules.extract_hla_bam.output.bam
    output:
        bam = 'work_snake/hla_bam/{batch}_{phenotype}/hla.downsampled.bam'
    params:
        target_read_cnt = 10_000,
        random_seed = 42,
    run:
        total_cnt = int(subprocess.check_output(f'samtools view -c {input.bam}', shell=True).strip())
        assert total_cnt > 0
        prop = params.target_read_cnt / total_cnt
        # The integer and fractional parts of the -s INT.FRAC option are used separately:
        # the part after the decimal point sets the fraction of templates/pairs to be kept,
        # while the integer part is used as a seed that influences which subset of reads is kept:
        param = str(params.random_seed + prop)
        shell('samtools view {input.bam} -Obam -s {param} -o {output.bam}')

rule merge_with_hla:
    input:
        full_bam = rules.remap_reads.output.bam,
        hla_bam = rules.subsample_hla_bam.output.bam,
    output:
        bam =  'work_snake/bam_remap_with_hla/{batch}_{phenotype}.bam'
    shell:
        'samtools merge {output.bam} {input.full_bam} {input.hla_bam}'

rule index_bam:
    input:
        bam = rules.merge_with_hla.output.bam
    output:
        bai = rules.merge_with_hla.output.bam + '.bai',
    shell:
        'samtools index {input}'

rule populate_downsampled:
    input:
        somatic_vcf    = rules.downsample_somatic.output[0],
        # germline_vcf   = rules.downsample_germline_random100.output[0],
        sv_vcf         = rules.downsample_sv.output[0],
        tumor_bam = 'work_snake/bam_remap_with_hla/{batch}_tumor.bam',
        normal_bam = 'work_snake/bam_remap_with_hla/{batch}_normal.bam',
        tumor_bai = 'work_snake/bam_remap_with_hla/{batch}_tumor.bam.bai',
        normal_bai = 'work_snake/bam_remap_with_hla/{batch}_normal.bam.bai',
        marker = join(PROJECT_NEW_PATH, '.populated_other_files_{batch}.done')  # need the bam_list.csv file to parse DragenProject
    output:
        marker = join(PROJECT_NEW_PATH, '.populated_downsampled_{batch}.done'),
        normal_bam = join(PROJECT_NEW_PATH, '{batch}.bam'),
        normal_bai = join(PROJECT_NEW_PATH, '{batch}.bam.bai'),
        tumor_bam  = join(PROJECT_NEW_PATH, '{batch}_tumor.bam'),
        tumor_bai  = join(PROJECT_NEW_PATH, '{batch}_tumor.bam.bai'),
    params:
        project_copy = PROJECT_NEW_PATH,
    run:
        batch = batch_by_name[wildcards.batch]

        new_project = DragenProject(input_dir=params.project_copy)
        new_batch = new_project.add_batch(batch.name)
        new_batch.add_tumor(batch.tumor.name)
        new_batch.add_normal(batch.normal.name)

        shell(f'bgzip -c {input.somatic_vcf} > {new_batch.somatic_vcf} && tabix {new_batch.somatic_vcf}')
        shell(f'bgzip -c {input.sv_vcf} > {new_batch.sv_vcf} && tabix {new_batch.sv_vcf}')

        shell(f'cp {input.tumor_bam} {new_batch.tumor.bam}')
        shell(f'cp {input.normal_bam} {new_batch.normal.bam}')
        shell(f'cp {input.tumor_bai} {new_batch.tumor.bam}.bai')
        shell(f'cp {input.normal_bai} {new_batch.normal.bam}.bai')

        shell('touch {output.marker}')


rule populate_other_files:
    input:
        qc_files = lambda wc: batch_by_name[wc.batch].all_qc_files(),
        replay_file = lambda wc: batch_by_name[wc.batch].replay_file,
    output:
        marker = join(PROJECT_NEW_PATH, '.populated_other_files_{batch}.done')
    params:
        project_copy = PROJECT_NEW_PATH,
    run:
        files = input.qc_files + [input.replay_file]
        for fpath in files:
            new_path = join(params.project_copy, basename(fpath))
            shell('cp {fpath} {new_path}')
        shell('touch {output.marker}')


######################################
#### Reference data ####
rule prep_gnomad:
    input:
        vcf = refdata.get_ref_file(GENOME, 'gnomad'),
        somatic_roi = rules.project_somatic_roi.output.bed
    output:
        vcf = f'data/genomes/{GENOME}/gnomad_genome.vcf.gz'
    shell:
        'bedtools intersect -header -a {input.vcf} -b {input.somatic_roi} | '
        'bgzip -c > {output.vcf} && tabix {output.vcf}'

rule prep_purple_gc:
    input:
        refdata.get_ref_file(GENOME, 'purple_gc')
    output:
        f'data/genomes/{GENOME}/hmf/GC_profile.1000bp.cnp'
    shell:
        'cp {input} {output}'

rule prep_germline_het:
    input:
        het_vcf = refdata.get_ref_file(GENOME, 'purple_het'),
        roi_bed = rules.project_roi.output[0]
    output:
        f'data/genomes/{GENOME}/hmf/germline_het_pon.vcf.gz'
    shell:
        'bedtools intersect -a {input.het_vcf} -b {input.roi_bed} | bgzip -c > {output} && tabix -p vcf {output}'

rule prep_hotspots:
    input:
        file = refdata.get_ref_file(GENOME, 'hotspots')
    output:
        f'data/genomes/{GENOME}/hmf/KnownHotspots.tsv.gz'
    shell:
        'cp {input} {output}'

rule prep_giab_conf:
    input:
        bed = refdata.get_ref_file(GENOME, 'hmf_giab_conf'),
        roi_bed = rules.project_roi.output[0]
    output:
        f'data/genomes/{GENOME}/hmf/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed.gz'
    shell:
        'bedtools intersect -a {input.bed} -b {input.roi_bed} > {output}'


######################################
#### Panel of normals ####
rule prep_pon_snps_vcf:
    input:
        pon_vcf = join(refdata.get_ref_file(GENOME, 'panel_of_normals_dir'), 'pon.snps.vcf.gz'),
        somatic_roi = rules.project_somatic_roi.output[0]
    output:
        vcf = join(Out_PON_PATH, 'panel_of_normals.snps.vcf.gz')
    shell:
        'bedtools intersect -header -a {input.pon_vcf} -b {input.somatic_roi} | ' 
        'bgzip -c > {output.vcf} && tabix {output.vcf}'

rule prep_pon_indels_vcf:
    input:
        pon_vcf = join(refdata.get_ref_file(GENOME, 'panel_of_normals_dir'), 'pon.indels.vcf.gz'),
        somatic_roi = rules.project_somatic_roi.output[0]
    output:
        vcf = join(Out_PON_PATH, 'panel_of_normals.indels.vcf.gz')
    shell:
        'bedtools intersect -header -a {input.pon_vcf} -b {input.somatic_roi} | ' 
        'bgzip -c > {output.vcf} && tabix {output.vcf}'


# rule prep_pon:
#     input:
#         pon_vcf = rules.prep_pon_vcf.output[0],
#         pon_snakefile = join(Spa_PON_DIR, 'Snakefile.prep_normals')
#     output:
#         pon_vcf = join(Out_PON_PATH, 'panel_of_normals.vcf.gz')
#     params:
#         pon = Out_PON_PATH,
#         vcf_path = lambda wc, input: relpath(input.pon_vcf, Out_PON_PATH)
#     shell:
#         'snakemake -p -s {input.pon_snakefile} --directory {params.pon} --config normal="{params.sname}:{params.vcf_path}"'


# rule prep_pon:
#     input:
#         pon_vcf = rules.prep_pon_vcf.output[0],
#         pon_snakefile = join(Spa_PON_DIR, 'Snakefile.prep_normals')
#     output:
#         pon_vcf = Out_PON_PATH + '/vcfs/' + PON_SAMPLE_NAME + '.vcf.gz',
#         pon_lua = Out_PON_PATH + '/code.lua',
#         pon_toml = Out_PON_PATH + '/annotate_normals_vcfanno.toml'
#     params:
#         sname = PON_SAMPLE_NAME,
#         pon = Out_PON_PATH,
#         vcf_path = lambda wc, input: relpath(input.pon_vcf, Out_PON_PATH)
#     shell:
#         'snakemake -p -s {input.pon_snakefile} --directory {params.pon} --config normal="{params.sname}:{params.vcf_path}"'














