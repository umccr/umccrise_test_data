import os
from ngs_utils.dragen import DragenProject
from os.path import join, abspath, relpath, basename, isfile, dirname
from ngs_utils.reference_data import get_key_genes_bed
from hpc_utils import hpc
from ngs_utils import logger as log

shell.executable(os.environ.get('SHELL', 'bash'))
shell.prefix('')  # Fixes Snakemake on Spartan

log.init(config.get('debug', False))
print('debug = ', log.is_debug)

GENOME = 'hg38'
#DRAGEN_DIR = config.get('dragen_dir', '/g/data/gx8/extras/umccrise_017_2020_Jan/umccrise_test_data/running_dragen/T_SRR7890902_20pc')
DRAGEN_DIR = config.get('dragen_dir', '/g/data/gx8/extras/umccrise_017_2020_Jan/umccrise_test_data/running_dragen/P025_with_sv')
run = DragenProject(DRAGEN_DIR)
included_names = [s.name for s in run.samples]
batch_by_name = {b.tumor.name: b for b in run.batch_by_name.values() if not b.is_germline()}


SORT_SEED = "--random-source <(echo 42)"


Out_PON_PATH = f'data/genomes/{GENOME}/panel_of_normals'
PROJECT_NEW_PATH = config.get('out', 'data/dragen_test_project')


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
        expand('work_snake/{batch}_tumor.bam{ext}' , batch=batch_by_name.keys(), ext=['', '.bai']),
        expand('work_snake/{batch}_normal.bam{ext}', batch=batch_by_name.keys(), ext=['', '.bai']),
        'work_snake/roi.bed',
        expand(join(PROJECT_NEW_PATH, '.populated_downsampled_{batch}.done'), batch=batch_by_name.keys()),
        expand(join(PROJECT_NEW_PATH, '.populated_other_files_{batch}.done')   , batch=batch_by_name.keys()),
        # reference data:
        pon_snps_vcf    = join(Out_PON_PATH, 'panel_of_normals.snps.vcf.gz'),
        pon_indels_vcf  = join(Out_PON_PATH, 'panel_of_normals.indels.vcf.gz'),
        gnomad          = f'data/genomes/{GENOME}/gnomad_genome.vcf.gz',
        purple_gc       = f'data/genomes/{GENOME}/hmf/GC_profile.1000bp.cnp',
        purple_het      = f'data/genomes/{GENOME}/hmf/germline_het_pon.vcf.gz',
        hmf_hotspot     = f'data/genomes/{GENOME}/hmf/KnownHotspots.tsv.gz',
        hmf_giab_conf   = f'data/genomes/{GENOME}/hmf/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed.gz'


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
        'work_snake/somatic/{batch}.bed'
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
        'work_snake/sv/{batch}.bed'
    shell:
        vcf_to_bed()

######################################
#### BAMS ####
rule batch_roi:
    input:
        rules.somatic_roi.output[0],
        # rules.germline_roi.output[0]
    output:
        'work_snake/{batch}-roi.bed'
    shell:
        'cat {input} | bedtools sort -i stdin | bedtools merge -i stdin > {output}'

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
        'data/conpair_markers/GRCh37.bed' if GENOME == 'GRCh37' else 'data/conpair_markers/hg38.liftover.bed'
    output:
        'work_snake/conpair_roi.bed'
    shell:
        'head -n100 {input} > {output}'

rule project_roi:
    input:
        expand(rules.batch_roi.output[0], batch=batch_by_name.keys()),
        rules.conpair_roi.output[0]
    output:
        'work_snake/roi.bed'
    shell:
        'cat {input} | bedtools sort -i stdin | bedtools merge -i stdin > {output}'

rule project_somatic_roi:
    input:
        expand(rules.somatic_roi.output[0], batch=batch_by_name.keys())
    output:
        'work_snake/somatic_roi.bed'
    shell:
        'cat {input} | bedtools sort -i stdin | bedtools merge -i stdin > {output}'

rule subset_bam:
    input:
        bam = lambda wc: getattr(batch_by_name[wc.batch], wc.phenotype).bam,
        roi_bed = rules.project_roi.output[0]
    output:
        'work_snake/{batch}_{phenotype}.bam'
    shell:
        'sambamba slice {input.bam} -L {input.roi_bed} | samtools sort -Obam -o {output} /dev/stdin'

rule index_bam:
    input:
        'work_snake/{batch}_{phenotype}.bam'
    output:
        'work_snake/{batch}_{phenotype}.bam.bai'
    shell:
        'samtools index {input}'


rule populate_downsampled:
    input:
        somatic_vcf    = rules.downsample_somatic.output[0],
        # germline_vcf   = rules.downsample_germline_random100.output[0],
        sv_vcf         = rules.downsample_sv.output[0],
        tumor_bam      = 'work_snake/{batch}_tumor.bam',
        normal_bam     = 'work_snake/{batch}_normal.bam',
        tumor_bai      = 'work_snake/{batch}_tumor.bam.bai',
        normal_bai     = 'work_snake/{batch}_normal.bam.bai',
    output:
        marker = join(PROJECT_NEW_PATH, '.populated_downsampled_{batch}.done')
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
        vcf = hpc.get_ref_file(GENOME, 'gnomad'),
        somatic_roi = rules.project_somatic_roi.output[0]
    output:
        vcf = f'data/genomes/{GENOME}/gnomad_genome.vcf.gz'
    shell:
        'bedtools intersect -header -a {input.vcf} -b {input.somatic_roi} | '
        'bgzip -c > {output.vcf} && tabix {output.vcf}'

rule prep_purple_gc:
    input:
        hpc.get_ref_file(GENOME, 'purple_gc')
    output:
        f'data/genomes/{GENOME}/hmf/GC_profile.1000bp.cnp'
    shell:
        'cp {input} {output}'

rule prep_germline_het:
    input:
        het_vcf = hpc.get_ref_file(GENOME, 'purple_het'),
        roi_bed = rules.project_roi.output[0]
    output:
        f'data/genomes/{GENOME}/hmf/germline_het_pon.vcf.gz'
    shell:
        'bedtools intersect -a {input.het_vcf} -b {input.roi_bed} | bgzip -c > {output} && tabix -p vcf {output}'

rule prep_hotspots:
    input:
        file = hpc.get_ref_file(GENOME, 'hotspots')
    output:
        f'data/genomes/{GENOME}/hmf/KnownHotspots.tsv.gz'
    shell:
        'cp {input} {output}'

rule prep_giab_conf:
    input:
        bed = hpc.get_ref_file(GENOME, 'hmf_giab_conf'),
        roi_bed = rules.project_roi.output[0]
    output:
        f'data/genomes/{GENOME}/hmf/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed.gz'
    shell:
        'bedtools intersect -a {input.bed} -b {input.roi_bed} > {output}'


######################################
#### Panel of normals ####
rule prep_pon_snps_vcf:
    input:
        pon_vcf = join(hpc.get_ref_file(GENOME, 'panel_of_normals_dir'), 'pon.snps.vcf.gz'),
        somatic_roi = rules.project_somatic_roi.output[0]
    output:
        vcf = join(Out_PON_PATH, 'panel_of_normals.snps.vcf.gz')
    shell:
        'bedtools intersect -header -a {input.pon_vcf} -b {input.somatic_roi} | ' 
        'bgzip -c > {output.vcf} && tabix {output.vcf}'

rule prep_pon_indels_vcf:
    input:
        pon_vcf = join(hpc.get_ref_file(GENOME, 'panel_of_normals_dir'), 'pon.indels.vcf.gz'),
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














