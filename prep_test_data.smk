import os
import yaml
from os.path import join, abspath, relpath, basename, isfile, dirname
from ngs_utils.bcbio import BcbioProject
from ngs_utils.file_utils import splitext_plus, safe_mkdir
from ngs_utils.reference_data import get_key_genes_bed, get_predispose_genes_bed
from ngs_utils import logger as log
from reference_data import api as refdata


shell.executable(os.environ.get('SHELL', 'bash'))
log.init(config.get('debug', False))

GENOME = 'hg38'
BCBIO_DIR = config.get('final', '/g/data/gx8/projects/Saveliev_SEQCII/bcbio/samples/final')
include_names = ['T_SRR7890936_50pc']
run = BcbioProject(BCBIO_DIR, include_samples=include_names)
included_names = [s.name for s in run.samples]
batch_by_name = {b.tumors[0].name: b for b in run.batch_by_name.values() if not b.is_germline()}

bcbio_copy_path = config.get('out', 'data/bcbio_test_project')
tsv_project_path = config.get('out', 'data/tsv_test_project')
bcbio_copy_final_dir = join(bcbio_copy_path, 'final')
bcbio_copy_work_dir = join(bcbio_copy_path, 'work')
bcbio_copy_date = join(bcbio_copy_final_dir, basename(run.date_dir))

mq_list_files = join(run.date_dir, 'multiqc', 'list_files_final.txt')
# if the file exists, we populate QC (because CWL output currently missing this file as longs as all qc)

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
        expand('work_snake/bam_remap/{batch}_tumor.bam{ext}', batch=batch_by_name.keys(), ext=['', '.bai']),
        expand('work_snake/bam_remap/{batch}_normal.bam{ext}', batch=batch_by_name.keys(), ext=['', '.bai']),
        'work_snake/roi.bed',
        expand(join(bcbio_copy_path, '.populated_batch_{batch}.done'), batch=batch_by_name.keys()),
        join(bcbio_copy_path, '.populated_qc.done') if isfile(mq_list_files) else [],
        join(bcbio_copy_date, '.rsynced.done'),
        join(tsv_project_path, 'input.tsv'),
        pon_snps_vcf   = join(Out_PON_PATH, 'panel_of_normals.snps.vcf.gz'),
        pon_indels_vcf = join(Out_PON_PATH, 'panel_of_normals.indels.vcf.gz'),
        # gnomad          = f'data/genomes/{GENOME}/gnomad_genome.vcf.gz',
        # purple_gc       = f'data/genomes/{GENOME}/hmf/GC_profile.1000bp.cnp',
        # purple_het      = f'data/genomes/{GENOME}/hmf/GermlineHetPon.hg19.vcf.gz',
        # hotspots        = f'data/genomes/{GENOME}/hotspots/merged.vcf.gz',
        # hmf_giab_conf   = f'data/genomes/{GENOME}/hmf/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed.gz'



######################################
#### SOMATIC ####
rule downsample_somatic:
    input:
        lambda wc: join(run.date_dir, f'{batch_by_name[wc.batch].name}-ensemble-annotated.vcf.gz')
    output:
        'work_snake/somatic/{batch}-ensemble.vcf'
    shell:
        extract_genes_from_vcf()

rule somatic_roi:
    input:
        rules.downsample_somatic.output[0]
    output:
        'work_snake/somatic/{batch}-ensemble.bed'
    shell:
        vcf_to_bed()

######################################
#### GERMLINE ####
rule peddy_bed:
    input:
        peddy_sites = f'data/peddy_GRCH38.sites',
    output:
        peddy_bed = 'work_snake/germline/peddy.bed'
    shell:
        """cat {input} | tr ':' '\t' | awk '{{ print "chr" $1 "\t" $2-1 "\t" $2 }}' > {output}"""

rule germline_bed:
    input:
        peddy_bed = 'work_snake/germline/peddy.bed',
        predispose_genes_bed = get_predispose_genes_bed(GENOME, coding_only=False),
    output:
        bed = 'work_snake/germline/target.bed',
    shell:
        'cat <(cat {input.predispose_genes_bed} | sort -R {SORT_SEED} | head -n10 | cut -f1-3)'
        ' {input.peddy_bed} | bedtools sort -i stdin | bedtools merge -i stdin > {output.bed}'

rule downsample_germline:
    input:
        vcf = lambda wc: join(run.date_dir, f'{batch_by_name[wc.batch].normal.name}-germline-ensemble-annotated.vcf.gz'),
        bed = rules.germline_bed.output
    output:
        vcf = 'work_snake/germline/{batch}-germline-ensemble-predispose-genes.vcf'
    shell:
        'bcftools view -T {input.bed} {input.vcf} > {output.vcf}'

# rule downsample_germline_random100:
#     input:
#         rules.downsample_germline.output[0]
#     output:
#         'work_snake/germline/{batch}-germline-ensemble-100.vcf'
#     shell:
#         downsample_vcf(100)

rule germline_roi:
    input:
        rules.downsample_germline.output[0]
    output:
        'work_snake/germline/{batch}-germline-ensemble.bed'
    shell:
        vcf_to_bed()


######################################
#### MANTA ####
rule downsample_manta:
    input:
        lambda wc: join(batch_by_name[wc.batch].tumors[0].dirpath, f'{batch_by_name[wc.batch].name}-manta.vcf.gz')
    output:
        'work_snake/sv/{batch}-manta.vcf'
    shell:
        'bcftools view -f.,PASS,REJECT,Intergenic,MissingAnn {input} -o {output}'

rule manta_roi:
    input:
        rules.downsample_manta.output[0]
    output:
        'work_snake/sv/{batch}-manta.bed'
    shell:
        vcf_to_bed()

######################################
#### BAMS ####
rule batch_roi:
    input:
        rules.somatic_roi.output[0],
        rules.germline_roi.output[0]
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
        bed = 'data/conpair_markers/GRCh37.bed' if GENOME == 'GRCh37' else \
        'data/conpair_markers/hg38.liftover.bed',
        fai = refdata.get_ref_file(GENOME, key='fa') + '.fai'
    output:
        'work_snake/conpair_roi.bed'
    shell:
        'head -n100 {input.bed} '
        ' | bedtools slop -b 1000 -i stdin -g {input.fai}'
        '> {output}'

rule project_roi:
    input:
        beds = expand(rules.batch_roi.output[0], batch=batch_by_name.keys()) +
               [rules.conpair_roi.output[0]],
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
        beds = expand(rules.somatic_roi.output[0], batch=batch_by_name.keys()),
        fai = refdata.get_ref_file(GENOME, key='fa') + '.fai',
    output:
        'work_snake/somatic_roi.bed'
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
        bam_namesorted = 'work_snake/bam_subset/{batch}_{phenotype}.bam'
    shell:
        'sambamba slice {input.bam} -L {input.roi_bed}'
        ' | samtools sort -n -Obam'
        ' -o {output} /dev/stdin'

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

rule index_bam:
    input:
        bam = rules.remap_reads.output.bam
    output:
        bai = rules.remap_reads.output.bam + '.bai',
    shell:
        'samtools index {input}'

# rule bcbio_csv:
#     input:
#         expand('work_snake/{batch}_tumor.bam{ext}', batch=batch_by_name.keys(), ext=['', '.bai']),
#         expand('work_snake/{batch}_normal.bam{ext}', batch=batch_by_name.keys(), ext=['', '.bai']),
#         roi_bed = rules.project_roi.output[0]
#     output:
#         'work_snake/bcbio.csv'
#     run:
#         sample_dicts = []
#         samples_added = set()
#         for s in run.samples:
#             if s.phenotype == 'germline':
#                 pass
#             elif s.name in samples_added:
#                 pass
#             else:
#                 samples_added.add(s.name)
#                 sample_dicts.append({
#                     'samplename': s.batch.tumor.name + '_' + s.phenotype,
#                     'description': s.name,
#                     'batch': s.batch.name,
#                     'phenotype': s.phenotype,
#                     'variant_regions': abspath(input.roi_bed),
#                 })
#         with open(output[0], 'w') as out:
#             writer = csv.DictWriter(out, fieldnames=sample_dicts[0].keys())
#             writer.writeheader()
#             for sd in sample_dicts:
#                 writer.writerow(sd)
#
# rule make_bcbio_yaml:
#     input:
#         bams = expand('work_snake/{batch}_{phenotype}.bam', batch=batch_by_name.keys(), phenotype=['tumor', 'normal']),
#         csv = rules.bcbio_csv.output[0],
#         template = 'bcbio_template.yaml'
#     output:
#         'bcbio/config/bcbio.csv'
#     shell:
#         '/home/vlad/bcbio/anaconda/bin/bcbio_nextgen.py -w template {input.template} {input.csv} {input.bams}'
#
# rule run_bcbio:
#     input:
#         rules.make_bcbio_yaml.output[0]
#     output:
#         'bcbio/final'
#     shell:
#         'cp run.sh bcbio/work; cd bcbio/work; sbatch run.sh'


rule rsync_project:
    input:
        run.date_dir,
        run.config_dir
    output:
        marker = join(bcbio_copy_date, '.rsynced.done')
    params:
        project_copy = bcbio_copy_path,
    run:
        exclude_samples = ''
        safe_mkdir(join(params.project_copy, 'config'))
        exclude_samples = [s.name for s in BcbioProject(BCBIO_DIR, include_samples=included_names).samples]
        exclude_samples = ' '.join(f'--exclude {sn}' for sn in exclude_samples)
        new_yaml = join(params.project_copy, 'config', basename(run.bcbio_yaml_fpath))
        with open(run.bcbio_yaml_fpath) as inp, open(new_yaml, 'w') as out:
            data = yaml.load(inp)
            data['details'] = [d for d in data['details'] if d['description'] in included_names]
            yaml.dump(data, out)

        shell(
            f'rsync -tavz'
            f' --exclude ".snakemake/"'
            f' --exclude "umccrised/"'
            f' --exclude "test/"'
            f' {exclude_samples}'
            f' --include "*/"'
            f' --include "data_versions.csv"'
            f' --include "programs.txt"'
            f' --include "multiqc_report.html"'
            f' --exclude "*"'
            f' {run.final_dir} {params.project_copy}'
        )
        shell('touch {output.marker}')


if isfile(mq_list_files):
    qc_files = [line.strip() for line in open(mq_list_files)
                if line.strip() and line.strip() != 'trimmed' and any(sn in line for sn in included_names)]
    inp_qc_files = [join(run.final_dir, qc_f) for qc_f in qc_files]
    out_qc_files = [join(bcbio_copy_final_dir, qc_f) for qc_f in qc_files]

    rule populate_qc:
        input:
            mq_list_files = join(run.date_dir, 'multiqc', 'list_files_final.txt'),
            mq_yaml = join(run.date_dir, 'multiqc', 'multiqc_config.yaml'),
            inp_qc_files = inp_qc_files
        output:
            mq_list_files = join(bcbio_copy_date, 'multiqc', 'list_files_final.txt'),
            mq_yaml= join(bcbio_copy_date, 'multiqc', 'multiqc_config.yaml'),
            out_qc_files = out_qc_files,
            marker = bcbio_copy_path + '/.populated_qc.done'
        run:
            shell(f'cp {input.mq_list_files} {output.mq_list_files}')
            shell(f'cp {input.mq_yaml} {output.mq_yaml}')
            for i_f, o_f in zip(input.inp_qc_files, output.out_qc_files):
                shell(f'cp {i_f} {o_f}')
            with open(output.mq_list_files, 'w') as out:
                for qf in qc_files:
                    out.write(qf + '\n')
            shell('touch {output.marker}')


rule populate_batch:
    input:
        somatic_vcf = rules.downsample_somatic.output[0],
        germline_vcf = rules.downsample_germline.output[0],
        manta_vcf = rules.downsample_manta.output[0],
        tumor_bam = 'work_snake/bam_remap/{batch}_tumor.bam',
        normal_bam = 'work_snake/bam_remap/{batch}_normal.bam',
        tumor_bai = 'work_snake/bam_remap/{batch}_tumor.bam.bai',
        normal_bai = 'work_snake/bam_remap/{batch}_normal.bam.bai',
        synced = join(bcbio_copy_date, '.rsynced.done'),
    output:
        marker = bcbio_copy_path + '/.populated_batch_{batch}.done'
    run:
        batch = batch_by_name[wildcards.batch]

        somatic_vcf = join(bcbio_copy_path, 'final', basename(run.date_dir), f'{batch.name}-ensemble-annotated.vcf.gz')
        shell(f'bcftools sort {input.somatic_vcf} | bgzip -c > {somatic_vcf} && tabix {somatic_vcf}')

        germline_vcf = join(bcbio_copy_path, 'final',
                            basename(run.date_dir),
                            f'{batch.normals[0].name}-germline-ensemble-annotated.vcf.gz')
        shell(f'bcftools sort {input.germline_vcf} | bgzip -c > {germline_vcf} && tabix {germline_vcf}')

        manta_vcf = join(bcbio_copy_path, 'final',
                         basename(batch.tumors[0].dirpath),
                         f'{batch.name}-manta.vcf.gz')
        shell(f'bcftools sort {input.manta_vcf} | bgzip -c > {manta_vcf} && tabix {manta_vcf}')

        tumor_bam_name = basename(batch.tumors[0].bam)
        normal_bam_name = basename(batch.normals[0].bam)
        tumor_dir = join(bcbio_copy_path, 'final', basename(batch.tumors[0].dirpath))
        normal_dir = join(bcbio_copy_path, 'final', basename(batch.normals[0].dirpath))
        shell(f'cp {input.tumor_bam} {tumor_dir}/{tumor_bam_name}')
        shell(f'cp {input.tumor_bai} {tumor_dir}/{tumor_bam_name}.bai')
        shell(f'cp {input.normal_bam} {normal_dir}/{normal_bam_name}')
        shell(f'cp {input.normal_bai} {normal_dir}/{normal_bam_name}.bai')

        shell('touch {output.marker}')


######################################
#### TSV test project ####
rule populate_tsv_project_batch:
    input:
        tumor_bam =  'work_snake/bam_remap/{batch}_tumor.bam',
        normal_bam = 'work_snake/bam_remap/{batch}_normal.bam',
        tumor_bai =  'work_snake/bam_remap/{batch}_tumor.bam.bai',
        normal_bai = 'work_snake/bam_remap/{batch}_normal.bam.bai',
    output:
        marker = tsv_project_path + '/.populated_batch_{batch}.done'
    run:
        batch = batch_by_name[wildcards.batch]
        tumor_bam_name = basename(batch.tumors[0].bam)
        normal_bam_name = basename(batch.normals[0].bam)
        shell(f'cp {input.tumor_bam} {tsv_project_path}/{tumor_bam_name}')
        shell(f'cp {input.tumor_bai} {tsv_project_path}/{tumor_bam_name}.bai')
        shell(f'cp {input.normal_bam} {tsv_project_path}/{normal_bam_name}')
        shell(f'cp {input.normal_bai} {tsv_project_path}/{normal_bam_name}.bai')
        shell('touch {output.marker}')

rule populate_tsv_projejct:
    input:
        expand(join(tsv_project_path, '.populated_batch_{batch}.done'), batch=batch_by_name.keys()),
    output:
        input_tsv = join(tsv_project_path, 'input.tsv')
    run:
        with open(output.input_tsv, 'w') as f:
            f.write('\t'.join(['sample', 'wgs', 'normal', 'exome', 'exome_normal',
                               'rna', 'rna_bcbio', 'rna_sample']) + '\n')
            for batch in batch_by_name.values():
                tumor_bam_name = basename(batch.tumors[0].bam)
                normal_bam_name = basename(batch.normals[0].bam)
                fields = [
                    f'{tumor_bam_name}',
                    f'{normal_bam_name}',
                    '.',
                    '.',
                ]
                # if patients.get(apgi_id, f'Reprocessed_RNA'):
                #     rna_bam = patients.find_rna_bams(apgi_id)[0]
                #     fields.extend([
                #         rna_bam,
                #         join(Patients.reprocessed_base_dir_raijin,
                #              patients.get(apgi_id, f'Reprocessed_RNA')),
                #         basename(rna_bam).replace('-ready.bam', ''),
                #     ])
                # else:
                fields.extend(['.', '.', '.'])
                f.write('\t'.join([batch.name] + fields) + '\n')

######################################
#### Reference data ####
rule prep_gnomad:
    input:
        vcf = refdata.get_ref_file(GENOME, 'gnomad'),
        somatic_roi = rules.project_somatic_roi.output[0]
    output:
        vcf = f'data/genomes/{GENOME}/gnomad_genome.vcf.gz'
    shell:
        'bedtools intersect -header -a {input.vcf} -b {input.somatic_roi} | ' \
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
        vcf = refdata.get_ref_file(GENOME, 'purple_het'),
        roi_bed = rules.project_roi.output[0]
    output:
        f'data/genomes/{GENOME}/hmf/GermlineHetPon.hg19.vcf.gz'
    shell:
        'bedtools intersect -a {input.vcf} -b {input.roi_bed} | '
        'bgzip -c > {output} && tabix -p vcf {output}'

rule prep_hotspots:
    input:
        file = refdata.get_ref_file(GENOME, 'hotspots')
    output:
        f'data/genomes/{GENOME}/hotspots/merged.vcf.gz'
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
        pon_vcf = join(refdata.get_ref_file(GENOME, 'panel_of_normals_dir'), 'panel_of_normals.snps.vcf.gz'),
        somatic_roi = rules.project_somatic_roi.output[0]
    output:
        vcf = join(Out_PON_PATH, 'panel_of_normals.snps.vcf.gz')
    shell:
        'bedtools intersect -header -a {input.pon_vcf} -b {input.somatic_roi} | ' \
        'bgzip -c > {output.vcf} && tabix {output.vcf}'

rule prep_pon_indels_vcf:
    input:
        pon_vcf = join(refdata.get_ref_file(GENOME, 'panel_of_normals_dir'), 'panel_of_normals.indels.vcf.gz'),
        somatic_roi = rules.project_somatic_roi.output[0]
    output:
        vcf = join(Out_PON_PATH, 'panel_of_normals.indels.vcf.gz')
    shell:
        'bedtools intersect -header -a {input.pon_vcf} -b {input.somatic_roi} | ' \
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














