import os
import sys
import traceback
from os.path import dirname, join, exists, isfile, splitext, basename, isdir, relpath, getctime, getsize, abspath, expanduser
import subprocess
from nose.plugins.attrib import attr

try:
    from ngs_utils.testing import BaseTestCase, echo, check_call, vcf_ignore_lines, swap_output
    from ngs_utils.utils import is_az, is_local, is_travis
    from ngs_utils.file_utils import safe_mkdir
    from ngs_utils.call_process import run_simple
    from hpc_utils.hpc import find_loc
except ImportError as e:
    traceback.print_exc()
    sys.stderr.write('\nUmccrise is not installed properly. Refer to the README.md for installation\n')
    sys.exit(1)

loc = find_loc()


TUMORS = ['cup_tissue']
BATCHES = ['cup']
PROJECT = 'cup_sc932'


class Test_umccrise(BaseTestCase):
    script = 'umccrise'

    test_data_clone = join(dirname(__file__), '.')
    data_dir = join(test_data_clone, BaseTestCase.data_dir)
    results_dir = join(test_data_clone, BaseTestCase.results_dir)
    gold_standard_dir = join(test_data_clone, BaseTestCase.gold_standard_dir)
    try:
        from hpc_utils import hpc
        loc = hpc.find_loc()
    except:  # When testing the dockerized version
        loc = None

    reuse = False  # Run on top of existing latest results. Also controlled with TEST_REUSE
    only_diff = False  # Do not run, just diff the latest results against the gold standard. Also controlled with TEST_ONLY_DIFF

    def setUp(self):
        assert os.system(f'which {self.script}') == 0, 'Umccrise is not installed. Refer to the README.md for installation'

        if not Test_umccrise.loc or Test_umccrise.loc.name == 'travis':
            echo('Server is not recognized, downloaded the reference data')
            ref_fasta_path = join(Test_umccrise.test_data_clone, 'data/genomes/Hsapiens/GRCh37/seq/GRCh37.fa')
            if not isfile(ref_fasta_path):
                print('Downloading GRCh37 genome...')
                run_simple(f'''wget -nv --no-check-certificate -c https://s3.amazonaws.com/biodata/genomes/GRCh37-seq.tar.gz && 
tar -xzvpf GRCh37-seq.tar.gz --directory {Test_umccrise.test_data_clone}/data/genomes/Hsapiens/GRCh37 && 
rm -f GRCh37-seq.tar.gz && 
gunzip {ref_fasta_path}.gz''')

        BaseTestCase.setUp(self)

    def _run_umccrise(self, bcbio_dirname, parallel=False, docker_wrapper_mode=False, pcgr=False):
        results_dir = join(self.results_dir, bcbio_dirname)
        bcbio_dir = join(self.data_dir, bcbio_dirname)
        cmdl = f'{self.script} {bcbio_dir} -o {results_dir} --no-s3'
        if not Test_umccrise.loc or docker_wrapper_mode:
            cmdl += f' --bcbio-genomes {Test_umccrise.test_data_clone}/data/genomes'
            cmdl += f' --pon {Test_umccrise.test_data_clone}/data/panel_of_normals'
        if parallel:
            cmdl += ' -j 10'
        if docker_wrapper_mode:
            cmdl += ' --docker'
        if pcgr and docker_wrapper_mode:
            cmdl += f' --pcgr-data {Test_umccrise.loc.pcgr_dir}/data'
        self._run_cmd(cmdl, bcbio_dir, results_dir)
        return results_dir

    def _check_file(self, diff_failed, path, ignore_matching_lines=None, wrapper=None, check_diff=True):
        try:
            self._check_file_throws(path, ignore_matching_lines=ignore_matching_lines, wrapper=wrapper, check_diff=check_diff)
        except subprocess.CalledProcessError as e:
            sys.stderr.write(f'Error: {e}\n')
            return True
        except AssertionError as e:
            sys.stderr.write(f'Error: {e}\n')
            return True
        return diff_failed

    @attr('normal')
    def test(self, docker_wrapper_mode=False, pcgr=False):
        results_dir = self._run_umccrise(bcbio_dirname='bcbio_test_project', parallel=False,
                                         docker_wrapper_mode=docker_wrapper_mode, pcgr=pcgr)

        failed = False
        failed = self._check_file(failed, f'{results_dir}/log/{PROJECT}-config/{PROJECT}-template.yaml'                                                 )
        failed = self._check_file(failed, f'{results_dir}/log/{PROJECT}-config/{PROJECT}.csv'                                                           )
        failed = self._check_file(failed, f'{results_dir}/log/{PROJECT}-config/{PROJECT}.yaml'                                                          )
        failed = self._check_file(failed, f'{results_dir}/log/{PROJECT}-data_versions.csv'                                                              )
        failed = self._check_file(failed, f'{results_dir}/log/{PROJECT}-programs.txt'                                                                   )
        failed = self._check_file(failed, f'{results_dir}/{PROJECT}-multiqc_report.html'                                              , check_diff=False)
        for T, B in zip(TUMORS, BATCHES):
            key = f'{B}__{T}'
            failed = self._check_file(failed, f'{results_dir}/{key}/coverage/{key}-indexcov/index.html'                               , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/coverage/{key}-normal.callable.bed'                                                 )
            failed = self._check_file(failed, f'{results_dir}/{key}/coverage/{key}-normal.depth.bed'                                                    )
            failed = self._check_file(failed, f'{results_dir}/{key}/coverage/{key}-tumor.depth.bed'                                                     )
            failed = self._check_file(failed, f'{results_dir}/{key}/igv/{key}-roi.bed'                                                , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/igv/{key}-normal.mini.bam'                                        , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/igv/{key}-normal.mini.bam.bai'                                    , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/igv/{key}-tumor.mini.bam'                                         , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/igv/{key}-tumor.mini.bam.bai'                                     , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/pcgr/input/{key}-normal.vcf.gz'                                   , vcf_ignore_lines, check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/pcgr/input/{key}-normal.vcf.gz.tbi'                               , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/pcgr/input/{key}-somatic.vcf.gz'                                  , vcf_ignore_lines, check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/pcgr/input/{key}-somatic.vcf.gz.tbi'                              , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/pcgr/input/{key}-somatic-cna.tsv'                                                   )
            if loc.name == 'spartan':
                failed = self._check_file(failed, f'{results_dir}/{key}/pcgr/{key}-somatic.pcgr_acmg.html'                            , check_diff=False)
                failed = self._check_file(failed, f'{results_dir}/{key}/pcgr/{key}-normal.pcgr_acmg.html'                             , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/{key}-rmd_report.html'                                            , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/work/{key}/rmd/afs/af_tumor.txt'                                        , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/work/{key}/rmd/afs/af_tumor_az300.txt'                                  , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/work/{key}/rmd/ensemble-with_chr_prefix.vcf'                            , vcf_ignore_lines, check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/small_variants/{key}-normal-ensemble-cancer_genes.vcf.gz'         , vcf_ignore_lines, check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/small_variants/{key}-normal-ensemble-cancer_genes.vcf.gz.tbi'     , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/small_variants/{key}-somatic-ensemble-pon_softfiltered.vcf.gz'    , vcf_ignore_lines, check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/small_variants/{key}-somatic-ensemble-pon_softfiltered.vcf.gz.tbi', check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/small_variants/{key}-somatic-ensemble-pon_hardfiltered.vcf.gz'    , vcf_ignore_lines, check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/small_variants/{key}-somatic-ensemble-pon_hardfiltered.vcf.gz.tbi', check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/work/{key}/structural/{key}-cnvkit-nolabels.cns'                                          )
            failed = self._check_file(failed, f'{results_dir}/{key}/structural/{key}-cnvkit-diagram.pdf'                              , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/structural/{key}-sv-prioritize-manta.bedpe'                                         )
            failed = self._check_file(failed, f'{results_dir}/{key}/structural/{key}-sv-prioritize-manta.ribbon.bed'                                    )
            failed = self._check_file(failed, f'{results_dir}/{key}/structural/{key}-sv-prioritize-manta-pass.tsv'                                      )
            failed = self._check_file(failed, f'{results_dir}/{key}/structural/{key}-sv-prioritize-manta.vcf'                         , vcf_ignore_lines)
            failed = self._check_file(failed, f'{results_dir}/{key}/purple/{key}.purple.cnv'                                          , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/purple/{key}.purple.circos.png'                                   , check_diff=False)

        assert not failed, 'some of file checks have failed'

    @attr('docker')
    def test_docker(self):
        self.test(docker=True)

    @attr('docker_with_pcgr')
    def test_docker(self):
        self.test(docker=True, pcgr=True)

