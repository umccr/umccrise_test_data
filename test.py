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
            ref_fasta_dir = join(Test_umccrise.test_data_clone, 'data/genomes/GRCh37')
            ref_fasta_path = join(ref_fasta_dir, 'GRCh37.fa')
            if not isfile(ref_fasta_path):
                print('Downloading GRCh37 genome...')
                run_simple(f'''wget -nv --no-check-certificate -c https://s3.amazonaws.com/biodata/genomes/GRCh37-seq.tar.gz && 
tar -xzvpf GRCh37-seq.tar.gz --directory {ref_fasta_dir} && 
mv {ref_fasta_dir}/seq/* {ref_fasta_dir}/ && 
rm -f GRCh37-seq.tar.gz && 
gunzip {ref_fasta_path}.gz''')

        BaseTestCase.setUp(self)

    def _run_umccrise(self, bcbio_dirname, parallel=False, docker_wrapper_mode=False, skip_pcgr=False):
        results_dir = join(self.results_dir, bcbio_dirname)
        bcbio_dir = join(self.data_dir, bcbio_dirname)
        cmdl = f'{self.script} {bcbio_dir} -o {results_dir} --no-s3'
        if not Test_umccrise.loc or docker_wrapper_mode:
            cmdl += f' --genomes {Test_umccrise.test_data_clone}/data/genomes'
        if parallel:
            cmdl += ' -j10'
        if docker_wrapper_mode:
            cmdl += ' --docker'
        if not skip_pcgr and docker_wrapper_mode:
            cmdl += f' --pcgr-data {Test_umccrise.loc.pcgr_data}'
        if skip_pcgr:
            cmdl += ' --no-pcgr'
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
    def test(self, docker_wrapper_mode=False, skip_pcgr=False):
        results_dir = self._run_umccrise(bcbio_dirname='bcbio_test_project', parallel=False,
                                         docker_wrapper_mode=docker_wrapper_mode, skip_pcgr=skip_pcgr)

        failed = False
        failed = self._check_file(failed, f'{results_dir}/log/{PROJECT}-config/{PROJECT}-template.yaml', check_diff=False)
        failed = self._check_file(failed, f'{results_dir}/log/{PROJECT}-config/{PROJECT}.csv', check_diff=False)
        failed = self._check_file(failed, f'{results_dir}/log/{PROJECT}-config/{PROJECT}.yaml', check_diff=False)
        failed = self._check_file(failed, f'{results_dir}/log/{PROJECT}-data_versions.csv', check_diff=False)
        failed = self._check_file(failed, f'{results_dir}/log/{PROJECT}-programs.txt', check_diff=False)
        for T, B in zip(TUMORS, BATCHES):
            batch = f'{B}__{T}'
            failed = self._check_file(failed, f'{results_dir}/{batch}/{batch}-multiqc_report.html', check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{batch}/coverage/{batch}-indexcov/index.html', check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{batch}/igv/{batch}-roi.bed', check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{batch}/pcgr/{batch}-somatic.pcgr.html', check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{batch}/pcgr/{batch}-normal.cpsr.html', check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{batch}/{batch}_book/purple-results.html', check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{batch}/{batch}_book/index.html', check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{batch}/structural/{batch}-sv-prioritize-manta.bedpe')
            failed = self._check_file(failed, f'{results_dir}/{batch}/structural/{batch}-sv-prioritize-manta.ribbon.bed')
            failed = self._check_file(failed, f'{results_dir}/{batch}/structural/{batch}-sv-prioritize-manta-pass.tsv')
            failed = self._check_file(failed, f'{results_dir}/{batch}/structural/{batch}-sv-prioritize-manta.vcf', vcf_ignore_lines)
            failed = self._check_file(failed, f'{results_dir}/{batch}/purple/{batch}.purple.cnv', wrapper='cut -f1,2,3,8,9,10')
            failed = self._check_file(failed, f'{results_dir}/{batch}/purple/{batch}.purple.circos.png', check_diff=False)

        assert not failed, 'some of file checks have failed'

    # @attr('skip_pcgr')
    # def test_no_pcgr(self):
    #     self.test(skip_pcgr=True)

    # @attr('docker')
    # def test_docker(self):
    #     self.test(docker=True, pcgr=False)
    #
    # @attr('docker_with_pcgr')
    # def test_docker(self):
    #     self.test(docker=True, pcgr=True)

