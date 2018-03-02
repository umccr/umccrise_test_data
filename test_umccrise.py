import traceback
import os
import sys
from os.path import dirname, join, exists, isfile, splitext, basename, isdir, relpath, getctime, getsize, abspath, expanduser
from datetime import datetime
import shutil
import subprocess

from ngs_utils.testing import BaseTestCase, info, check_call
from ngs_utils.utils import is_az, is_local, is_travis
from ngs_utils.file_utils import safe_mkdir


# Run on top of existing latest results
REUSE = False
# Do not run, just diff the latest results against the gold standard
ONLY_DIFF = True


BATCHES = ['cup_tissue']
PROJECT = 'cup_sc932'


class Test_umccrise(BaseTestCase):
    script = 'umccrise'

    data_dir = join(dirname(__file__), BaseTestCase.data_dir)
    results_dir = join(dirname(__file__), BaseTestCase.results_dir)
    gold_standard_dir = join(dirname(__file__), BaseTestCase.gold_standard_dir)

    def setUp(self):
        BaseTestCase.setUp(self)

    def _run_umccrise(self, bcbio_dirname, parallel=False):
        results_dir = join(self.results_dir, bcbio_dirname)
        ran_with_error = False

        if not ONLY_DIFF:
            bcbio_dir = join(self.data_dir, bcbio_dirname)
            assert isdir(bcbio_dir), f'Data directory {bcbio_dir} not found'

            if not REUSE:
                if exists(results_dir):
                    last_changed = datetime.fromtimestamp(getctime(results_dir))
                    prev_run = results_dir + '_' + last_changed.strftime('%Y_%m_%d_%H_%M_%S')
                    os.rename(results_dir, prev_run)
            safe_mkdir(results_dir)

            cmdl = f'{self.script} {bcbio_dir} -o {results_dir}'
            if parallel:
                cmdl += ' -j 10'

            ran_with_error = False
            info('-' * 100)
            try:
                check_call(cmdl)
            except subprocess.CalledProcessError:
                sys.stderr.write(f'{self.script} finished with error:\n')
                sys.stderr.write(f'{traceback.format_exc()}\n')
                ran_with_error = True
            info('-' * 100)
            info('')

        return results_dir, ran_with_error

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

    def test_one(self):
        results_dir, ran_with_error = self._run_umccrise(bcbio_dirname='bcbio_test_project', parallel=False)

        failed = False
        failed = self._check_file(failed, f'{results_dir}/log/{PROJECT}-config/{PROJECT}-template.yaml'                                    )
        failed = self._check_file(failed, f'{results_dir}/log/{PROJECT}-config/{PROJECT}.csv'                                              )
        failed = self._check_file(failed, f'{results_dir}/log/{PROJECT}-config/{PROJECT}.yaml'                                             )
        failed = self._check_file(failed, f'{results_dir}/log/{PROJECT}-data_versions.csv'                                                 )
        failed = self._check_file(failed, f'{results_dir}/log/{PROJECT}-programs.txt'                                                      )
        failed = self._check_file(failed, f'{results_dir}/{PROJECT}-multiqc_report.html'                                                   )
        for batch in BATCHES:
            failed = self._check_file(failed, f'{results_dir}/{batch}/coverage/{batch}-indexcov/index.html'                                )
            failed = self._check_file(failed, f'{results_dir}/{batch}/coverage/{batch}-normal.callable.bed'                                )
            failed = self._check_file(failed, f'{results_dir}/{batch}/coverage/{batch}-normal.depth.bed'                                   )
            failed = self._check_file(failed, f'{results_dir}/{batch}/coverage/{batch}-tumor.depth.bed'                                    )
            failed = self._check_file(failed, f'{results_dir}/{batch}/igv/{batch}-normal_mini.bam'                                         )
            failed = self._check_file(failed, f'{results_dir}/{batch}/igv/{batch}-normal_mini.bam.bai'                                     )
            failed = self._check_file(failed, f'{results_dir}/{batch}/igv/{batch}-tumor_mini.bam'                                          )
            failed = self._check_file(failed, f'{results_dir}/{batch}/igv/{batch}-tumor_mini.bam.bai'                                      )
            failed = self._check_file(failed, f'{results_dir}/{batch}/pcgr/input/{batch}-*-normal.tar.gz'                          , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{batch}/pcgr/input/{batch}-*-somatic.tar.gz'                           , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{batch}/{batch}-rmd_report.html'                                            , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/work/{batch}/rmd/afs/af_tumor.txt'                                   )
            failed = self._check_file(failed, f'{results_dir}/work/{batch}/rmd/afs/af_tumor_az300.txt'                             )
            failed = self._check_file(failed, f'{results_dir}/work/{batch}/rmd/ensemble-hg19.vcf'                                 , vcf_ignore_lines)
            failed = self._check_file(failed, f'{results_dir}/{batch}/small_variants/{batch}-normal-ensemble-cancer_genes.vcf.gz'       , vcf_ignore_lines)
            failed = self._check_file(failed, f'{results_dir}/{batch}/small_variants/{batch}-normal-ensemble-cancer_genes.vcf.gz.tbi'   , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{batch}/small_variants/{batch}-somatic-ensemble-pon_softfiltered.vcf.gz'    , vcf_ignore_lines)
            failed = self._check_file(failed, f'{results_dir}/{batch}/small_variants/{batch}-somatic-ensemble-pon_softfiltered.vcf.gz.tbi', check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{batch}/small_variants/{batch}-somatic-ensemble-pon_hardfiltered.vcf.gz'    , vcf_ignore_lines)
            failed = self._check_file(failed, f'{results_dir}/{batch}/small_variants/{batch}-somatic-ensemble-pon_hardfiltered.vcf.gz.tbi', check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{batch}/structural/{batch}-cnvkit-diagram.pdf'                              , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/work/{batch}/structural/{batch}-cnvkit-nolabels.cns'                         )
            failed = self._check_file(failed, f'{results_dir}/{batch}/structural/{batch}-sv-prioritize-manta-pass.bedpe'                   )
            failed = self._check_file(failed, f'{results_dir}/{batch}/structural/{batch}-sv-prioritize-manta-pass.ribbon.bed'              )
            failed = self._check_file(failed, f'{results_dir}/{batch}/structural/{batch}-sv-prioritize-manta-pass.tsv'                     )
            failed = self._check_file(failed, f'{results_dir}/{batch}/structural/{batch}-sv-prioritize-manta-pass.vcf'                    , vcf_ignore_lines)

        assert not ran_with_error, 'umccrise finished with error'
        assert not failed, 'some of file checks have failed'

        # failed = self._check_file(failed, join(datestamp_dir, 'cnv', 'cnvkit.filt.tsv'), wrapper=['sort'])


vcf_ignore_lines = [
    '^##bcftools_',
    '^##INFO=',
    '^##FILTER=',
    '^##contig=',
]
