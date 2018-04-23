import traceback
import os
import sys
from os.path import dirname, join, exists, isfile, splitext, basename, isdir, relpath, getctime, getsize, abspath, expanduser
from datetime import datetime
import shutil
import subprocess

from ngs_utils.testing import BaseTestCase, info, check_call, vcf_ignore_lines, swap_output
from ngs_utils.utils import is_az, is_local, is_travis
from ngs_utils.file_utils import safe_mkdir


TUMORS = ['cup_tissue']
BATCHES = ['cup']
PROJECT = 'cup_sc932'


class Test_umccrise(BaseTestCase):
    script = 'umccrise'

    data_dir = join(dirname(__file__), BaseTestCase.data_dir)
    results_dir = join(dirname(__file__), BaseTestCase.results_dir)
    gold_standard_dir = join(dirname(__file__), BaseTestCase.gold_standard_dir)

    reuse = False  # Run on top of existing latest results. Also controlled with TEST_REUSE
    only_diff = False  # Do not run, just diff the latest results against the gold standard. Also controlled with TEST_ONLY_DIFF

    def setUp(self):
        BaseTestCase.setUp(self)

    def _run_umccrise(self, bcbio_dirname, parallel=False):
        results_dir = join(self.results_dir, bcbio_dirname)
        bcbio_dir = join(self.data_dir, bcbio_dirname)
        cmdl = f'{self.script} {bcbio_dir} -o {results_dir}'
        if parallel:
            cmdl += ' -j 10'
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

    def test_one(self):
        results_dir = self._run_umccrise(bcbio_dirname='bcbio_test_project', parallel=False)

        failed = False
        failed = self._check_file(failed, f'{results_dir}/log/{PROJECT}-config/{PROJECT}-template.yaml'                                                     )
        failed = self._check_file(failed, f'{results_dir}/log/{PROJECT}-config/{PROJECT}.csv'                                                               )
        failed = self._check_file(failed, f'{results_dir}/log/{PROJECT}-config/{PROJECT}.yaml'                                                              )
        failed = self._check_file(failed, f'{results_dir}/log/{PROJECT}-data_versions.csv'                                                                  )
        failed = self._check_file(failed, f'{results_dir}/log/{PROJECT}-programs.txt'                                                                       )
        failed = self._check_file(failed, f'{results_dir}/{PROJECT}-multiqc_report.html'                                                                    )
        for T, B in zip(TUMORS, BATCHES):
            key = f'{B}__{T}'
            failed = self._check_file(failed, f'{results_dir}/{key}/coverage/{key}-indexcov/index.html'                               , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/coverage/{key}-normal.callable.bed'                                                 )
            failed = self._check_file(failed, f'{results_dir}/{key}/coverage/{key}-normal.depth.bed'                                                    )
            failed = self._check_file(failed, f'{results_dir}/{key}/coverage/{key}-tumor.depth.bed'                                                     )
            failed = self._check_file(failed, f'{results_dir}/{key}/igv/{key}-roi.bed'                                                                  )

            failed = self._check_file(failed, f'{results_dir}/{key}/igv/{key}-normal.mini.bam'                                        , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/igv/{key}-normal.mini.bam.bai'                                    , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/igv/{key}-tumor.mini.bam'                                         , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/igv/{key}-tumor.mini.bam.bai'                                     , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/pcgr/input/{key}-normal.tar.gz'                                   , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/pcgr/input/{key}-somatic.tar.gz'                                  , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/{key}-rmd_report.html'                                            , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/work/{key}/rmd/afs/af_tumor.txt'                                                            )
            failed = self._check_file(failed, f'{results_dir}/work/{key}/rmd/afs/af_tumor_az300.txt'                                                      )
            failed = self._check_file(failed, f'{results_dir}/work/{key}/rmd/ensemble-hg19.vcf'                                         , vcf_ignore_lines)
            failed = self._check_file(failed, f'{results_dir}/{key}/small_variants/{key}-normal-ensemble-cancer_genes.vcf.gz'         , vcf_ignore_lines)
            failed = self._check_file(failed, f'{results_dir}/{key}/small_variants/{key}-normal-ensemble-cancer_genes.vcf.gz.tbi'     , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/small_variants/{key}-somatic-ensemble-pon_softfiltered.vcf.gz'    , vcf_ignore_lines)
            failed = self._check_file(failed, f'{results_dir}/{key}/small_variants/{key}-somatic-ensemble-pon_softfiltered.vcf.gz.tbi', check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/small_variants/{key}-somatic-ensemble-pon_hardfiltered.vcf.gz'    , vcf_ignore_lines)
            failed = self._check_file(failed, f'{results_dir}/{key}/small_variants/{key}-somatic-ensemble-pon_hardfiltered.vcf.gz.tbi', check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/work/{key}/structural/{key}-cnvkit-nolabels.cns'                                          )
            failed = self._check_file(failed, f'{results_dir}/{key}/structural/{key}-cnvkit-diagram.pdf'                              , check_diff=False)
            failed = self._check_file(failed, f'{results_dir}/{key}/structural/{key}-sv-prioritize-manta-pass.bedpe'                                    )
            failed = self._check_file(failed, f'{results_dir}/{key}/structural/{key}-sv-prioritize-manta-pass.ribbon.bed'                               )
            failed = self._check_file(failed, f'{results_dir}/{key}/structural/{key}-sv-prioritize-manta-pass.tsv'                                      )
            failed = self._check_file(failed, f'{results_dir}/{key}/structural/{key}-sv-prioritize-manta-pass.vcf'                    , vcf_ignore_lines)

        assert not failed, 'some of file checks have failed'

        # failed = self._check_file(failed, join(datestamp_dir, 'cnv', 'cnvkit.filt.tsv'), wrapper=['sort'])
