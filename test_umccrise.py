import traceback
import os
import sys
from os.path import dirname, join, exists, isfile, splitext, basename, isdir, relpath, getctime, getsize, abspath, expanduser
from datetime import datetime
import shutil
import subprocess

from ngs_utils.testing import BaseTestCase, info, check_call
from ngs_utils.utils import is_az, is_local, is_travis


REUSE = False      # Run on top of existing latest results
ONLY_DIFF = False   # Do not run, just diff the latest results against the gold standard


class Test_umccrize(BaseTestCase):
    script = 'umccrise'

    data_dir = join(dirname(__file__), BaseTestCase.data_dir, script)
    results_dir = join(dirname(__file__), BaseTestCase.results_dir, script)
    gold_standard_dir = join(dirname(__file__), BaseTestCase.gold_standard_dir, script)

    def setUp(self):
        # if is_local():
        #     info('Running locally: setting up PATH')
        #     os.environ['PATH'] = '/Users/vlad/miniconda3/envs/ngs_reporting/bin:' + expanduser('~/bin') + ':/usr/local/bin:/usr/bin:/bin:/usr/sbin:' + os.environ['PATH']
        #     info('PATH = ' + os.environ['PATH'])
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

    # def _check_var_in_datestamp(self, failed, datestamp_dir, caller):
    #     failed = self._check_file(failed, join(datestamp_dir, caller + '.PASS.txt'))
    #     failed = self._check_file(failed, join(datestamp_dir, 'var', caller + '.PASS.txt'))
    #     failed = self._check_file(failed, join(datestamp_dir, 'var', caller + '.txt'))
    #     failed = self._check_file(failed, join(datestamp_dir, 'var', caller + '.REJECT.txt'))
    #     failed = self._check_file(failed, join(datestamp_dir, 'var', caller + '.PASS.json'))
    #     if caller != 'freebayes':  # cannot merge using `bcftools merge`: > Incorrect number of AD fields (3) at chr21:11049225, cannot merge.
    #         failed = self._check_file(failed, join(datestamp_dir, 'var', caller + '.vcf.gz'), vcf_ignore_lines)
    #     return failed
    #
    # def _check_var_in_sample(self, failed, sample_dir, sample, caller):
    #     failed = self._check_file(failed, join(sample_dir, 'varAnnotate', sample + '-' + caller + '.anno.vcf.gz'), vcf_ignore_lines)
    #     failed = self._check_file(failed, join(sample_dir, 'varFilter', sample + '-' + caller + '.anno.filt.vcf.gz'), vcf_ignore_lines)
    #     failed = self._check_file(failed, join(sample_dir, 'varFilter', sample + '-' + caller + '.anno.filt.PASS.vcf.gz'), vcf_ignore_lines)
    #     failed = self._check_file(failed, join(sample_dir, 'varFilter', caller + '.PASS.json'), check_diff=False)
    #     failed = self._check_file(failed, join(sample_dir, 'varFilter', caller + '.PASS.txt'))
    #     failed = self._check_file(failed, join(sample_dir, 'varFilter', caller + '.txt'))
    #     failed = self._check_file(failed, join(sample_dir, 'varFilter', caller + '.REJECT.txt'))
    #     return failed

    def test_one(self):
        results_dir, ran_with_error = self._run_umccrise(bcbio_dirname="bcbio_test_project", parallel=False)

        # datestamp_name = '2014-08-13_' + name  # TODO: change after run is finished
        # datestamp_dir = join(bcbio_proj_dir, 'final', datestamp_name)

        failed = False
        # failed = self._check_file(failed, join(bcbio_proj_dir, 'config', 'run_info_ExomeSeq.yaml'))
        # failed = self._check_file(failed, join(datestamp_dir, 'NGv3.chr21.4col.clean.sorted.bed'))
        # failed = self._check_file(failed, join(datestamp_dir, 'report.html'), check_diff=False)
        # failed = self._check_file(failed, join(datestamp_dir, 'reports', 'call_vis.html'), wrapper=html_wrapper, check_diff=False)
        # failed = self._check_file(failed, join(datestamp_dir, 'cnv', 'seq2c.tsv'), wrapper=['sort'])
        # failed = self._check_file(failed, join(datestamp_dir, 'cnv', 'seq2c.filt.tsv'), wrapper=['sort'])
        # failed = self._check_file(failed, join(datestamp_dir, 'cnv', 'seq2c-coverage.tsv'), wrapper=['sort'])
        # failed = self._check_file(failed, join(datestamp_dir, 'cnv', 'cnvkit.tsv'), wrapper=['sort'])
        # failed = self._check_file(failed, join(datestamp_dir, 'cnv', 'cnvkit.filt.tsv'), wrapper=['sort'])
        # for caller, samples in callers.items():
        #     if all(s.endswith('-germline') for s in samples):
        #         failed = self._check_var_in_datestamp(failed, datestamp_dir, caller + '-germline')
        #     else:
        #         failed = self._check_var_in_datestamp(failed, datestamp_dir, caller)
        #     for sample in samples:
        #         failed = self._check_file(failed, join(datestamp_dir, 'reports', sample + '.html'), wrapper=html_wrapper, check_diff=False)
        #         sample_dir = join(bcbio_proj_dir, 'final', sample)
        #         failed = self._check_var_in_sample(failed, sample_dir, sample, caller)

        assert not ran_with_error, 'umccrise finished with error'
        assert not failed, 'some of file checks have failed'

        # if exists(join(bcbio_proj_dir, 'work')):
        #     shutil.rmtree(join(bcbio_proj_dir, 'work'))


# Find and parse all elements containing json data, put data into a list and dumps the result.
# The resulting text is unique per json data, so we can run simple `diff` on them.
html_wrapper = [
    'grep', '-A1', '<div id=".*_json">', '|', 'grep', '-v', '<div id=".*_json">', '|',
    'python', '-c',
        'import sys, json; '
        'sys.stdout.write(json.dumps([json.loads(el) for el in sys.stdin.read().split(\'--\')], '
                                     'indent=2, sort_keys=True))'
]

vcf_ignore_lines = [
    'bcftools_annotateVersion',
    'bcftools_annotateCommand',
    '^##INFO=',
    '^##FILTER=',
    '^##contig=',
]
