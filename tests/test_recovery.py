import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


def load(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).resolve().parents[1] / 'bin' / (name + '.py'))
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


support = load('recovery_support')
outcomes = load('report_recovery_outcomes')


class Recovery(unittest.TestCase):
    def fixture(self):
        meta = {'sample': 'TEST', 'recovery_profile': 'recovery',
                **{k + '_sha256': '0' * 64 for k in ['cobalt','sv','cnv','purity']}}
        spec = {'sample': 'TEST', 'profile': 'recovery', 'mode': 'NATIVE_WARMSTART',
                'cpus': 12, 'memory_GiB': 64, 'solver_memory_GB': 48,
                'tilim': 21600, 'walltime_hours': 48,
                'original_input_sha256': {k: '0' * 64 for k in ['cobalt','sv','cnv','purity']},
                'checkpoint_files': [{'file': n, 'uri': 's3://test/' + n, 'sha256': '0' * 64, 'bytes': 1}
                                     for n in support.CHECKPOINTS]}
        return meta, spec

    def test_input_identity_cannot_change_during_recovery(self):
        meta, spec = self.fixture()
        support.validate_spec(meta, spec)
        spec['original_input_sha256']['sv'] = '1' * 64
        with self.assertRaises(AssertionError):
            support.validate_spec(meta, spec)

    def test_partial_or_preprocessing_checkpoint_rejected(self):
        meta, spec = self.fixture()
        spec['checkpoint_files'][0]['file'] = 'mipstart.rds'
        with self.assertRaises(AssertionError):
            support.validate_spec(meta, spec)

    def test_resource_change_must_match_named_profile(self):
        meta, spec = self.fixture()
        spec['cpus'] = 64
        with self.assertRaisesRegex(AssertionError, 'Unapproved resource'):
            support.validate_spec(meta, spec)

    def test_missing_output_remains_failed_not_zero_or_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / 'samples.csv'
            with manifest.open('w') as f:
                w = csv.DictWriter(f, fieldnames=['sample']); w.writeheader()
                w.writerows([{'sample': 'COMPLETE'}, {'sample': 'LOST'}])
            d = root / 'COMPLETE'; d.mkdir()
            (d / 'execution_status.json').write_text(json.dumps({'sample': 'COMPLETE',
                'status': 'OUTPUTS_COMPLETE_REVIEW_PENDING', 'accepted': False}))
            rows, summary = outcomes.report(root, manifest, 2)
            self.assertEqual((summary['review_pending'], summary['failed'], summary['accepted']), (1, 1, 0))
            self.assertEqual(summary['missing_worker_records'], ['LOST'])
            self.assertFalse(next(r for r in rows if r['sample'] == 'LOST')['accepted'])
