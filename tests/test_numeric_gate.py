import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

spec=importlib.util.spec_from_file_location('runner',Path(__file__).resolve().parents[1]/'bin/run_sample.py')
runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)

class NumericGate(unittest.TestCase):
    def check(self,status,gap,cn='2',balance=0):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            (p/'optimization_report.tsv').write_text(f'status\tepgap\n{status}\t{gap}\n')
            (p/'nodes.tsv').write_text(f'start\tend\tcn\n1\t100\t{cn}\n')
            (p/'sample_QC.json').write_text(json.dumps({'max_abs_internal_balance_residual':balance}))
            return runner.strict_numeric(p)
    def test_107_with_zero_gap_rejected(self):
        with self.assertRaises(AssertionError):self.check(107,0)
    def test_no_tolerance_above_original_gap(self):
        with self.assertRaises(AssertionError):self.check(101,0.000100001)
    def test_unknown_and_fractional_CN(self):
        self.assertEqual(self.check(102,0.0001,'NA')['unknown_node_bp'],100)
        with self.assertRaises(AssertionError):self.check(101,0,'1.5')
    def test_internal_balance_must_be_exact_zero(self):
        with self.assertRaises(AssertionError):self.check(101,0,balance=1e-9)
    def test_valid_final_status(self):self.assertEqual(self.check(102,0.000099)['final_components'],1)

if __name__=='__main__':unittest.main()
