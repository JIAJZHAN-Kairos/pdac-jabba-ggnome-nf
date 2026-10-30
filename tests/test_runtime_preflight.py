import contextlib
import hashlib
import importlib.util
import io
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('runner', Path(__file__).resolve().parents[1] / 'bin/run_sample.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


class RuntimePreflight(unittest.TestCase):
    def fixture(self, root):
        pins = {}
        for name in ['cplex/include/ilcplex/cplex.h', 'cplex/lib/x86-64_linux/static_pic/libcplex.a', 'cplex/bin/x86-64_linux/cplex']:
            p = root / name
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(b'test-only fixture')
            pins[name] = hashlib.sha256(p.read_bytes()).hexdigest()
        return pins

    def test_bin_only_deployment_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            (p / 'cplex/bin/x86-64_linux').mkdir(parents=True)
            with self.assertRaisesRegex(AssertionError, 'Incomplete CPLEX installation'):
                runner.verify_cplex(p, {'cplex/include/ilcplex/cplex.h': '0' * 64})

    def test_empty_directories_do_not_pass(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            for name in ['include', 'lib', 'bin/x86-64_linux']:
                (p / 'cplex' / name).mkdir(parents=True)
            with self.assertRaisesRegex(AssertionError, 'Missing CPLEX file'):
                runner.verify_cplex(p, {'cplex/include/ilcplex/cplex.h': '0' * 64})

    def test_original_hashes_required(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            pins = self.fixture(p)
            runner.verify_cplex(p, pins)
            (p / 'cplex/lib/x86-64_linux/static_pic/libcplex.a').write_bytes(b'changed')
            with self.assertRaisesRegex(AssertionError, 'CPLEX source hash mismatch'):
                runner.verify_cplex(p, pins)

    def test_failed_subprocess_keeps_log_and_raises(self):
        with tempfile.TemporaryDirectory() as d:
            log = Path(d) / 'model.log'
            with contextlib.redirect_stdout(io.StringIO()) as live:
                with self.assertRaises(subprocess.CalledProcessError) as caught:
                    runner.run_logged([sys.executable, '-c', 'print("native error retained"); raise SystemExit(7)'], log, d)
            self.assertEqual(caught.exception.returncode, 7)
            self.assertEqual(log.read_text(), 'native error retained\n')
            self.assertEqual(live.getvalue(), log.read_text())


if __name__ == '__main__':
    unittest.main()
