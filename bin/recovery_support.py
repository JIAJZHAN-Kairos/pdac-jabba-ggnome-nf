"""Stage a provenance-checked native recovery without changing original inputs."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

PROFILES = {
    'standard': (4, 16, 12, 6000, 8),
    'recovery': (12, 64, 48, 21600, 48),
    'hard': (24, 96, 72, 43200, 48),
}
CHECKPOINTS = {'karyograph.rds', 'jabba.raw.rds', 'edge.nudge.rds',
               'ab.force.rds', 'ab.exclude.rds'}


def digest(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def validate_spec(meta, spec):
    assert spec['sample'] == meta['sample'], 'Recovery sample mismatch'
    assert spec['profile'] == meta['recovery_profile'], 'Recovery profile mismatch'
    assert spec['mode'] in {'FIRST_RUN', 'LONG_FRESH', 'NATIVE_WARMSTART'}
    actual = tuple(spec[k] for k in ['cpus', 'memory_GiB', 'solver_memory_GB',
                                    'tilim', 'walltime_hours'])
    assert actual == PROFILES[spec['profile']], 'Unapproved resource profile'
    assert (spec['mode'] == 'FIRST_RUN') == (spec['profile'] == 'standard')
    for key in ['cobalt', 'sv', 'cnv', 'purity']:
        assert spec['original_input_sha256'][key] == meta[key + '_sha256']
    files = spec['checkpoint_files']
    if spec['mode'] == 'NATIVE_WARMSTART':
        assert len(files) == 5 and {x['file'] for x in files} == CHECKPOINTS
        assert all(x['uri'].startswith('s3://') and len(x['sha256']) == 64
                   and x['bytes'] > 0 for x in files)
    else:
        assert not files, 'Fresh attempts must not masquerade as warm starts'


def prepare_recovery(meta, root, spec_path, script_path):
    assert digest(spec_path) == meta['recovery_spec_sha256'], 'Recovery spec hash mismatch'
    spec = json.loads(Path(spec_path).read_text())
    validate_spec(meta, spec)
    assert digest(script_path) == spec['recovery_adapter_sha256'], 'Recovery adapter hash mismatch'
    shutil.copyfile(spec_path, root / 'recovery_spec.json')
    target = root / 'scripts' / 'run_sample_recovery.R'
    shutil.copyfile(script_path, target)
    if spec['checkpoint_files']:
        previous = root / 'previous_native' / 'graphs' / meta['sample']
        previous.mkdir(parents=True)
        for item in spec['checkpoint_files']:
            path = previous / item['file']
            subprocess.run(['aws', 's3', 'cp', '--only-show-errors', item['uri'], str(path)], check=True)
            assert path.stat().st_size == item['bytes'], 'Checkpoint size mismatch'
            assert digest(path) == item['sha256'], 'Checkpoint hash mismatch'
    return spec, target
