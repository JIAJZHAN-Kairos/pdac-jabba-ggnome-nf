#!/usr/bin/env python3
"""Portable orchestration only; pinned Gadi scientific scripts remain unchanged."""
import csv
import datetime
import hashlib
import json
import math
import os
import resource
from pathlib import Path
import shutil
import subprocess
import sys


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def strict_numeric(out):
    with (out / 'optimization_report.tsv').open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    assert rows and all(int(x['status']) in (101, 102) and math.isfinite(float(x['epgap']))
                        and 0 <= float(x['epgap']) <= 1e-4 for x in rows), 'Final status/gap gate failed'
    with (out / 'nodes.tsv').open() as f:
        nodes = list(csv.DictReader(f, delimiter='\t'))
    unknown_bp = 0
    unknown_ranges = []
    for n in nodes:
        try:
            cn = float(n['cn'])
        except ValueError:
            cn = math.nan
        if math.isfinite(cn):
            assert cn >= 0 and cn.is_integer(), 'Known CN must be nonnegative integer'
        else:
            unknown_bp += int(n['end']) - int(n['start']) + 1
            unknown_ranges.append(n)
    qc = json.loads((out / 'sample_QC.json').read_text())
    assert qc['max_abs_internal_balance_residual'] == 0, 'Internal native balance must be zero'
    return {'final_components': len(rows), 'max_gap': max(float(x['epgap']) for x in rows),
            'unknown_node_bp': unknown_bp, 'unknown_node_ranges': unknown_ranges, 'max_abs_internal_balance_residual': 0}


def main():
    meta_path, scripts_path, pins_path, root_path, *inputs = sys.argv[1:]
    meta = json.loads(Path(meta_path).read_text())
    scripts_path = Path(scripts_path).resolve()
    pins = json.loads(Path(pins_path).read_text())
    sample = meta['sample']
    root = Path(root_path).resolve()
    root.mkdir(parents=True, exist_ok=False)
    state = {'sample': sample, 'wave_id': meta['wave_id'], 'accepted': False,
             'status': 'RUNNING', 'started_at_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
             'AWS_BATCH_JOB_ID': os.environ.get('AWS_BATCH_JOB_ID'), 'resource_profile': '4CPU16GiB_CPLEX12GB',
             'review_required': ['original numerical QC', 'all target PNGs', 'actual coverage support', 'architecture limitations']}
    try:
        for name, expected in pins.items():
            assert sha(scripts_path / name) == expected, 'Scientific source hash mismatch: ' + name
        shutil.copytree(scripts_path, root / 'scripts')
        for name in ['inputs','logs','prepared','graph_inputs','graphs']:
            (root / name).mkdir()
        suffixes = ['cobalt.ratio.tsv.gz','esvee.somatic.vcf.gz','purple.cnv.somatic.tsv','purple.purity.tsv']
        keys = ['cobalt','sv','cnv','purity']
        audit = []
        for source, suffix, key in zip(inputs, suffixes, keys):
            source = Path(source).resolve()
            assert source.stat().st_size == int(meta[key + '_bytes']), 'Input size mismatch'
            assert sha(source) == meta[key + '_sha256'], 'Input SHA256 mismatch'
            target = root / 'inputs' / (sample + '.' + suffix)
            os.symlink(source, target)
            audit.append({'file': target.name, 'bytes': source.stat().st_size, 'sha256': meta[key + '_sha256']})
        with (root / 'input_sha256.tsv').open('w') as f:
            w = csv.DictWriter(f, fieldnames=['file','bytes','sha256'], delimiter='\t'); w.writeheader(); w.writerows(audit)
        with (root / 'inputs' / (sample + '.purple.purity.tsv')).open() as f:
            purity = list(csv.DictReader(f, delimiter='\t'))
        assert len(purity) == 1
        metadata = dict(sample=sample, cohort=meta['cohort'], purity=purity[0]['purity'],
                        ploidy=purity[0]['ploidy'], sex=purity[0]['gender'])
        with (root / 'pilot_samples.tsv').open('w') as f:
            w = csv.DictWriter(f, fieldnames=list(metadata), delimiter='\t'); w.writeheader(); w.writerow(metadata)
        state['input_provenance'] = audit
        state['scientific_source_sha256'] = pins
        state['stage'] = 'runtime'
        with (root / 'logs' / 'runtime.log').open('w') as f:
            subprocess.run(['Rscript','--vanilla','-e','print(R.version.string); print(.libPaths()); library(JaBbA); library(gGnome); print(packageVersion("JaBbA")); print(packageVersion("gGnome"))'], stdout=f, stderr=subprocess.STDOUT, check=True)
        stages = [('input_QC', ['python3','scripts/check_inputs.py']),
                  ('input', ['Rscript','--vanilla','scripts/build_graph_inputs.R']),
                  ('model', ['Rscript','--vanilla','scripts/run_sample_16gb.R']),
                  ('postfit', ['Rscript','--vanilla','scripts/finish_sample.R']),
                  ('CN_review', ['python3','scripts/review_cn.py']),
                  ('native_plot', ['Rscript','--vanilla','scripts/plot_native.R'])]
        for stage, command in stages:
            state['stage'] = stage
            (root / 'execution_status.json').write_text(json.dumps(state, indent=2) + '\n')
            with (root / 'logs' / (sample + '.' + stage + '.log')).open('w') as f:
                subprocess.run(command + [str(root), sample], cwd=root, stdout=f, stderr=subprocess.STDOUT, check=True)
        state['stage'] = 'strict_numeric_gate'
        state['strict_numeric'] = strict_numeric(root / 'graphs' / sample)
        state['visual_artifact_sha256'] = {str(p.relative_to(root)): sha(p) for p in sorted((root / 'graphs' / sample).rglob('*.png'))}
        assert state['visual_artifact_sha256']
        state['status'] = 'OUTPUTS_COMPLETE_REVIEW_PENDING'
    except Exception as exc:
        state.update(status='FAILED_REQUIRES_REVIEW', error_type=type(exc).__name__, error=str(exc),
                     subprocess_exit_code=getattr(exc,'returncode',None))
    finally:
        state['resource_observations'] = {'max_child_RSS_KiB': resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss, 'note': 'Child RSS is not total container peak; raw cgroup readings retained separately'}
        for path in ['/sys/fs/cgroup/memory.max','/sys/fs/cgroup/memory.peak','/sys/fs/cgroup/cpu.max','/sys/fs/cgroup/memory/memory.limit_in_bytes','/sys/fs/cgroup/memory/memory.max_usage_in_bytes','/sys/fs/cgroup/cpu/cpu.cfs_quota_us','/sys/fs/cgroup/cpu/cpu.cfs_period_us']:
            if Path(path).is_file():
                state['resource_observations'][path] = Path(path).read_text().strip()
        state['finished_at_UTC'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        (root / 'execution_status.json').write_text(json.dumps(state, indent=2) + '\n')
        manifest = {str(p.relative_to(root)): sha(p) for p in sorted(root.rglob('*')) if p.is_file() and not p.is_symlink()}
        (root / 'output_sha256.json').write_text(json.dumps(manifest, indent=2) + '\n')
    # The following report and execution gate mark scientific failures after their
    # original files/checkpoints have been staged out and published by Nextflow.
    print(json.dumps(state), flush=True)


if __name__ == '__main__':
    main()
