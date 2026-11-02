#!/usr/bin/env python3
"""Count every requested sample, including missing outputs, without accepting QC."""
import collections
import csv
import json
from pathlib import Path
import sys


def report(root, manifest, expected):
    with Path(manifest).open() as handle:
        samples = [r['sample'] for r in csv.DictReader(handle)]
    assert len(samples) == expected and len(set(samples)) == expected
    rows = [json.loads(p.read_text()) for p in Path(root).glob('*/execution_status.json')]
    assert len({r['sample'] for r in rows}) == len(rows)
    assert {r['sample'] for r in rows} <= set(samples)
    missing = sorted(set(samples) - {r['sample'] for r in rows})
    rows += [{'sample': s, 'status': 'NO_PUBLISHED_RESULT_REQUIRES_SCHEDULER_REVIEW',
              'accepted': False, 'stage': 'unknown',
              'error': 'No worker execution record published; inspect scheduler/logs. Not a biological negative.'}
             for s in missing]
    rows.sort(key=lambda r: r['sample'])
    complete = sum(r['status'] == 'OUTPUTS_COMPLETE_REVIEW_PENDING' for r in rows)
    summary = {'samples': expected, 'accepted': 0, 'review_pending': complete,
               'failed': expected - complete, 'missing_worker_records': missing,
               'status_counts': dict(collections.Counter(r['status'] for r in rows)),
               'note': 'All failures and missing outputs remain explicit; no scientific QC acceptance.'}
    return rows, summary


if __name__ == '__main__':
    rows, summary = report(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    Path('cohort_outcomes.json').write_text(json.dumps(summary, indent=2) + '\n')
    with Path('cohort_outcomes.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=['sample','status','accepted','stage','error'],
                                delimiter='\t', extrasaction='ignore')
        writer.writeheader(); writer.writerows(rows)
