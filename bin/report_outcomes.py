#!/usr/bin/env python3
import collections
import csv
import json
from pathlib import Path
import sys

paths = sorted(Path(sys.argv[1]).glob('*/execution_status.json'))
rows = [json.loads(p.read_text()) for p in paths]
assert len(rows) == int(sys.argv[2])
assert len({r['sample'] for r in rows}) == len(rows)
summary = {'samples': len(rows), 'accepted': 0, 'review_pending': sum(r['status'] == 'OUTPUTS_COMPLETE_REVIEW_PENDING' for r in rows),
           'failed': sum(r['status'] != 'OUTPUTS_COMPLETE_REVIEW_PENDING' for r in rows),
           'status_counts': dict(collections.Counter(r['status'] for r in rows)),
           'note': 'Execution is not actual numerical/image QC acceptance. Failed samples are not biological negatives.'}
Path('cohort_outcomes.json').write_text(json.dumps(summary, indent=2) + '\n')
with Path('cohort_outcomes.tsv').open('w') as f:
    w = csv.DictWriter(f, fieldnames=['sample','status','accepted','stage','error'], delimiter='\t', extrasaction='ignore')
    w.writeheader(); w.writerows(rows)
