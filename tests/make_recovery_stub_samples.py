from pathlib import Path
import csv
import hashlib
import json
import sys

root = Path(sys.argv[1]).resolve(); root.mkdir(parents=True, exist_ok=True)
rows = []
for profile in ['standard', 'recovery', 'hard']:
    sample = 'STUB_' + profile.upper()
    row = {'sample': sample, 'cohort': 'TEST', 'recovery_profile': profile}
    for column, key in [('cobalt_ratio','cobalt'),('sv_vcf','sv'),('purple_cnv','cnv'),('purple_purity','purity')]:
        path = root / (sample + '.' + column); path.write_text('STUB_ONLY\n')
        row[column] = str(path); row[key + '_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
        row[key + '_bytes'] = path.stat().st_size
    spec = root / (sample + '.recovery_spec.json')
    spec.write_text(json.dumps({'sample': sample, 'stub_only': True}))
    row['recovery_spec'] = str(spec); row['recovery_spec_sha256'] = hashlib.sha256(spec.read_bytes()).hexdigest()
    rows.append(row)
with (root / 'samplesheet.csv').open('w') as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
