from pathlib import Path
import csv
import sys
root=Path(sys.argv[1]).resolve();root.mkdir(parents=True,exist_ok=True);rows=[]
for sample in ['STUB_A','STUB_B']:
 row={'sample':sample,'cohort':'TEST'}
 for col,key in [('cobalt_ratio','cobalt'),('sv_vcf','sv'),('purple_cnv','cnv'),('purple_purity','purity')]:
  p=root/(sample+'.'+col);p.write_text('STUB_ONLY\n');row[col]=str(p);row[key+'_sha256']='0'*64;row[key+'_bytes']=p.stat().st_size
 rows.append(row)
with (root/'samplesheet.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
