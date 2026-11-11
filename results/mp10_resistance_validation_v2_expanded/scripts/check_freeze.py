from env import *
p=B/'00_plan_audit/frozen_config.json';assert hashlib.sha256(p.read_bytes()).hexdigest()==(B/'00_plan_audit/frozen_config.sha256').read_text().strip();cfg=json.loads(p.read_text())
for filename,digest in cfg['locked_files'].items():assert hashlib.sha256(Path(filename).read_bytes()).hexdigest()==digest,filename
for _,r in read('00_plan_audit/input_manifest.tsv').iterrows():assert hashlib.sha256(Path(r.path).read_bytes()).hexdigest()==r.sha256,r.path
print('FROZEN_INPUTS_VERIFIED')

p=B/'00_plan_audit/amendment_manifest.tsv'
if p.exists():
 for _,r in pd.read_csv(p,sep='\t').iterrows():assert hashlib.sha256(Path(r.path).read_bytes()).hexdigest()==r.sha256,r.path

p=B/'00_plan_audit/epi_input_manifest.tsv'
if p.exists():
 for _,r in pd.read_csv(p,sep='\t').iterrows():
  h=hashlib.sha256()
  with Path(r.path).open('rb') as f:
   while block:=f.read(2**20):h.update(block)
  assert h.hexdigest()==r.sha256,r.path
