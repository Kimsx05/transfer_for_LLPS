import runpy
from pathlib import Path
ns=runpy.run_path(str(Path(__file__).with_name('01_inputs.py')));globals().update({k:v for k,v in ns.items() if not k.startswith('__')})
jobs=[]
for p in S.glob('GSE*_s3.xml'):
 for k in re.findall('<Key>(.*?)</Key>',p.read_text()):
  if k.endswith(('.txt','.xlsx','.docx')):jobs.append((k.split('/')[-1],'https://pmc-oa-opendata.s3.amazonaws.com/'+k))
for i in [7,8]:jobs.append((f'DTP_mmc{i}.xlsx',f'https://ars.els-cdn.com/content/image/1-s2.0-S009286742031535X-mmc{i}.xlsx'))
for item in jobs:get(item)
