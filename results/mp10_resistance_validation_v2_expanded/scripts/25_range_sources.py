from env import *
import urllib.request,time,zipfile

def fetch(name,url):
 p=S/name
 if p.exists():return p
 req=urllib.request.Request(url,method='HEAD');h=urllib.request.urlopen(req,timeout=40);n=int(h.headers['Content-Length']);out=p.with_name(p.name+'.range_part');print(name,n,flush=True)
 with out.open('wb') as f:
  for a in range(0,n,262144):
   b=min(n-1,a+262143)
   for attempt in range(8):
    try:
     req=urllib.request.Request(url+'?part='+str(a),headers={'Range':f'bytes={a}-{b}'});r=urllib.request.urlopen(req,timeout=40);data=r.read();assert len(data)==b-a+1;f.write(data);break
    except Exception as e:
     print('retry',a,str(e)[:90],flush=True);time.sleep(1)
   else:raise RuntimeError('range unrecoverable')
 if name.endswith('.xlsx'):
  with zipfile.ZipFile(out) as z:assert z.testzip() is None
 out.replace(p);return p
p=fetch('PMC13042512.1.xml','https://pmc-oa-opendata.s3.amazonaws.com/PMC13042512.1/PMC13042512.1.xml');t=p.read_text()
for fig in re.findall(r'<fig .*?</fig>',t,re.S):
 label=re.search('<label>(.*?)</label>',fig)
 if label and re.sub('<[^>]+>','',label[1]).strip()=='Fig. 1':
  keys=re.findall('xlink:href="([^"]+)"',fig);print('ACTUAL FIG1',keys,flush=True)
  if keys:fetch('GSE285648_IC50_Figure1.jpg','https://pmc-oa-opendata.s3.amazonaws.com/PMC13042512.1/'+keys[-1]+('.jpg' if not keys[-1].endswith('.jpg') else ''))
fetch('can-23-1976_supplementary_tables_suppst1-7.xlsx','https://pmc-oa-opendata.s3.amazonaws.com/PMC11094426.1/can-23-1976_supplementary_tables_suppst1-7.xlsx')
