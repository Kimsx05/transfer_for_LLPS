exec(open(__file__.replace('01_inputs.py','00_acquire.py')).read().split("if __name__")[0])
import subprocess,urllib.parse,zipfile
# Downloads are retried and validated; partial bytes never become final inputs.
def get(item):
 name,url=item;p=S/name
 if p.exists():return dict(file=name,url=url,status='EXISTS',bytes=p.stat().st_size)
 tmp=p.with_name(p.name+'.part')
 r=subprocess.run(['curl','-fL','--retry','3','--retry-all-errors','--connect-timeout','20','--max-time','240','-o',str(tmp),url],capture_output=True)
 ok=r.returncode==0
 if ok:
  try:
   if name.endswith('.gz'):
    with gzip.open(tmp,'rb') as f:
     while f.read(2**20):pass
   if name.endswith('.zip') or name.endswith('.xlsx'):
    with zipfile.ZipFile(tmp) as z:assert z.testzip() is None
  except Exception:ok=False
 if ok:tmp.replace(p)
 out=dict(file=name,url=url,status='OK' if ok else ('FAILED '+r.stderr.decode()[-300:]),bytes=p.stat().st_size if p.exists() else 0);print(name,out['status'],out['bytes'],flush=True);return out
if __name__=='__main__':
 jobs=[]
 for g in cell+clinical:
  root=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{g[:-3]}nnn/{g}'
  jobs.append((g+'_metadata.soft',f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={g}&targ=all&view=full&form=text'))
  p=S/(g+'_suppl.html')
  if p.exists():
   for name in re.findall(r'href="([^"]+)"',p.read_text()):
    if name.startswith(g) and ('RAW' not in name or g in ['GSE231835','GSE313089']):jobs.append((name,root+'/suppl/'+name))
  p=S/(g+'_matrix.html')
  if p.exists() and g in ['GSE58624','GSE77883','GSE52219','GSE244266','GSE87304','GSE69795','GSE5287','GSE124305','GSE48276']:
   for name in re.findall(r'href="([^"]+)"',p.read_text()):
    if name.endswith('.gz'):jobs.append((name,root+'/matrix/'+name))
 for screen in ['GDSC1','GDSC2']:jobs.append((screen+'_fitted_dose_response_27Oct23.xlsx',f'https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/{screen}_fitted_dose_response_27Oct23.xlsx'))
 jobs += [('Cell_line_RMA_proc_basalExp.txt.zip','https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Data/preprocessed/Cell_line_RMA_proc_basalExp.txt.zip'),('Cell_Lines_Details.xlsx','https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/Cell_Lines_Details.xlsx')]
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:res=list(ex.map(get,jobs))
 (B/'00_plan_audit/acquisition_processed.json').write_text(json.dumps(res,indent=2))
