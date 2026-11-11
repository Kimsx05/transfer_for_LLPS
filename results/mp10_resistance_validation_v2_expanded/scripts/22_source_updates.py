from env import *
import importlib,urllib.request,zipfile
get=importlib.import_module('01_inputs').get
p=S/'GSE285648_Figure1.jpg'
if p.exists() and p.stat().st_size<100000:p.rename(S/'GSE285648_Figure1_INCOMPLETE.jpg')
get(('GSE285648_Figure1.jpg','https://pmc-oa-opendata.s3.amazonaws.com/PMC13042512.1/MOL2-20-779-g001.jpg'))
for pid in ['PMC11094426','PMC10591055']:
 u='https://pmc-oa-opendata.s3.amazonaws.com/?list-type=2&prefix='+pid+'.1/'
 try:
  xml=urllib.request.urlopen(u,timeout=40).read().decode();(S/(pid+'_assets.xml')).write_text(xml)
  for key in re.findall('<Key>(.*?)</Key>',xml):
   if 'suppst1-7' in key or key.endswith('SourceData.xlsx'):get((key.split('/')[-1],'https://pmc-oa-opendata.s3.amazonaws.com/'+key))
 except Exception as e:print(pid,e)
# DTP publisher supplemental tables have .xls as well as .xlsx in older Cell archives.
for i in [7,8]:
 name=f'DTP_alternate_mmc{i}.xls';url=f'https://ars.els-cdn.com/content/image/1-s2.0-S009286742031535X-mmc{i}.xls';get((name,url))
