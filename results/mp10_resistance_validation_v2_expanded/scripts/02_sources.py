import runpy
globals().update(runpy.run_path(str(__file__).replace('02_sources.py','01_inputs.py')))
jobs=[]
for name in ['CisSig_article.html','Taber2020.html']:
 for u in re.findall('href="([^"]+)"',(S/name).read_text()):
  if 'static-content.springer' in u:jobs.append((u.split('/')[-1],u))
for pmc in ['PMC5039423','PMC5085609','PMC8437243','PMC10992486']:
 jobs.append((pmc+'.html',f'https://pmc.ncbi.nlm.nih.gov/articles/{pmc}/?pdf=1'))
 jobs.append((pmc+'_s3.xml',f'https://pmc-oa-opendata.s3.amazonaws.com/?list-type=2&prefix={pmc}.'))
for item in jobs:get(item)
