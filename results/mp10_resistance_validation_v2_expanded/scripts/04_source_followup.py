import runpy
from pathlib import Path
ns=runpy.run_path(str(Path(__file__).with_name('01_inputs.py')));globals().update({k:v for k,v in ns.items() if not k.startswith('__')})
jobs=[('hgnc_complete_set.txt','https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt'),('DTP_TableS6.xlsx','https://pmc.ncbi.nlm.nih.gov/articles/instance/8437243/bin/NIHMS1732990-supplement-14.xlsx'),('Tanaka.html','https://www.oncotarget.com/article/11427/text/'),('Tanaka_PMC.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5085609/fullTextXML'),('CrM_elsevier.xml','https://api.elsevier.com/content/article/doi/10.1016/j.canlet.2024.217339?httpAccept=text/xml'),('ITGB4_article.html','https://www.nature.com/articles/s41416-026-03364-7'),('T24R2_article.html','https://academic.oup.com/narcancer/article/5/4/zcad054/7424730'),('S1314_github_tree.json','https://api.github.com/repos/ZB-WCM/GMLF_response_NAC_MIBC/git/trees/main?recursive=1'),('v1_reference_tree.json','https://api.github.com/repos/Kimsx05/transfer_for_neibor/git/trees/a1aaeba7cdd2b7398c6798a295253b3062668e49?recursive=1')]
for p in S.glob('PMC*_s3.xml'):
 for k in re.findall('<Key>(.*?)</Key>',p.read_text()):
  if k.endswith(('.txt','.docx')):jobs.append((k.split('/')[-1],'https://pmc-oa-opendata.s3.amazonaws.com/'+k))
for u in re.findall('href="([^"]+)"',(S/'PMC10992486.html').read_text()):
 if '/bin/' in u:jobs.append((u.split('/')[-1],'https://pmc.ncbi.nlm.nih.gov'+u))
for item in dict(jobs).items():get(item)
