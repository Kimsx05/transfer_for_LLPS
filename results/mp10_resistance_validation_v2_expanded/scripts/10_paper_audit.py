from env import *
import urllib.request,time,zipfile
papers={'CrM':'39608442','GSE231835':'10.1093/narcan/zcad054','GSE247185':'39119802','GSE124305':'30224344','GSE69795':'26343003','GSE5287':'17671123','GSE285648':'41216708','GSE313089':'42045965','GSE190636':'36344487','GSE268855':'10.1038/s41416-026-03364-7'}
out=[]
for name,key in papers.items():
 try:
  url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?query='+urllib.request.quote(key)+'&format=json';data=urllib.request.urlopen(url,timeout=25).read();(S/(name+'_paper_lookup.json')).write_bytes(data);j=json.loads(data);r=j['resultList']['result'][0];pmc=r.get('pmcid','');out.append(dict(dataset=name,pmid=r.get('id'),pmcid=pmc,title=r.get('title'),doi=r.get('doi'),isOpenAccess=r.get('isOpenAccess')));print(name,pmc,flush=True)
  if pmc:
   for suffix,url in [('.xml',f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML'),('_s3.xml',f'https://pmc-oa-opendata.s3.amazonaws.com/?list-type=2&prefix={pmc}.')]:
    p=S/(name+suffix)
    if not p.exists():
     try:p.write_bytes(urllib.request.urlopen(url,timeout=25).read())
     except Exception as e:print(name,suffix,str(e),flush=True)
 except Exception as e:out.append(dict(dataset=name,error=str(e)));print(name,str(e),flush=True)
save(pd.DataFrame(out),'00_plan_audit/paper_lookup.tsv')
