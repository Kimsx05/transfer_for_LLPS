from pathlib import Path
import urllib.request,gzip,json,re,time,hashlib,csv,concurrent.futures,datetime
P=Path('/home/data/t070721/codex_workspace/Bladder Metabolism');B=P/'results/mp10_resistance_validation_v2_expanded';V=P/'results/mp10_resistance_validation_v1';S=B/'00_plan_audit/sources'
def fetch(item):
 name,url=item;p=S/name
 try:
  if not p.exists():
   with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0'}),timeout=90) as r: data=r.read()
   p.write_bytes(data)
  data=p.read_bytes();r=dict(file=name,url=url,status='OK',bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),time=datetime.datetime.now(datetime.timezone.utc).isoformat())
 except Exception as e:r=dict(file=name,url=url,status=str(e),bytes=0)
 print(name,r['status'],r['bytes'],flush=True);return r
cell=['GSE171537','GSE217413','GSE285648','GSE313089','GSE268855','GSE231835','GSE58624','GSE190636','GSE210954','GSE210921','GSE210953','GSE77883']
clinical=['GSE169455','GSE212810','GSE52219','GSE244266','GSE87304','GSE69795','GSE247185','GSE5287','GSE124305','GSE235067','GSE48276','GSE70691','GSE83586','GSE222073']
if __name__=='__main__':
 jobs=[]
 for g in cell+clinical:
  root=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{g[:-3]}nnn/{g}'
  old=V/'00_audit/sources'/f'{g}_family.soft.gz'
  if old.exists() and not (S/old.name).exists():(S/old.name).symlink_to(old)
  jobs += [(f'{g}_family.soft.gz',root+'/soft/'+g+'_family.soft.gz'),(f'{g}_suppl.html',root+'/suppl/'),(f'{g}_matrix.html',root+'/matrix/')]
 for pre in ['LI','KANG','WHITESIDE']:
  for direct in ['UP','DN']:
   sig=f'{pre}_CISPLATIN_RESISTANCE_{direct}';jobs += [(sig+'.json',f'https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName={sig}&fileType=json'),(sig+'.html',f'https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/{sig}.html')]
 for sig in ['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION','HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT','HALLMARK_P53_PATHWAY','HALLMARK_HYPOXIA','HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY','HALLMARK_DNA_REPAIR','HALLMARK_CHOLESTEROL_HOMEOSTASIS','BENPORATH_ES_CORE_NINE']:
  jobs.append((sig+'.json',f'https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName={sig}&fileType=json'))
 for pmc in ['PMC5039423','PMC11641020','PMC5085609','PMC8437243','PMC10992486']:
  jobs += [(pmc+'.xml',f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML'),(pmc+'_supp.zip',f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/supplementaryFiles')]
 jobs += [('cissig_plot_sig.Rmd','https://raw.githubusercontent.com/jessicascarborough/cissig/c8317dc0a6ab6b8ee66dc98f2e817014f1a3c2a5/Scripts/plot_sig.Rmd'),('CisSig_article.html','https://www.nature.com/articles/s41698-023-00375-y'),('Taber2020.html','https://www.nature.com/articles/s41467-020-18640-0')]
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:res=list(ex.map(fetch,jobs))
 (B/'00_plan_audit/acquisition_initial.json').write_text(json.dumps(res,indent=2))
