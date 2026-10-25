from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import urllib.request, hashlib, csv
root=Path('results/fig2_like_C0/phase2_external/resources')
urls={f'{g}_series_matrix.txt.gz':f'https://ftp.ncbi.nlm.nih.gov/geo/series/{g[:-3]}nnn/{g}/matrix/{g}_series_matrix.txt.gz' for g in ['GSE31684','GSE48075','GSE13507']}
urls.update({'TCGA_BLCA_HiSeqV2.gz':'https://tcga.xenahubs.net/download/TCGA.BLCA.sampleMap/HiSeqV2.gz','TCGA_BLCA_clinicalMatrix.gz':'https://tcga.xenahubs.net/download/TCGA.BLCA.sampleMap/BLCA_clinicalMatrix.gz','TCGA_BLCA_survival.txt.gz':'https://tcga.xenahubs.net/download/survival/BLCA_survival.txt.gz','TCGA_BLCA_HiSeqV2.json':'https://tcga.xenahubs.net/download/TCGA.BLCA.sampleMap/HiSeqV2.json','MuSiC_construct.R':'https://raw.githubusercontent.com/xuranw/MuSiC/master/R/construct.R','MuSiC_utils.R':'https://raw.githubusercontent.com/xuranw/MuSiC/master/R/utils.R','MuSiC_LICENSE':'https://raw.githubusercontent.com/xuranw/MuSiC/master/LICENSE','MuSiC_commit.json':'https://api.github.com/repos/xuranw/MuSiC/commits/master'})
def get(item):
 name,url=item;p=root/name
 try:
  if not p.exists():
   with urllib.request.urlopen(url,timeout=90) as response, (root/(name+'.part')).open('wb') as f:
    while True:
     b=response.read(1024*1024)
     if not b:break
     f.write(b)
   (root/(name+'.part')).rename(p)
  row=[name,url,p.stat().st_size,hashlib.sha256(p.read_bytes()).hexdigest(),'OK']
 except Exception as e:row=[name,url,0,'',str(e)]
 print(row[0],row[-1],flush=True);return row
with ThreadPoolExecutor(max_workers=3) as ex: rows=list(ex.map(get,urls.items()))
with (root/'download_manifest.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['file','source_url','bytes','sha256','status']);w.writerows(rows)
