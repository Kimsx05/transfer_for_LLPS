exec(open('results/fig2_like_C0/phase2_external/scripts/01_download_public.py').read().split('with ThreadPoolExecutor')[0])
urls={'TCGA_BLCA_clinicalMatrix.tsv':'https://tcga.xenahubs.net/download/TCGA.BLCA.sampleMap/BLCA_clinicalMatrix','TCGA_BLCA_survival.tsv':'https://tcga.xenahubs.net/download/survival/BLCA_survival.txt'}
urls.update({f'{g}.annot.gz':f'https://ftp.ncbi.nlm.nih.gov/geo/platforms/{g[:-3]}nnn/{g}/annot/{g}.annot.gz' for g in ['GPL570','GPL6947','GPL6102']})
urls['GPL570.annot.gz']='https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz'
with ThreadPoolExecutor(max_workers=3) as ex: rows=list(ex.map(get,urls.items()))
with (root/'download_manifest_additional.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['file','source_url','bytes','sha256','status']);w.writerows(rows)
