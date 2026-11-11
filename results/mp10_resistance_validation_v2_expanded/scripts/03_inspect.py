from env import *
import tarfile,zipfile
for p in S.glob('*.xlsx'):
 try:
  x=pd.ExcelFile(p);print('\nXLSX',p.name,x.sheet_names)
  for sh in x.sheet_names[:4]:print(sh,pd.read_excel(x,sheet_name=sh,header=None,nrows=5).iloc[:,:10].to_string(index=False,header=False))
 except Exception as e:print('ERROR',p.name,e)
for name in ['GSE231835_RAW.tar','GSE313089_RAW.tar']:
 with tarfile.open(S/name) as t:
  ns=t.getnames();print(name,ns[:15]);f=t.extractfile(ns[0]);d=f.read(3000);print(d[:600])
for name in ['GSE217413_Gene_cpm_WT_CR.xls.gz','GSE268855_counts.csv.gz','GSE190636_T24.GRCh38.CPM.log2.NT.txt.gz']:
 print(name,gzip.open(S/name,'rt').read(1200))
