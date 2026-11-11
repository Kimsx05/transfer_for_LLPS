from env import *
import io,tarfile,zipfile
D=B/'00_plan_audit/prepared';D.mkdir(exist_ok=True)
HG=pd.read_csv(S/'hgnc_complete_set.txt',sep='\t',dtype=str).fillna('');approved=set(HG.symbol)
def index_multicol(cols):
 d={}
 for _,r in HG.iterrows():
  for c in cols:
   for k in r.get(c,'').split('|'):
    if k:d.setdefault(k,set()).add(r.symbol)
 return d
AL=index_multicol(['alias_symbol','prev_symbol']);EN=index_multicol(['ensembl_gene_id']);ENT=index_multicol(['entrez_id'])
def mapgene(original,symbol):
 original=str(original);symbol='' if pd.isna(symbol) else str(symbol);base=original.split('__')[0].split('.')[0]
 if symbol in approved:return symbol,'HGNC_exact'
 if base in EN and len(EN[base])==1:return next(iter(EN[base])),'HGNC_Ensembl'
 if base in ENT and len(ENT[base])==1:return next(iter(ENT[base])),'HGNC_Entrez'
 if symbol not in approved and len(AL.get(symbol,[]))==1:return next(iter(AL[symbol])),'HGNC_unique_previous_or_alias'
 return 'UNMAPPED|'+original,'ambiguous' if len(AL.get(symbol,[]))>1 else 'unmapped_unique_original'
def annotated_matrix(g,df,meta,scale,original=None,symbol=None,note=''):
 original=df.index.astype(str).tolist() if original is None else list(original);symbol=original if symbol is None else list(symbol)
 mapped=[mapgene(o,s) for o,s in zip(original,symbol)];mapping=pd.DataFrame(dict(original_id=original,author_symbol=symbol,mapped_symbol=[m[0] for m in mapped],mapping=[m[1] for m in mapped]));save(mapping,'00_plan_audit/'+g+'_gene_mapping.tsv')
 xx=df.astype(float).copy();assert len(xx)==len(mapping);xx.index=mapping.mapped_symbol
 # Only explicitly mapped gene entities are aggregated. Unmapped original IDs remain distinct.
 if xx.index.duplicated().any():xx=xx.groupby(level=0,sort=False).sum() if scale=='counts' else xx.groupby(level=0,sort=False).median()
 x=xx.reset_index(names='gene');save(x,'00_plan_audit/prepared/'+g+'_expression.tsv.gz');save(meta,'00_plan_audit/prepared/'+g+'_samples.tsv')
 js=dict(dataset=g,scale=scale,n_features_original=len(df),n_background=len(xx),n_assays=df.shape[1],missing_values=int(xx.isna().sum().sum()),min=float(xx.min().min()),median=float(np.nanmedian(xx)),max=float(xx.max().max()),note=note)
 (D/(g+'_input.json')).write_text(json.dumps(js,indent=2));print('PREPARED',js,flush=True)
def excel(name,header=0):return pd.read_csv(S/'excel_tsv'/name,sep='\t',skiprows=1+header)
def array(g):
 p=S/(g+'_series_matrix.txt.gz');lines=gzip.open(p,'rt');buf=[];on=False
 for l in lines:
  if l.startswith('!series_matrix_table_begin'):on=True
  elif l.startswith('!series_matrix_table_end'):break
  elif on:buf.append(l)
 a=pd.read_csv(io.StringIO(''.join(buf)),sep='\t',index_col=0);m=read('00_plan_audit/'+g+'_samples.tsv');platform=m.platform_id.iloc[0];anno=read('00_plan_audit/'+platform+'_annotation.tsv');anno['ID']=anno['ID'].astype(str);anno=anno.set_index('ID').reindex(a.index.astype(str));sym=None
 for c in ['Gene Symbol','GENE_SYMBOL','Symbol','ORF']:
  if c in anno:sym=anno[c];break
 if 'gene_assignment' in anno:
  sym=anno.gene_assignment.fillna('').map(lambda s:';'.join(sorted(set(x.split(' // ')[1] for x in s.split(' /// ') if len(x.split(' // '))>1))))
 sym=sym.fillna('');sym=sym.map(lambda s: '' if '///' in str(s) or ';' in str(s) else s)
 if platform=='GPL22995':sym=pd.Series(a.index,index=a.index)
 m['matrix_id']=m.GSM
 q=np.nanquantile(a.values,[.01,.5,.99]);scale='author_log2' if q[2]<30 else 'linear_normalized'
 # SCAN z-like expression has no absolute linear abundance; ratio specifically blocked later.
 if platform=='GPL22995':scale='SCAN'
 if g=='GSE77883':scale='author_normalized_unknown_linear_basis'
 annotated_matrix(g,a,m,scale,symbol=sym.to_list(),note=m.data_processing.iloc[0])
if __name__=='__main__':
 for g in ['GSE58624','GSE77883','GSE52219','GSE244266','GSE87304','GSE69795','GSE5287','GSE124305']:
  if not (D/(g+'_input.json')).exists():array(g)
 for g,fn,skip,scale in [('GSE285648','GSE285648_Experiment_processed_data.xlsx__1.tsv',0,'linear_RPKM'),('GSE217413','GSE217413_Gene_cpm_WT_CR.xls__1.tsv',0,'linear_CPM'),('GSE210921','GSE210921_5637GRC1.xlsx__1.tsv',0,'author_log2CPM'),('GSE210953','GSE210953_5637GRC2.xlsx__1.tsv',0,'author_log2CPM'),('GSE247185','GSE247185_BladderProject_RawData.xlsx__1.tsv',1,'counts')]:
  if (D/(g+'_input.json')).exists():continue
  a=excel(fn,skip).set_index(excel(fn,skip).columns[0]);m=read('00_plan_audit/'+g+'_samples.tsv')
  if g=='GSE217413':m['matrix_id']=m.title.map(lambda s:'UM.UC3.CisR_CR'+s[-1] if 'resistant' in s else 'UM.UC3.WT_WT'+s[-1]);sy=[str(x).split('__')[-1] for x in a.index]
  elif g=='GSE285648':
   m['matrix_id']=m.description.str.extract(r'(1016-\d+)')[0]+'_RPKM';sy=a.index
  elif g=='GSE247185':m['matrix_id']=m.description;sy=a.index
  else:m['matrix_id']=m.title.str.replace(r' GRC[12]$','',regex=True);sy=a.index
  assert set(m.matrix_id)==set(a.columns),(g,m.matrix_id.tolist(),a.columns.tolist())
  annotated_matrix(g,a,m,scale,symbol=sy)
 for g in ['GSE231835','GSE313089']:
  if (D/(g+'_input.json')).exists():continue
  aa=[];sy=None
  with tarfile.open(S/(g+'_RAW.tar')) as t:
   for name in t.getnames():
    z=pd.read_csv(gzip.open(t.extractfile(name),'rt'),sep='\t',comment='#');gsm=name.split('_')[0];aa.append(z.set_index(z.columns[0]).iloc[:,-1].rename(gsm))
    if 'gene_name' in z:sy=z.gene_name.to_list()
  a=pd.concat(aa,axis=1);m=read('00_plan_audit/'+g+'_samples.tsv');m['matrix_id']=m.GSM;annotated_matrix(g,a,m,'counts',symbol=sy)
 for g,fn,scale in [('GSE268855','GSE268855_counts.csv.gz','counts'),('GSE190636','GSE190636_T24.GRCh38.CPM.log2.NT.txt.gz','author_log2CPM')]:
  if (D/(g+'_input.json')).exists():continue
  a=pd.read_csv(S/fn,sep=',' if 'csv' in fn else '\t',index_col=0);m=read('00_plan_audit/'+g+'_samples.tsv');print(g,a.columns.tolist(),flush=True)
  if g=='GSE268855':m=m[m.GSM.isin(['GSM919284'+str(i) for i in range(6)])].copy();m['matrix_id']=m.title.str.replace('RT4R-','R',regex=False).str.replace('RT4-','NC',regex=False)
  else:m['matrix_id']=m.title.str.replace('-','.',regex=False)
  assert set(m.matrix_id)==set(a.columns),(g,m.matrix_id.to_list(),a.columns.tolist());annotated_matrix(g,a,m,scale)
 # Reuse v1 public inputs without altering them.
 g='GSE171537';a=pd.read_csv(V/'00_audit/sources/GSE171537_RNAseq_ensembl_raw_counts.txt.gz',sep='\t');m=pd.read_csv(V/'02_cell_lines/sample_mapping.tsv',sep='\t');annotated_matrix(g,a.iloc[:,3:].set_axis(a.iloc[:,0]),m,'counts',symbol=a.iloc[:,1])
 g='GSE169455';a=pd.read_csv(V/'00_audit/sources/GSE169455_normalized_by_gene.txt.gz',sep='\t',skiprows=6,index_col=0);m=read('00_plan_audit/'+g+'_samples.tsv');m['matrix_id']=m.title;annotated_matrix(g,a,m,'author_log2',note='author RMA/quantile/ComBat; probe median; no repeat transformation')
 g='GSE212810';aa=[]
 with tarfile.open(V/'00_audit/sources/GSE212810_RAW.tar') as t:
  for name in t.getnames():
   z=pd.read_csv(gzip.open(t.extractfile(name),'rt'),sep='\t',header=0,index_col=0);aa.append(z.iloc[:,0].rename(name.split('_')[0]))
 m=read('00_plan_audit/'+g+'_samples.tsv');m['matrix_id']=m.GSM;annotated_matrix(g,pd.concat(aa,axis=1),m,'counts')
 # Public Taber normalized expression; labels are joined separately, not inferred from order.
 a=excel('41467_2020_18640_MOESM7_ESM.xls__1.tsv');a=a.set_index(a.columns[0]);m=pd.DataFrame({'matrix_id':a.columns,'dataset':'Taber2020','patient_id':a.columns});annotated_matrix('Taber2020',a,m,'author_log2CPM',note='Supplementary Data4 public 96-patient normalized expression; original study RNA121 distinct scope')
