from env import *
import importlib,tarfile,zipfile,io
mod=importlib.import_module('06_prepare');annotated_matrix=mod.annotated_matrix
# Finish v1/Taber data pathway (all previous datasets preserved).
g='GSE212810';aa=[]
with tarfile.open(V/'00_audit/sources/GSE212810_RAW.tar') as t:
 for name in t.getnames():
  z=pd.read_csv(gzip.open(t.extractfile(name),'rt'),sep='\t',index_col=0);aa.append(z.iloc[:,0].rename(name.split('_')[0]))
m=read('00_plan_audit/'+g+'_samples.tsv');m['matrix_id']=m.GSM;annotated_matrix(g,pd.concat(aa,axis=1),m,'counts')
a=mod.excel('41467_2020_18640_MOESM7_ESM.xls__1.tsv');a=a.set_index(a.columns[0]);m=pd.DataFrame({'matrix_id':a.columns,'dataset':'Taber2020','patient_id':a.columns});annotated_matrix('Taber2020',a,m,'author_log2CPM',note='Public Supplementary Data4; fresh frozen96 of RNA121; edgeR normalized expression')
# Independent baseline models: real GDSC2/1 fits joined by SIDM, not by guessed names.
models=pd.read_csv(S/'CMP_models.csv');dr=[]
for screen in ['GDSC1','GDSC2']:
 d=mod.excel(screen+'_fitted_dose_response_27Oct23.xlsx__1.tsv');d=d[(d.TCGA_DESC=='BLCA')&(d.DRUG_NAME.str.lower()=='cisplatin')].copy();dr.append(d)
dr=pd.concat(dr);save(dr,'04_measured_drug_response/GDSC_cisplatin_measured.tsv');ids=sorted(dr.SANGER_MODEL_ID.unique())
with zipfile.ZipFile(S/'CMP_RNAseq.zip') as z:
 with z.open('rnaseq_tpm_20220624.csv') as f:
  header=f.readline().decode().strip().split(',');cols=[i for i,h in enumerate(header) if h in ids or i<2]
 with z.open('rnaseq_tpm_20220624.csv') as f:a=pd.read_csv(f,skiprows=[1,2,3,4],usecols=cols)
sym=a.iloc[:,1];a=a.set_index(a.columns[0]).iloc[:,1:];m=models[models.model_id.isin(a.columns)].copy();m['matrix_id']=m.model_id
annotated_matrix('GDSC_BLCA',a,m,'linear_TPM',symbol=sym,note='CMP RNASeq20220624 TPM; stable SIDM linked to real release8.5 cisplatin fits; CisSig-GDSC training overlap')
# Explicit limited panel: preserve endogenous genes, no full-transcriptome scoring.
a=mod.excel('GSE235067_Raw_data_AFont_et_al_FINAL_v3.xlsx__1.tsv');a=a[a['Class Name']=='Endogenous'];xx=a.set_index('Gene Name').iloc[:,2:];m=read('00_plan_audit/GSE235067_samples.tsv');m['matrix_id']=m.title;annotated_matrix('GSE235067',xx,m,'panel_counts',note='41-gene NanoString endogenous panel; no ssGSEA rank scoring')
