from pathlib import Path
import hashlib,json,subprocess,struct
import pandas as pd,numpy as np
from scipy.stats import false_discovery_control
O=Path('results/C5_Endothelial_GSEA_nominalNicheNet_v1');OLD=Path('results/C5_Endothelial_NicheNet_spatial_v1');checks=[]
def ck(name,value,detail=''):checks.append(dict(check=name,status='PASS' if value else 'FAIL',detail=detail))
def r(p):return pd.read_csv(O/p,sep='\t')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
for x in r('00_input_audit/input_hashes.tsv').itertuples():ck('unchanged '+x.path,sha(x.path)==x.sha256)
for x in r('00_input_audit/resource_verification.tsv').itertuples():ck('resource '+x.path,sha(x.path)==x.sha256)
g=r('01_Endothelial_Hallmark_GSEA/all_models.tsv');ck('250 Hallmark tests',len(g)==250)
for m,d in g.groupby('model'):
 ck('GSEA BH '+m,np.allclose(d.padj,false_discovery_control(d.pval,method='bh')));ck('finite '+m,np.isfinite(d[['NES','ES','pval','padj']]).all().all());rr=r('01_Endothelial_Hallmark_GSEA/ranking_'+m+'.tsv');ck('rank descending '+m,(np.diff(rr.t)<=0).all());ck('gene uniqueness '+m,not rr.gene.duplicated().any())
defs=r('02_nominal_response_genesets/frozen_genesets.tsv');a=r('03_nominal_NicheNet/all_activities.tsv');bg=pd.read_csv(OLD/'02_receiver_geneset/expressed_background.tsv',sep='\t');ck('fixed background',len(bg)==7070);ck('genesets background',set(defs.gene)<=set(bg.gene));ck('no duplicated sets',not defs.duplicated().any())
for p,d in a.groupby('program'):ck('all246 '+p,len(d)==246 and d.test_ligand.nunique()==246 and sorted(d['rank'])==list(range(1,247)))
for s,p in [('A',.01),('B',.05),('C',.1)]:
 d=r('02_nominal_response_genesets/Set_'+s+'_metrics.tsv');ck('nominal rules '+s,((d['P.Value']<p)&(d.logFC>0)&(d.logFC_InhouseData>0)&(d.logFC_GSE222315>0)&d.in_background).all());ck('no primary FDR DEG '+s,(d['adj.P.Val']>=.1).all())
lo=r('04_rank_sensitivity/LOPO_status.tsv');ck('21 LOPO fits',len(lo)==21);ck('21 LOPO activities',a.program.str.startswith('LOPO_').sum()==21*246)
rg=r('04_rank_sensitivity/random_genesets.tsv');ra=r('04_rank_sensitivity/random_activities.tsv.gz');ck('50 random programs',ra.program.nunique()==50 and len(ra)==50*246);b=set(defs.loc[defs.program.eq('Set_B'),'gene']);bins=pd.cut(bg.pooled_pct,np.unique(np.quantile(bg.pooled_pct,np.linspace(0,1,11))),include_lowest=True);binmap=dict(zip(bg.gene,bins.astype(str)));exp=pd.Series([binmap[x] for x in b]).value_counts().sort_index()
for p,d in rg.groupby('program'):
 obs=pd.Series([binmap[x] for x in d.gene]).value_counts().sort_index();ck('matched random '+p,len(d)==len(b) and d.gene.nunique()==len(b) and obs.equals(exp))
sp=r('06_CellChat_spatial_integration/new_axes_expression_audit.tsv');ck('63x3 light audit',len(sp)==189 and sp.section_id.nunique()==63);known=r('06_CellChat_spatial_integration/new_axes_known_patient_description.tsv');ck('29 verified spatial patients',known.patient_uid.nunique()==29)
e=r('06_CellChat_spatial_integration/SetB_complete_evidence.tsv');ck('all246 retained integration',e.loc[e['rank'].notna(),'test_ligand'].nunique()==246)
for f in (O/'figures').glob('*.pdf'):
 out=subprocess.run(['pdfimages','-list',str(f)],capture_output=True,text=True,check=True).stdout;ck('vector '+f.name,len(out.strip().splitlines())==2)
for f in (O/'figures').glob('*.png'):
 data=f.read_bytes();i=data.index(b'pHYs');x,y,unit=struct.unpack('>IIB',data[i+4:i+13]);ck('300dpi '+f.name,unit==1 and abs(x*.0254-300)<1 and abs(y*.0254-300)<1)
ck('comprehensive report',(O/'COMPREHENSIVE_REPORT.md').stat().st_size>10000)
pd.DataFrame(checks).to_csv(O/'00_input_audit/final_validation.tsv',sep='\t',index=False)
files=[f for f in O.rglob('*') if f.is_file() and f.name not in ['output_manifest.tsv','08_validate.log']];pd.DataFrame([dict(path=str(f.relative_to(O)),bytes=f.stat().st_size,sha256=sha(f)) for f in files]).to_csv(O/'00_input_audit/output_manifest.tsv',sep='\t',index=False)
failed=[c for c in checks if c['status']=='FAIL'];print('CHECKS',len(checks),'FAIL',len(failed));print(failed);assert not failed
