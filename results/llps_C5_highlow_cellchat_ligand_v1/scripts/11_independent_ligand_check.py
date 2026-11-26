import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
import numpy as np,pandas as pd
from scipy.stats import pearsonr,spearmanr,t
from pathlib import Path
O=Path('/home/jinshengxi/bladder_metabolism/results/llps_C5_highlow_cellchat_ligand_v1');D=O/'05_LLPS_ligand_correlation'
p=pd.read_csv(D/'patient_pseudobulk_metadata.csv');e=pd.read_csv(D/'pseudobulk_ligand_expression.csv.gz');r=pd.read_csv(D/'pseudobulk_continuous_associations.csv.gz');focus=pd.read_csv(O/'00_input_audit/candidate_overlap.csv').gene.tolist();checks=[]
for _,row in r[(r.gene.isin(focus))&(r.stratum=='ALL')].iterrows():
 d=e[(e.gene==row.gene)&(e.min_cells==row.min_cells)].merge(p,on='patient_id',validate='one_to_one');d=d[d.n_cells>=row.min_cells];sc=d[row.score_type].to_numpy();v=d.logCPM.to_numpy();x=np.column_stack([np.ones(len(d)),sc,pd.get_dummies(d.dataset,drop_first=True,dtype=float).to_numpy()]);beta=np.linalg.lstsq(x,v,rcond=None)[0];df=len(v)-np.linalg.matrix_rank(x);res=v-x@beta;se=np.sqrt(np.sum(res*res)/df*np.linalg.inv(x.T@x).diagonal());pv=2*t.sf(abs(beta[1]/se[1]),df)
 ok=len(d)==row.n_patients and np.isclose(pearsonr(sc,v).statistic,row.pearson,atol=1e-8) and np.isclose(spearmanr(sc,v).statistic,row.spearman,atol=1e-8)
 if np.isfinite(row.adjusted_beta):ok=ok and np.isclose(beta[1],row.adjusted_beta,atol=1e-6) and np.isclose(pv,row.adjusted_p,atol=1e-8)
 checks.append(dict(gene=row.gene,min_cells=row.min_cells,score_type=row.score_type,pass_=bool(ok)))
 assert ok,(row.gene,row.min_cells,row.score_type)
pd.DataFrame(checks).to_csv(O/'00_input_audit/independent_pseudobulk_validation.csv',index=False);print('INDEPENDENT PSEUDOBULK PASS',len(checks))
