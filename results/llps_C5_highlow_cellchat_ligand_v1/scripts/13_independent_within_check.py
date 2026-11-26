import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']:os.environ[k]='1'
import pandas as pd,numpy as np
from scipy.stats import wilcoxon
from pathlib import Path
O=Path('/home/jinshengxi/bladder_metabolism/results/llps_C5_highlow_cellchat_ligand_v1');D=O/'05_LLPS_ligand_correlation';r=pd.read_csv(D/'cell_and_specimen_correlations.csv.gz');su=pd.read_csv(D/'patient_effect_tests.csv.gz');focus=pd.read_csv(O/'00_input_audit/candidate_overlap.csv').gene.tolist();r=r[(r.scope!='pooled_tumor')&(r.clinical_group!='Normal')&(r.n_cells>=5)&(r.sensitivity=='original')&(r.gene.isin(focus))];metrics=['pearson','spearman','adjusted_r','expression_difference','detection_difference'];p=r.groupby(['gene','patient_id'])[metrics].median();checks=[]
for _,z in su[(su.sensitivity=='original')&(su.stratum=='ALL')&(su.gene.isin(focus))].iterrows():
 v=p.loc[z.gene,z.metric].dropna().to_numpy();expected=wilcoxon(v,correction=True,method='approx').pvalue if len(v)>=5 and np.any(v!=0) else 1 if len(v)>=5 else np.nan;ok=len(v)==z.n and np.isclose(np.median(v),z['median'],equal_nan=True) and np.isclose(expected,z.p,equal_nan=True,atol=1e-8);checks.append(dict(gene=z.gene,metric=z.metric,pass_=bool(ok)));assert ok,(z.gene,z.metric,expected,z.p)
pd.DataFrame(checks).to_csv(O/'00_input_audit/independent_within_patient_validation.csv',index=False);print('INDEPENDENT WITHIN PATIENT PASS',len(checks))
