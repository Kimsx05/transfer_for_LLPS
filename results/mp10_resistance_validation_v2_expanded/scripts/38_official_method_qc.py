from env import *
os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:4]))
f=read('02_epi/cell_average_tie_scores.tsv.gz')
o=read('02_epi/cell_official_scores.tsv.gz')
sets=read('00_plan_audit/epi_coverage.tsv').query('status=="READY"').signature.tolist()
assert len(f)==123699 and f.cell_id.is_unique
assert set(o.method)=={'official_AUCell','UCell'}
assert len(o)==len(f)*len(sets)*2
assert not o.duplicated(['cell_id','signature','method']).any()
assert set(o.cell_id)==set(f.cell_id) and set(o.signature)==set(sets)
assert np.isfinite(o.score).all() and o.score.between(-1e-12,1+1e-12).all()
meta=f[['cell_id','dataset','sample']].set_index('cell_id')
wide={m:d.pivot(index='cell_id',columns='signature',values='score').reindex(f.cell_id) for m,d in o.groupby('method')}
wide['average_tie']=f.set_index('cell_id')[sets]
rows=[]
for sample,idx in f.groupby('sample',sort=True).indices.items():
 ids=f.iloc[idx].cell_id
 for s in sets:
  for ma,mb in [('average_tie','official_AUCell'),('average_tie','UCell'),('official_AUCell','UCell')]:
   a=wide[ma].loc[ids,s].to_numpy();b=wide[mb].loc[ids,s].to_numpy();valid=np.isfinite(a)&np.isfinite(b)
   rho=stats.spearmanr(a[valid],b[valid]).statistic if valid.sum()>=3 and np.ptp(a[valid])>0 and np.ptp(b[valid])>0 else np.nan
   rows.append(dict(dataset=f.iloc[idx[0]].dataset,sample=sample,signature=s,method_a=ma,method_b=mb,n_cells=len(idx),Spearman=rho,mean_a=np.mean(a),mean_b=np.mean(b),SD_a=np.std(a),SD_b=np.std(b),reason='' if np.isfinite(rho) else 'insufficient cells or constant score; no inferential P',interpretation='within-sample descriptive method agreement; absolute score scales not interchangeable'))
save(pd.DataFrame(rows),'02_epi/method_score_concordance.tsv')
qcrows=[]
cv=read('00_plan_audit/epi_coverage.tsv').set_index('signature')
for ds,idx in f.groupby('dataset',sort=True).indices.items():
 ids=f.iloc[idx].cell_id
 for method,w in wide.items():
  for sn in sets:
   v=w.loc[ids,sn].to_numpy()
   qcrows.append(dict(module='scoring_QC',dataset=ds,signature=sn,endpoint='score',method=method,n_cells=len(v),mean=float(np.mean(v)),SD=float(np.std(v)),n_zero=int(sum(v==0)),coverage=cv.loc[sn,'coverage'],n_present=cv.loc[sn,'n_present'],n_original=cv.loc[sn,'n_original'],status='COMPLETE',P=np.nan,reason='complete finite full-background scoring; cell count is not patient sample size'))
save(pd.DataFrame(qcrows),'02_epi/all_methods_score_QC.tsv')
cy=read('02_epi/cell_cycle_scores.tsv.gz').set_index('cell_id');checks=[]
for p in sorted((B/'02_epi/specimen_checkpoints').glob('*/cycle_scores.tsv')):
 d=pd.read_csv(p,sep='\t').set_index('cell_id');a=cy.loc[d.index,['S.Score','G2M.Score']];diff=float(np.max(np.abs(a.to_numpy()-d[['S.Score','G2M.Score']].to_numpy())));assert diff<1e-8
 checks.append(dict(sample=p.parent.name,n_cells=len(d),max_absolute_cycle_difference=diff,passed=True))
assert sum(d['n_cells'] for d in checks)==len(f)
save(pd.DataFrame(checks),'00_plan_audit/official_cycle_recomputation_QC.tsv')
(B/'00_plan_audit/official_full_QC.json').write_text(json.dumps(dict(complete=True,n_cells=len(f),n_samples=len(checks),n_sets=len(sets),n_methods=2,n_score_rows=len(o),all_IDs_unique_and_complete=True,finite_and_unit_interval=True,max_cycle_recomputation_error=max(d['max_absolute_cycle_difference'] for d in checks),date=datetime.datetime.now().isoformat()),indent=2))
import shutil
shutil.copy2(B/'02_epi/specimen_checkpoints/CNP0000460__CNP0000460_P01N/sessionInfo.txt',B/'00_plan_audit/R_epi_full_sessionInfo.txt')
if (B/'logs/resource_snapshots.jsonl').exists():shutil.copy2(B/'logs/resource_snapshots.jsonl',B/'00_plan_audit/full_execution_resource_snapshots.jsonl')
job=json.loads((B/'00_plan_audit/night_job.json').read_text());execution=dict(cpu_affinity=sorted(os.sched_getaffinity(0)),R_workers=job.get('R_workers'),compute_CPU_budget=4,thread_env={k:os.environ.get(k) for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','RCPP_PARALLEL_NUM_THREADS','R_PROFILE','R_PROFILE_USER']},authorization=job.get('authorization_record'))
assert len(execution['cpu_affinity'])<=4
(B/'00_plan_audit/full_execution_resource_verification.json').write_text(json.dumps(execution,indent=2))
print('OFFICIAL_FULL_QC_PASS',len(o),'rows',len(f),'cells',len(sets),'sets',flush=True)
