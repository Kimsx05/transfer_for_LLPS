from env import *
from pandas.testing import assert_frame_equal
old=read('07_summary/before_official_full/evidence_matrix.tsv.gz');new=read('07_summary/evidence_matrix.tsv.gz')
keys=['module','dataset','signature','endpoint','method','adjustment','subset','branch','background','state','effect_type','family']
vals=['n','effect','CI_low','CI_high','P','FDR']
new=new[~new.method.isin(['official_AUCell','UCell'])]
cols=[c for c in keys+vals if c in old.columns and c in new.columns]
a=old[cols].copy();b=new[cols].copy()
for c in keys:
 if c in cols:a[c]=a[c].fillna('NA').astype(str);b[c]=b[c].fillna('NA').astype(str)
a=a.sort_values(cols,na_position='last').reset_index(drop=True);b=b.sort_values(cols,na_position='last').reset_index(drop=True)
assert_frame_equal(a,b,check_dtype=False,check_exact=False,rtol=1e-10,atol=1e-12)
s=read('00_plan_audit/analysis_status.tsv.gz');assert not s.status.isin(['RUNNING','DEFERRED_RESOURCE']).any()
q=json.loads((B/'00_plan_audit/official_full_QC.json').read_text());assert q['complete'] and q['n_cells']==123699
r=dict(complete=True,original_result_rows_unchanged=len(a),compared_columns=cols,all_preexisting_effects_CIs_P_family_FDR_unchanged=True,official_score_rows=q['n_score_rows'],no_remaining_resource_deferred_or_running_rows=True,remaining_blocked_status_counts=s[s.status.str.startswith('BLOCKED',na=False)].status.value_counts().to_dict(),timestamp=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat())
(B/'00_plan_audit/full_completion_check.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
