from pathlib import Path
import csv,gzip,re,collections,math,json
P=Path('results/fig2_like_C0/phase2_external');R=P/'resources';T=P/'tables'
def write(name,rows):
 rows=list(rows)
 if not rows:return
 fields=list(dict.fromkeys(k for r in rows for k in r))
 with (T/name).open('w') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
def num(x):
 try:return float(x)
 except:return None
def t_stage(s):
 # Preserve ambiguous ranges and mixed Ta/Tis as unclassified, not upgraded/downstaged.
 if not s or re.search(r'Ta[-/]|T[0-4][-/]T',s,re.I):return ''
 m=re.match(r'^(?:c|p)?(T(?:is|a|[0-4]))(?:[ab])?(?:N|$)',s,re.I)
 return m[1] if m else ''
summ=[];issues=[]
for cohort,platform in [('GSE31684','GPL570'),('GSE48075','GPL6947'),('GSE13507','GPL6102')]:
 meta=[];probes=[];nfeatures=0;lo=float('inf');hi=-lo
 with gzip.open(R/f'{cohort}_series_matrix.txt.gz','rt') as f:
  for l in f:
   if l.startswith('!Sample_'):meta.append(next(csv.reader([l],delimiter='\t')))
   if l.startswith('!series_matrix_table_begin'):
    head=next(csv.reader([next(f)],delimiter='\t'));samples=head[1:];break
  for l in f:
   if l.startswith('!'):break
   a=next(csv.reader([l],delimiter='\t'));probes.append(a[0]);nfeatures+=1
   vals=[float(x) for x in a[1:] if x not in ('null','NA','')];lo=min(lo,min(vals));hi=max(hi,max(vals))
 rows=[dict(sample=s,cohort=cohort,platform=platform) for s in samples];raw=[]
 for a in meta:
  assert len(a)==len(samples)+1
  for i,v in enumerate(a[1:]):
   if a[0]=='!Sample_geo_accession':assert v==samples[i]
   if a[0] in ['!Sample_title','!Sample_source_name_ch1','!Sample_organism_ch1','!Sample_data_processing']:
    k=a[0].replace('!Sample_','');rows[i][k]=(rows[i].get(k,'')+' | '+v).strip(' |')
   if a[0]=='!Sample_characteristics_ch1' and ':' in v:
    k,value=v.split(':',1);k=k.strip().lower();value=value.strip()
    if k in rows[i] and rows[i][k]!=value:issues.append(dict(cohort=cohort,sample=samples[i],issue='conflicting characteristic '+k))
    rows[i][k]=value
    raw.append(dict(sample=samples[i],field=k,value=value))
 for d in rows:
  d['patient_unit']=d['title'].split('_')[0] if cohort=='GSE31684' else d['title']
  d['primary_tumor_eligible']=cohort!='GSE13507' or d.get('biological source','')=='primary bladder cancer'
  if cohort=='GSE31684':
   d['stage_raw']=d.get('rc_stage','');d['stage_definition']='pathological T at radical cystectomy';d['grade_raw']=d.get('rc grade','');d['os_months']=num(d.get('survival.months'));d['os_event']={'DOC':1,'DOD':1,'NED':0}.get(d.get('last known status'));d['age_years']=num(d.get('age at rc'))
  elif cohort=='GSE48075':
   d['stage_raw']=d.get('cstage','');d['stage_definition']='pathological T in mixed cstage field' if d['stage_raw'].startswith('p') else 'clinical T';d['grade_raw']='';d['os_months']=num(d.get('survival (mo)'));d['os_event']={'uncensored':1,'censored':0}.get(d.get('os censor'));d['age_years']=num(d.get('age (at specimen collection)'))
  else:
   d['stage_raw']=d.get('stage','');d['stage_definition']='reported TNM T component; clinical/pathological basis unspecified';d['grade_raw']=d.get('grade','');d['os_months']=num(d.get('survival month'));d['os_event']={'death':1,'survival':0}.get(d.get('overall survival'));d['age_years']=num(d.get('age'))
  d['T_component']=t_stage(d['stage_raw']);d['stage_group']='' # T is never converted to overall stage
  d['os_endpoint']='OS';d['source_time_unit']='months';d['os_usable']=d['primary_tumor_eligible'] and d['os_months'] is not None and d['os_months']>0 and d['os_event'] is not None
  if d.get('organism_ch1')!='Homo sapiens' or 'Multipotent Progenitor' in d.get('source_name_ch1',''):issues.append(dict(cohort=cohort,sample=d['sample'],issue='source/organism inconsistent with bladder cohort; exclude from clinical primary-tumor set'))
 write(f'{cohort}_clinical_harmonized.tsv',rows);write(f'{cohort}_characteristics_long.tsv',raw)
 # Use only unambiguously single-symbol annotated probes represented in the downloaded matrix.
 mapping=[];ps=set(probes)
 with gzip.open(R/f'{platform}.annot.gz','rt') as f:
  for l in f:
   if l.startswith('ID\t'):fields=l.rstrip('\n').split('\t');break
  for l in f:
   if l.startswith('!'):break
   a=l.rstrip('\n').split('\t');d=dict(zip(fields,a));sym=d.get('Gene symbol','').strip()
   if d.get('ID') in ps and sym and not re.search(r'///|;|\s',sym):mapping.append(dict(probe=d['ID'],symbol=sym))
 write(f'{platform}_probe_symbol.tsv',mapping)
 (T/f'{platform}_covered_genes.txt').write_text('\n'.join(sorted(set(d['symbol'] for d in mapping)))+'\n')
 eligible=[d for d in rows if d['primary_tumor_eligible']]
 summ.append(dict(cohort=cohort,platform=platform,n_expression_samples=len(samples),n_features=nfeatures,expression_min=lo,expression_max=hi,n_primary_tumor=len(eligible),n_distinct_patient_labels=len(set(d['patient_unit'] for d in eligible)),patient_identity_note='study specimen identifiers; no duplicate labels; explicit cross-specimen patient links absent',n_valid_T=sum(bool(d['T_component']) for d in eligible),n_valid_overall_stage=0,n_valid_OS=sum(d['os_usable'] for d in eligible),OS_events=sum(d['os_event']==1 and d['os_usable'] for d in eligible),time_unit='months',n_mapped_genes=len(set(d['symbol'] for d in mapping))))
# TCGA identifiers and curated Xena survival table; retain both raw clinical and curated endpoint fields.
with gzip.open(R/'TCGA_BLCA_HiSeqV2.gz','rt') as f:samples=next(csv.reader(f,delimiter='\t'))[1:]
with (R/'TCGA_BLCA_clinicalMatrix.tsv').open() as f:clin={d['sampleID']:d for d in csv.DictReader(f,delimiter='\t')}
with (R/'TCGA_BLCA_survival.tsv').open() as f:surv={d['sample']:d for d in csv.DictReader(f,delimiter='\t')}
rows=[]
for s in samples:
 c=clin.get(s,{});v=surv.get(s,{});days=num(v.get('OS.time'));event=num(v.get('OS'));stage=c.get('pathologic_stage','')
 d=dict(sample=s,cohort='TCGA-BLCA',patient_unit=s[:12],primary_tumor_eligible=s[13:15]=='01',age_years=num(c.get('age_at_initial_pathologic_diagnosis')),stage_raw=stage,stage_definition='pathological overall AJCC stage',stage_group=stage if re.match(r'^Stage (I|II|III|IV)$',stage) else '',T_raw=c.get('pathologic_T',''),clinical_T_raw=c.get('clinical_T',''),grade_raw=c.get('neoplasm_histologic_grade',''),os_days=days,os_months=days/30.4375 if days is not None else None,os_event=event,source_time_unit='days',OS_source='Xena curated BLCA_survival.txt; OS=1 death,0 censored',purity=None)
 d['os_usable']=d['primary_tumor_eligible'] and days is not None and days>0 and event in [0,1];rows.append(d)
write('TCGA_clinical_harmonized.tsv',rows)
eligible=[d for d in rows if d['primary_tumor_eligible']]
summ.append(dict(cohort='TCGA-BLCA',platform='Illumina HiSeq RNASeqV2',n_expression_samples=len(samples),n_primary_tumor=len(eligible),n_distinct_patient_labels=len(set(d['patient_unit'] for d in eligible)),patient_identity_note='TCGA 12-character participant barcode',n_valid_overall_stage=sum(bool(d['stage_group']) for d in eligible),n_valid_OS=sum(d['os_usable'] for d in eligible),OS_events=sum(d['os_event']==1 and d['os_usable'] for d in eligible),time_unit='days; months=days/30.4375',expression_type='log2(RSEM normalized counts+1)'))
write('cohort_inventory.tsv',summ);write('metadata_issues.tsv',issues)
for s in summ:print(s)
print('ISSUES',issues)
