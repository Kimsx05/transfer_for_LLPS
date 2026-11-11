from env import *
import importlib
# Clinical raw metadata and independently evaluable endpoints.
def endpoint(t,n,surgery=True):
 if not surgery:return np.nan,np.nan
 t=str(t).upper();n=str(n).upper();nz=n in ['0','N0'];np_=n in ['1','2','3','N1','N2','N3'];tz=t in ['0','T0'];tl=t in ['A','IS','1','TA','TIS','T1','NONINVASIVE_OR_CIS'];th=bool(re.match(r'^(T)?[234]',t))
 return (0 if np_ or th else 1 if (tz or tl) and nz else np.nan,0 if np_ or th or tl else 1 if tz and nz else np.nan)
def tn(s):
 s=str(s).upper().replace('NOMO','N0M0').replace('NOM0','N0M0');t=re.search(r'T(0|IS|A|[1234][AB]?|X)',s);n=re.search(r'N([0123X])',s);return(t[1] if t else '',n[1] if n else '')
cc=[]
v=pd.read_csv(V/'03_clinical/patient_inclusion_endpoint_audit.tsv',sep='\t')
for _,r in v.iterrows():
 x=dict(dataset=r.cohort,matrix_id=r.title if r.cohort=='GSE169455' else r.GSM,GSM=r.GSM,patient_id=r.patient_id,pretreatment=True,regimen=r.regimen,cycles=r.cycles,cT=r.cT,subtype=r.subtype,surgery=r.surgery,ypT=r.ypT,ypN=r.ypN,pathology_raw=r.pathology_raw,good_response_ypT_le1_N0=r.good,pCR_ypT0N0=r.pCR,author_response=r.author_good,strict_inclusion=r.strict_inclusion,fallback_inclusion=r.fallback_inclusion,therapy_intent=r.therapy_intent,evidence='v1 audited original GEO/S6A; encoding retained and checked',raw_ypT_code=r.ypT_code,raw_ypN_code=r.ypN_code,procedure=r.procedure)
 cc.append(x)
for g in ['GSE52219','GSE244266','GSE87304','GSE69795','GSE247185','GSE5287','GSE124305','GSE235067']:
 d=read('00_plan_audit/'+g+'_samples.tsv')
 for _,r in d.iterrows():
  x=dict(dataset=g,matrix_id=r.GSM,GSM=r.GSM,patient_id=r.title,pretreatment=True,regimen='unknown',cycles=np.nan,cT=np.nan,subtype=np.nan,surgery=np.nan,ypT='',ypN='',pathology_raw=r.get('pstage',''),good_response_ypT_le1_N0=np.nan,pCR_ypT0N0=np.nan,author_response=np.nan,strict_inclusion=False,fallback_inclusion=True,therapy_intent='NAC',evidence='GEO full SOFT original fields')
  if g=='GSE52219':
   x.update(regimen='MVAC',cT=tn(r.cstage)[0],subtype=r['subset type'],surgery=True,strict_inclusion=True,author_response=1 if r['response to mvac (downstaged to p0 or t1 w/ nac)']=='yes' else 0);x['ypT'],x['ypN']=tn(r.pstage);x['good_response_ypT_le1_N0'],x['pCR_ypT0N0']=endpoint(x['ypT'],x['ypN'])
  elif g=='GSE244266':x.update(patient_id=re.sub('_F2$','',r.title),regimen=r['treatment arm'],cT='2' if r['clinical_stage_strat_factor']=='Clinical Stage -T2' else '3/4',strict_inclusion=True,evidence='GEO189 assays; duplicate base titles median-merge independent of scores/outcomes; ypTN missing')
  elif g=='GSE87304':x.update(regimen='cisplatin_based_NAC_author_cohort_individual_unknown',cT=r.cstage,subtype=r.bca_subtype,evidence='GEO patient305 cohort; no verified per-patient response crosswalk')
  elif g=='GSE69795':
   pre=r['sample procured']=='TUR';surgery=r.get('cystectomy')=='Yes';blca=str(r.get('source_name_ch1')).upper()=='BLADDER' and r.get('cystectomy')!='Upper tract';ct=pd.to_numeric(r.get('cstage'),errors='coerce');p=str(r.get('pstage'));author=1 if p in ['0','1','TisN0'] else 0 if p in ['2','3','4','4*'] else np.nan
   x.update(pretreatment=pre,regimen='ddMVAC_bevacizumab',cT=ct,surgery=surgery,strict_inclusion=pre and surgery and blca and ct>=2,fallback_inclusion=pre,author_response=author,therapy_intent='NAC_plus_antiVEGF',subtype=r.get('subtype (tur) 1-basal;2-p53 like;3-luminal'),pair_pre_raw=r.get('pretreatment sample record-pair'),evidence='GEO stage GROUP is author endpoint only; numeric pstage never converted to ypT; exact TisN0 retained')
   if p=='TisN0' and surgery:x.update(ypT='IS',ypN='0',good_response_ypT_le1_N0=1,pCR_ypT0N0=0)
  elif g=='GSE247185':x.update(matrix_id=r.description,patient_id='Patient'+re.search(r'number (\d+)',r.title)[1],regimen=r['treatment.'],author_response=1 if r['response to_tratment']=='pCR' else 0,evidence='GEO author pCR includes no muscle-invasive disease and biopsy evaluations; not strict ypT0N0; Gem footnote gentamicin preserved',strict_inclusion=True)
  elif g=='GSE5287':x.update(therapy_intent='advanced_first_line',regimen='cisplatin_containing_author_cohort',evidence='GEO survival only; actual radiologic response and patient-linked treatment timing not recovered')
  elif g=='GSE124305':x.update(patient_id=r.description,pretreatment=False,therapy_intent='post_NAC_state',subtype=r.cluster,fallback_inclusion=False,evidence='post-NAC RC/scar samples; no verified pre-NAC patient crosswalk yet')
  cc.append(x)
# Taber original clinical table; fresh-frozen expression96, not all300/121.
t=pd.read_csv(S/'excel_tsv/41467_2020_18640_MOESM4_ESM.xlsx__1.tsv',sep='\t').iloc[4:].reset_index(drop=True);exprids=read('00_plan_audit/prepared/Taber2020_samples.tsv').matrix_id
for _,r in t.iterrows():
 vals=r.tolist();pid=vals[0]
 if pid not in set(exprids):continue
 nac=bool(pd.notna(vals[11]) and vals[11]!='-');fl=bool(pd.notna(vals[14]) and vals[14]!='-');surg=str(vals[7]).upper().split('|')[0]=='YES';pt,pn=tn(str(vals[9])+str(vals[10]));pt='IS' if 'CIS' in str(vals[9]).upper() else pt;good,pcr=endpoint(pt,pn,surg);author=1 if vals[16] in ['CR','PR'] else 0 if vals[16] in ['SD','PD'] else np.nan
 cc.append(dict(dataset='Taber2020',matrix_id=pid,GSM='',patient_id=pid,pretreatment=True,regimen=vals[11] if nac else vals[14],cycles=vals[12] if nac else vals[15],cT=tn(vals[4])[0],subtype=np.nan,surgery=surg,ypT=pt,ypN=pn,pathology_raw=str(vals[9])+' '+str(vals[10]),good_response_ypT_le1_N0=good if nac else np.nan,pCR_ypT0N0=pcr if nac else np.nan,author_response=author if not nac else (1 if str(vals[13]).lower()=='yes' else 0 if str(vals[13]).lower()=='no' else np.nan),strict_inclusion=nac and str(vals[11]) in ['GC','MVAC','ddMVAC'],fallback_inclusion=True,therapy_intent='NAC_and_FL' if nac and fl else 'NAC' if nac else 'advanced_first_line',evidence='Taber2020 public SupplementaryData1+4 joined NewID; author response cohort; individual biopsy dates absent, contemporaneous pre-first-line timing not independently verified; 96 fresh frozen',raw_NAC_regimen=vals[11],raw_FL_regimen=vals[14],raw_FL_response=vals[16],raw_surgery=vals[7]))
cl=pd.DataFrame(cc);save(cl,'05_clinical/patient_endpoint_audit.tsv')
# Uniform model identity audit; states are not independently induced RNA replicates.
models=[]
for g in importlib.import_module('00_acquire').cell:
 if g=='GSE210954':continue # parent accession only;31unique GSM supplied in children
 p=B/'00_plan_audit/prepared'/f'{g}_samples.tsv'
 if not p.exists():continue
 d=pd.read_csv(p,sep='\t')
 for _,r in d.iterrows():
  title=r.get('title',r.get('GEO_title',''));line=r.get('cell_line',r.get('cell line',''));state='';drug='cisplatin';include=True;exposure='not_documented';process='single derivative; RNA/culture replicates, independent inductions unproven'
  if g=='GSE171537':line=r.cell_line;state=r.state;drug={'CisR':'cisplatin','GemR':'gemcitabine','GemCisR':'gemcitabine_cisplatin','parental':'none'}[state];exposure='several drug-free passages; sampled untreated (GSM protocol)'
  elif g=='GSE217413':state='CisR' if 'resistant' in title else 'parental';line='UM-UC-3';exposure='maintained in CDDP5uM; RNA washout not explicitly documented; human Ensembl and hg19 agree; mm9 protocol text contradictory'
  elif g=='GSE285648':state='CisR' if 'LTT' in str(line) else 'parental';line=str(line).replace('LTT','');include=r.treatment=='DMSO';exposure='DMSO72h baseline arm; long-term cisplatin background; PARPi arms separately scored'
  elif g=='GSE313089':state='CisR' if 'N/P' in str(line) else 'parental';line='NTUB1';exposure='untreated RNA; historical induction separate'
  elif g=='GSE268855':state='CisR' if 'RT4R' in title else 'parental';line='RT4';exposure='GEO Cisplatin treatment for RT4R vs PBS for RT4; exposure confounded with resistance'
  elif g=='GSE231835':state='CisR' if 'T24R2' in title else 'parental';line='T24';exposure='PMC10653028: T24R2 maintained6.6uM cisplatin; RNA3days confluence; RNA-specific washout not stated'
  elif g=='GSE58624':state='CisR' if 'PR' in title else 'parental';line='5637' if '5637' in title else 'T24';process='one expression assay/state; no within-state variance inference'
  elif g=='GSE190636':m=re.search('P(3|7|15)',title);state='GemR_P'+m[1] if m else 'parental';line='T24';drug='gemcitabine';exposure='NT (not treated) sampling; P denotes phase, not days'
  elif g in ['GSE210921','GSE210953']:m=re.search(r'P(\d+)',title);state='parental' if m[1]=='0' else 'GemR_P'+m[1];line='5637';drug='gemcitabine';process='GRC1/2 establishment processes; near-identical same-name P0; physical RNA identity not proven; treat controls as dependent'
  elif g=='GSE77883':state='GemR' if 'resistant' in title else 'parental';line='T24';drug='gemcitabine'
  models.append(dict(dataset=g,study='GSE210954' if g in ['GSE210921','GSE210953'] else g,GSM=r.GSM,matrix_id=r.matrix_id,title=title,background=line,state=state,drug=drug,baseline_include=include,RNA_replicate=title,derived_model=g+'|'+str(line)+'|'+state,replicate_type=process,exposure=exposure))
mo=pd.DataFrame(models);save(mo,'03_cell_models/model_identity.tsv')
# Apply the explicit user confirmation only to the exact12 Inhouse samples from the same123699-cell source.
e=pd.read_csv(V/'config/frozen_v1.json') if False else pd.read_csv(P/'mp_ucell_state_v1/00_input_audit/sample_patient_mapping.tsv',sep='\t');ih=e[e.dataset=='InhouseData'];assert len(ih)==12 and set(['InhouseData|T4463209a','InhouseData|T4463209b']).issubset(set(ih['sample']))
e['patient']=e['patient'].astype('object');e['verified']=False;e['evidence_v2']='unverified local rename; source identity correspondence still needed';mask=e.dataset=='InhouseData';e.loc[mask,'patient']=e.loc[mask,'sample'].str.replace('|','::',regex=False).str.replace(r'T4463209[ab]$','T4463209',regex=True);e.loc[mask,'verified']=True;e.loc[mask,'evidence_v2']='explicit user confirmation in v2 task; exact12 sample IDs match v1/source cohort';save(e,'00_plan_audit/epi_patient_identity.tsv')
# Explicit duplicate accessions, sharing and training relationships.
dup=[]
for g in ['GSE169455','GSE69795','GSE70691','GSE222073']:
 d=read('00_plan_audit/'+g+'_samples.tsv')
 for _,r in d.iterrows():
  refs=re.findall(r'GSM\d+',str(r.get('relation',''))+' '+str(r.get('overlap with previous submisson gse83586',''))+' '+str(r.title))
  for ref in set(refs)-{r.GSM}:dup.append(dict(dataset=g,GSM=r.GSM,related_GSM=ref,relationship='source/reanalysis reference',evidence='GEO original metadata'))
save(pd.DataFrame(dup),'00_plan_audit/sample_deduplication.tsv');save(pd.DataFrame([dict(source='GSE210954',target='GSE210921;GSE210953',relationship='parent-child same31GSM; counted once'),dict(source='GSE52219',target='GSE70691',relationship='reanalysis; not independent NAC validation'),dict(source='GSE87304',target='GSE124305',relationship='pre/post study; patient crosswalk required'),dict(source='CisSig19',target='GDSC_BLCA',relationship='training source overlap; internal consistency'),dict(source='NPEPPS_VRAC_ratio',target='GSE171537',relationship='related model study system'),dict(source='Tanaka49',target='GSE58624',relationship='signature derivation source'),dict(source='CrM63_weighted_z',target='GSE210921;GSE210953',relationship='signature derivation source'),dict(source='CLSPN_derived_program',target='GSE217413',relationship='derived on same data; internal only')]),'00_plan_audit/independence_matrix.tsv')
print('DESIGN',len(cl),'clinical assay rows',len(mo),'cell-model assays','verified Epi',e.loc[e.verified,'patient'].nunique())
