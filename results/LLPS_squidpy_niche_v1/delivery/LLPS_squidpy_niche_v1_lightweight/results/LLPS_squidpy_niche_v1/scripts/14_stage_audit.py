from pathlib import Path
from zipfile import ZipFile
import xml.etree.ElementTree as E
import re,hashlib,json
import pandas as pd
O=Path(__file__).resolve().parents[1];meta=pd.read_csv(O/'00_input_audit/section_metadata.tsv',sep='\t');source=Path(meta.source.dropna().iloc[0]);digest=hashlib.sha256(source.read_bytes()).hexdigest();assert set(meta.source_sha256.dropna())=={digest};ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'};rows={}
with ZipFile(source) as z:
 ss=[''.join(t.itertext()) for t in E.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns)]
 for row in E.fromstring(z.read('xl/worksheets/sheet1.xml')).findall('s:sheetData/s:row',ns):
  vals={}
  for c in row.findall('s:c',ns):
   v=c.find('s:v',ns);x=v.text if v is not None else '';vals[re.sub(r'\d','',c.attrib['r'])]=ss[int(x)] if c.get('t')=='s' and x else x
  rows[int(row.attrib['r'])]=vals
assert rows[1]['B']=='Sample' and rows[1]['J']=='T stage'
meta['T_stage_raw']=None;meta['T_stage']=None;meta['T_stage_status']='NO_VERIFIED_STAGE';audit=[]
for i,r in meta.iterrows():
 if r.mapping_rule!='exact_source_accession' or pd.isna(r.workbook_row):continue
 row=rows[int(r.workbook_row)];assert row['B']==r.source_sample_id
 raw=row.get('J','');stage=re.match(r'^(?:p|c)?(T(?:[0-4](?:[ab])?|a|is))(?=[(N\s]|$)',raw)
 meta.loc[i,'T_stage_raw']=raw
 if stage:meta.loc[i,'T_stage']=stage.group(1);meta.loc[i,'T_stage_status']='EXACT_ACCESSION_AND_EXPLICIT_STAGE_TEXT'
 audit.append(dict(section_id=r.section_id,source_sample_id=r.source_sample_id,workbook_row=r.workbook_row,T_stage_raw=raw,T_stage=stage.group(1) if stage else None,source=str(source),sha256=digest))
meta.to_csv(O/'00_input_audit/section_metadata.tsv',sep='\t',index=False);pd.DataFrame(audit).to_csv(O/'00_input_audit/clinical_stage_source_audit.tsv',sep='\t',index=False)
stagecols=['T_stage_raw','T_stage','T_stage_status'];d=pd.read_parquet(O/'04_niche_characterization/spot_metadata.parquet');d=d.drop(columns=[c for c in stagecols if c in d]).merge(meta[['section_id']+stagecols],on='section_id',validate='many_to_one');d.to_parquet(O/'04_niche_characterization/spot_metadata.parquet',index=False,compression='zstd')
valid=meta[meta.T_stage.notna()];valid.groupby(['dataset','T_stage']).agg(n_sections=('section_id','size'),n_verified_patients=('patient_uid','nunique')).reset_index().to_csv(O/'08_clinical_association/stage_sample_sizes.tsv',sep='\t',index=False)
p=pd.read_csv(O/'08_clinical_association/section_niche_proportions_complete.tsv',sep='\t').merge(meta[['section_id']+stagecols],on='section_id',validate='many_to_one');p[p.T_stage.notna()].to_csv(O/'08_clinical_association/stage_niche_proportions_descriptive.tsv',sep='\t',index=False)
report=f'Original workbook checked by SHA256, exact accession and prior verified worksheet row. {len(valid)} sections have explicit T-stage text. Normal rows labeled MIBC-Normal are not assigned tumor stage. GSE171351 has T2 three and T3 one confirmed patients; GSE285715 has T3 two and T4 one sections with unresolved patient mapping. No sufficiently sized within-dataset independent-patient comparison; stage results descriptive only. Earlier preliminary statement of no available T stage is corrected by this source audit. No stage entered clustering.\n';(O/'08_clinical_association/STAGE_AUDIT.md').write_text(report)
manifest=pd.read_csv(O/'00_input_audit/input_manifest.tsv',sep='\t');manifest=manifest[manifest.path!=str(source)];manifest=pd.concat([manifest,pd.DataFrame([dict(path=str(source),role='verified_clinical_source',bytes=source.stat().st_size,sha256=digest)])]);manifest.to_csv(O/'00_input_audit/input_manifest.tsv',sep='\t',index=False)
