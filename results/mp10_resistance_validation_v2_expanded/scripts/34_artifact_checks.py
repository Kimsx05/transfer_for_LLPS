from env import *
import zipfile
from pypdf import PdfReader
from PIL import Image
checks=[]
def ck(name,passed,detail=''):
 checks.append(dict(check=name,passed=bool(passed),detail=detail));assert passed,(name,detail)
fig=read('07_summary/figure_index.tsv');ck('figure_index_nonempty',len(fig)>0,str(len(fig)))
for _,r in fig.iterrows():
 stem=B/'07_summary/figures'/r.figure;pdf=Path(str(stem)+'.pdf');png=Path(str(stem)+'.png');ck(r.figure+'_PDF',pdf.exists() and len(PdfReader(pdf).pages)>=1)
 im=Image.open(png);dpi=im.info.get('dpi',(0,0));ck(r.figure+'_PNG300',all(abs(v-300)<1 for v in dpi),str(dpi))
 for suffix in ['_legend_only.pdf','_panel_only.pdf']:
  p=Path(str(stem)+suffix)
  if p.exists():ck(r.figure+suffix,len(PdfReader(p).pages)>=1)
report=B/'07_summary/中文结果报告.html';text_=report.read_text();links=re.findall(r'(?:href|src)="([^\"]+)"',text_)
for link in links:
 if link.startswith(('http:','https:','#')):continue
 ck('report_link', (report.parent/link).exists(),link)
pack=B/'mp10_resistance_validation_v2_expanded_review.zip'
with zipfile.ZipFile(pack) as z:
 ck('ZIP_CRC',z.testzip() is None);ck('no_large_raw_matrix',not any('/prepared/' in n or '/sources/' in n or '/runtime/' in n for n in z.namelist()));ck('self_contained_report', '07_summary/中文结果报告.html' in z.namelist());ck('review_size_small',pack.stat().st_size<25*1024**2,str(pack.stat().st_size))
 if '00_plan_audit/amendments/GSE285648_DMSO_reference.json' not in z.namelist():print('AMENDMENT missing in zip; packaging should recurse audit amendments')
q=read('05_clinical/MP10_encoding_flip_check.tsv');ck('all25_outcome_flips',q.passed.all(),str(len(q)));e=read('07_summary/evidence_matrix.tsv.gz');ck('P_range',e.P.dropna().between(0,1).all());ck('FDR_range',e.FDR.dropna().between(0,1).all());ck('all_signatures_in_evidence',set(read('01_signatures/signature_manifest.tsv').signature)<=set(e.signature));ck('no_zero_for_blocked_P',e.loc[e.status.str.startswith(('BLOCKED','DEFERRED'),na=False),'P'].isna().all())
p=B/'02_epi/cell_official_scores.tsv.gz'
if p.exists():
 d=pd.read_csv(p,sep='\t');nsets=read('00_plan_audit/epi_coverage.tsv').status.eq('READY').sum();ck('official_full_count',len(d)==123699*nsets*2,str(len(d)));ck('official_unique_cells',d.cell_id.nunique()==123699);ck('official_score_finite',np.isfinite(d.score).all())
night=json.loads((B/'logs/night_status.json').read_text());ck('real_night_supervisor',night.get('supervisor_pid') is not None,str(night));save(pd.DataFrame(checks),'00_plan_audit/artifact_checks.tsv');print('ARTIFACT_CHECKS_PASSED',len(checks),'figures',len(fig),'review_MB',pack.stat().st_size/1024**2)
