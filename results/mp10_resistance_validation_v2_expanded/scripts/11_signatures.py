from env import *
import importlib
mod=importlib.import_module('06_prepare');rows=[];manifest=[]
def add(name,genes,weights=None,method='ssGSEA',direction=1,source='',property='published_members_adapted_score',training='',expected=None):
 weights=weights or [1]*len(genes)
 for i,(g,w) in enumerate(zip(genes,weights)):
  sym,how=mod.mapgene(g,g);rows.append(dict(signature=name,original_member=g,gene='' if sym.startswith('UNMAPPED|') else sym,weight=w,position=i+1,mapping=how))
 manifest.append(dict(signature=name,n_original=len(genes),method=method,expected_resistance_direction=direction,source=source,property=property,training=training,status='READY' if genes else 'BLOCKED_DATA',reason='' if genes else str(expected)))
cfg=json.loads((V/'config/frozen_v1.json').read_text())
for k,gs in cfg['genes'].items():
 name='MP10' if k=='MP10' else ('TSUNODA_'+k)
 add(name,gs,direction=(-1 if k.startswith('DN') else 1),source='v1 frozen Top20' if k=='MP10' else 'MSigDB official Tsunoda snapshot + original Tables1/2',training='BLCA discovery Epi' if k=='MP10' else 'T24/KK47 cisplatin models')
add('CisSig19','ADAT2 ATP1B3 CDIN1 C1QBP CDC7 CDCA7 FKBP14 KRT5 LRRC8C LY6K MMP10 NPM3 PSAT1 RIOK1 SLFN11 STOML2 USP31 WDR3 ZNF750'.split(),method='median_z',direction=-1,source='s41698-023-00375-y; cissig c8317dc0a6ab6b8ee66dc98f2e817014f1a3c2a5',property='native',training='GDSC2 pan-cancer cisplatin')
add('CisSig_gene_set_adapted','ADAT2 ATP1B3 CDIN1 C1QBP CDC7 CDCA7 FKBP14 KRT5 LRRC8C LY6K MMP10 NPM3 PSAT1 RIOK1 SLFN11 STOML2 USP31 WDR3 ZNF750'.split(),direction=-1,source='CisSig19 published members; Epi enrichment only',property='single_cell_adapted_not_native',training='GDSC2 pan-cancer cisplatin')
add('Xylinas9','EPAS1 ADD1 TLR4 ZNF582 GCNT4 HOXA9 RAPGEF5 DBNDD2 TSTD1'.split(),weights=[1]*5+[-1]*4,method='upmean_minus_dnmean_z',source='PMC5039423',training='BLCA intrinsic sensitivity cell-line comparison')
add('NPEPPS_VRAC_ratio','NPEPPS LRRC8A LRRC8B LRRC8C LRRC8D LRRC8E'.split(),method='linear_ratio',source='PMC11641020',property='native_linear_ratio',training='GSE171537-related resistant model system')
for p in sorted(S.glob('*CISPLATIN_RESISTANCE_*.json'))+sorted(S.glob('HALLMARK_*.json'))+list(S.glob('BENPORATH_*.json')):
 j=json.loads(p.read_text());name=next(iter(j));v=j[name];gs=v.get('geneSymbols',[])
 if not gs:continue
 add(name,gs,direction=0 if name.startswith(('HALLMARK','BENPORATH')) else (-1 if name.endswith('_DN') else 1),source='MSigDB official snapshot '+str(v.get('pmid',''))+' '+str(v.get('exactSource','')),property='state_control' if name.startswith(('HALLMARK','BENPORATH')) else 'published_members_adapted_score',training=str(v.get('description',v.get('descriptionBrief',''))))
t=(S/'Tanaka_supp_actual.txt').read_text();up=t.split('Up-regulated (25)')[1].split('Down-regulated (24)')[0];dn=t.split('Down-regulated (24)')[1].split('Supplementary Table S2')[0]
clean=lambda t:[g.strip().rstrip('.') for g in re.sub(r'\s+',' ',t).strip().split(',')]
up,dn=clean(up),clean(dn);assert len(up)==25 and len(dn)==24
add('Tanaka49',up+dn,weights=[1]*25+[-1]*24,method='upmean_minus_dnmean_z',source='PMC5085609 Supplementary Table S1',training='GSE58624 T24PR/5637PR; internal source comparison only')
c=read('01_signatures/CrM_table_1.tsv');c=c[c['0'].str.fullmatch(r'\d+',na=False)].copy();assert len(c)==63
crgenes=c['1'].str.replace(r'a\s*\[\d+\]','',regex=True).str.strip().to_list();cw=pd.to_numeric(c['3']).to_list()
add('CrM63_weighted_z',crgenes,cw,method='weighted_z',source='DOI10.1016/j.canlet.2024.217339 Supplementary Table2; z-score x published motility weight',property='adapted_weighted_mean_preserving_published_weights',training='GSE210921/GSE210953; GSE190636 comparison not independent training; gemcitabine/motility')
add('DTP124',[],method='signed_mean_z',source='PMC8437243 TableS6 NIHMS1732990-supplement-14.xlsx',property='native_BLOCKED',training='mouse embryonic pausing and CRC irinotecan DTP',expected='TableS6 download returns HTML/403/404; full signed member list and cross-species map not yet recovered')
add('CLSPN_published_program',[],source='PMC10992486 Supplemental TableS1 / PDF2',property='published_program_BLOCKED',training='GSE217413',expected='PDF2 contains full13743-gene CPM table, no identified2284-UP membership; cannot call all rows UP')
# Derived program: fixed before inspecting associations, no DESeq2 on CPM.
spec=B/'00_plan_audit/CLSPN_derivation_spec.json';sp=json.loads(spec.read_text());sp['timestamp_Asia_Shanghai']=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat();spec.write_text(json.dumps(sp,indent=2))
a=read('00_plan_audit/prepared/GSE217413_expression.tsv.gz').set_index('gene');cr=a.filter(like='CisR').mean(axis=1);wt=a.filter(like='WT').mean(axis=1);fc=np.log2((cr+.1)/(wt+.1));sel=pd.DataFrame({'gene':a.index,'log2ratio':fc.values,'meanCPM':a.mean(axis=1).values});sel=sel[(sel.meanCPM>=1)&(sel.log2ratio.abs()>=2)&~sel.gene.str.startswith('UNMAPPED|')];up=sel[sel.log2ratio>0].sort_values(['log2ratio','gene'],ascending=[False,True]).head(100);dn=sel[sel.log2ratio<0].sort_values(['log2ratio','gene']).head(100);save(pd.concat([up,dn]),'01_signatures/CLSPN_derived_selection.tsv');add('CLSPN_derived_program',up.gene.to_list()+dn.gene.to_list(),[1]*len(up)+[-1]*len(dn),method='upmean_minus_dnmean_z',source='CLSPN_derivation_spec.json',property='derived',training='GSE217413, internal only')
add('Choi_p53_like_classifier',[],method='classifier',source='Choi2014 original classifier required',property='classifier_BLOCKED',expected='Full original fitted classifier/centroids not recovered; existing author subtype labels retained where available; Hallmark P53 is distinct')
# Prespecified overlap-removal sensitivity for enrichment sets, never substitute originals.
mp=set(cfg['genes']['MP10']);initial=pd.DataFrame(rows)
for m in list(manifest):
 if m['signature']=='MP10' or m['method']!='ssGSEA':continue
 rr=initial[initial.signature==m['signature']];inter=set(rr.gene)&mp
 if inter:add(m['signature']+'_noMP10_overlap',rr.loc[~rr.gene.isin(mp),'original_member'].to_list(),direction=m['expected_resistance_direction'],source=m['source'],property='sensitivity_remove_MP10_overlap',training=m['training'])
genes=pd.DataFrame(rows);save(genes,'01_signatures/members.tsv');mf=pd.DataFrame(manifest);mf['MP10_intersection']=mf.signature.map(lambda x:';'.join(sorted(set(genes[genes.signature==x].gene)&mp)));save(mf,'01_signatures/signature_manifest.tsv')
with (B/'01_signatures/sets.gmt').open('w') as f:
 for m in manifest:
  rr=genes[(genes.signature==m['signature'])&(genes.gene!='')];f.write('\t'.join([m['signature'],m['source']]+rr.gene.to_list())+'\n')
print('SIGNATURES',len(mf),'members',len(genes),mf[['signature','n_original','status']].to_string(index=False))
