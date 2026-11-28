from pathlib import Path
import os,json,itertools
import numpy as np,pandas as pd,scipy.sparse as sp
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import font_manager
from PIL import Image
from scipy import stats
O=Path(__file__).resolve().parents[1];F=O/'09_figures';R=Path('/home/jinshengxi/LLPS/spatial_data/LLPS_spatial_transfer_manual');d=pd.read_parquet(O/'04_niche_characterization/spot_metadata.parquet');meta=pd.read_csv(O/'00_input_audit/section_metadata.tsv',sep='\t').set_index('section_id');niches=sorted(d.niche.unique(),key=int)
# Arial is available via system fontconfig; register exact files if matplotlib does not discover them.
import subprocess
font=subprocess.check_output(['fc-match','-f','%{file}','Arial']).decode();font_manager.fontManager.addfont(font)
plt.rcParams.update({'font.family':'Arial','font.size':12,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':1,'lines.linewidth':1.1,'pdf.fonttype':42,'figure.facecolor':'white','savefig.facecolor':'white'})
all_niches=[str(i) for i in range(1+max(int(d[c].astype(int).max()) for c in ['recommended_niche','major_best_niche','subtype_best_niche','kmeans_niche']))];colors={n:matplotlib.colors.to_hex(plt.get_cmap('tab20')(int(n)%20)) for n in all_niches};pd.DataFrame({'niche':list(colors),'color':list(colors.values())}).to_csv(F/'niche_colors.tsv',sep='\t',index=False)
def export(fig,path,legend=None,colorbar=None):
 path=F/path;path.parent.mkdir(exist_ok=True,parents=True);fig.savefig(str(path)+'.pdf',bbox_inches='tight');fig.savefig(str(path)+'.png',dpi=300,bbox_inches='tight')
 if legend is not None:
  fig.canvas.draw();bbox=legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted());fig.savefig(str(path)+'_legend_only.pdf',bbox_inches=bbox.expanded(1.05,1.05));legend.remove();fig.savefig(str(path)+'_panel_only.pdf',bbox_inches='tight')
 if colorbar is not None:
  fig.canvas.draw();bbox=colorbar.ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted());fig.savefig(str(path)+'_legend_only.pdf',bbox_inches=bbox.expanded(1.05,1.05));colorbar.remove();fig.savefig(str(path)+'_panel_only.pdf',bbox_inches='tight')
 plt.close(fig)
def heat(tab,path,label='Value',cmap='viridis',center=False):
 fig,ax=plt.subplots(figsize=(max(5,len(tab.columns)*.5),max(4,len(tab)*.32)));v=tab.to_numpy(float);limit=np.nanmax(np.abs(v)) if center else None;im=ax.imshow(v,aspect='auto',cmap=cmap,vmin=-limit if center else None,vmax=limit);ax.set_xticks(range(len(tab.columns)),tab.columns,rotation=60,ha='right');ax.set_yticks(range(len(tab.index)),tab.index);cb=fig.colorbar(im,ax=ax,label=label);export(fig,path,colorbar=cb)
comp=pd.read_csv(O/'04_niche_characterization/niche_composition.tsv',sep='\t');naming=[]
for n,g in comp[(comp.kind=='major')&(comp.scale=='center')&(comp.balance=='section')].groupby('niche'):
 g=g.sort_values('fraction',ascending=False);name=f'{g.iloc[0].celltype}-rich niche' if g.iloc[0].fraction>=.35 else f'Niche {n}';naming.append(dict(niche=str(n),name=name,top_celltype=g.iloc[0].celltype,top_fraction=g.iloc[0].fraction))
pd.DataFrame(naming).to_csv(O/'04_niche_characterization/niche_names.tsv',sep='\t',index=False)
for kind,scale,balance in itertools.product(['major','subtype'],['center','neighborhood'],['pooled','section']):
 t=comp[(comp.kind==kind)&(comp.scale==scale)&(comp.balance==balance)].pivot(index='niche',columns='celltype',values='fraction');heat(t,f'composition/{kind}_{scale}_{balance}_heatmap','Cell fraction')
# Representative sections: best entropy per source, independent of scores.
rank=[]
for sid,g in d.groupby('section_id'):
 p=g.niche.value_counts(normalize=True).to_numpy();rank.append(dict(section_id=sid,dataset=g.dataset.iloc[0],entropy=-np.sum(p*np.log(p)),n_spots=len(g),coverage=g.R2_coverage.mean()))
rank=pd.DataFrame(rank);representatives=rank.sort_values(['entropy','coverage'],ascending=False).groupby('dataset',sort=False).head(1);representatives.to_csv(F/'representative_sections.tsv',sep='\t',index=False)
scorefields=['five_UCell','LLPS77_UCell','Fibroblast_abundance','Endothelial_abundance'];limits={f:(float(d[f].min()),float(d[f].max())) for f in scorefields};(F/'continuous_color_limits.json').write_text(json.dumps(limits,indent=2))
models=['recommended_niche','major_best_niche','subtype_best_niche'];unique=[];aliases=[]
for model in models:
 same=next((old for old in unique if d[old].equals(d[model])),None)
 aliases.append(dict(model=model,spatial_folder=same if same else model,identical_to=same))
 if same:continue
 unique.append(model)
pd.DataFrame(aliases).to_csv(F/'spatial_map_aliases.tsv',sep='\t',index=False)
for sid,g in d.groupby('section_id',sort=False):
 imagepath=R/meta.loc[sid,'histology_relative_path'];img=np.array(Image.open(imagepath));xy=g[['x','y']].to_numpy();(F/'spatial/HE').mkdir(parents=True,exist_ok=True);link=F/'spatial/HE'/f'{sid}.png'
 if not link.exists():link.symlink_to(imagepath)
 for model in unique:
  labs=sorted(g[model].unique(),key=int);palette={v:matplotlib.colors.to_hex(plt.get_cmap('tab20')(int(v)%20)) for v in labs}
  for overlay in [False,True]:
   fig,ax=plt.subplots(figsize=(6.5,5.4))
   if overlay:ax.imshow(img,origin='upper')
   ax.scatter(xy[:,0],xy[:,1],s=4,c=[palette[v] for v in g[model]],alpha=.65 if overlay else 1,linewidths=0)
   ax.set_xlim(0,img.shape[1]);ax.set_ylim(img.shape[0],0);ax.set_aspect('equal');ax.axis('off');leg=ax.legend(handles=[Line2D([],[],marker='o',ls='',color=palette[v],label=f'Niche {v}') for v in labs],loc='center left',bbox_to_anchor=(1,0.5),frameon=False)
   export(fig,f'spatial/{model}/{sid}_{"overlay" if overlay else "niche"}',legend=leg)
 # Full section R1 graph QC, every edge verified in analysis.
 A=sp.load_npz(O/'01_spatial_graph'/f'{sid}_R1.npz');u,v=sp.triu(A,k=1).nonzero();fig,ax=plt.subplots(figsize=(6,5));from matplotlib.collections import LineCollection
 ax.add_collection(LineCollection(np.stack([xy[u],xy[v]],axis=1),linewidths=1,colors='#999999',alpha=.3));ax.scatter(xy[:,0],xy[:,1],s=2,c='black');ax.autoscale();ax.invert_yaxis();ax.set_aspect('equal');ax.axis('off');export(fig,f'graph_QC/{sid}')
 if sid in set(representatives.section_id):
  for f in scorefields:
   for normalized in [False,True]:
    field=f+'__z' if normalized else f
    fig,ax=plt.subplots(figsize=(6.8,5.4));lo,hi=(-3,3) if normalized else limits[f]
    # Standardized maps use full global range too, avoiding silent clipping.
    if normalized:lo,hi=float(d[field].min()),float(d[field].max())
    im=ax.scatter(xy[:,0],xy[:,1],s=6,c=g[field],cmap='magma',vmin=lo,vmax=hi,linewidths=0);ax.set_xlim(0,img.shape[1]);ax.set_ylim(img.shape[0],0);ax.set_aspect('equal');ax.axis('off');label={'five_UCell':'C5-associated score','LLPS77_UCell':'LLPS score','Fibroblast_abundance':'Fibroblast abundance','Endothelial_abundance':'Endothelial abundance'}[f];cb=fig.colorbar(im,ax=ax,label=label+(' (section z)' if normalized else ''));export(fig,f'representative/{sid}_{field}',colorbar=cb)
 print('figures',sid,flush=True)
# Neighborhood radius diagram from real representative section, no synthetic missing-tissue interpolation.
sid=representatives.iloc[0].section_id;g=d[d.section_id==sid];A=sp.load_npz(O/'01_spatial_graph'/f'{sid}_R1.npz');xy=g[['x','y']].to_numpy();center=int(np.argmax(g.R3_coverage));mask=np.zeros(len(g),bool);mask[center]=True
for r in [1,2,3]:
 mask=mask|(np.asarray(A@mask)>0);fig,ax=plt.subplots(figsize=(5,5));ax.scatter(xy[:,0],xy[:,1],s=10,c='#dddddd');ax.scatter(xy[mask,0],xy[mask,1],s=40,c='#637BA8');ax.scatter(*xy[center],s=70,c='#B2182B');span=60;ax.set_xlim(xy[center,0]-span,xy[center,0]+span);ax.set_ylim(xy[center,1]+span,xy[center,1]-span);ax.set_aspect('equal');ax.axis('off');export(fig,f'graph_QC/R{r}_example')
# Enrichment heatmaps, all sections plus source-level and overall summaries.
p=pd.read_csv(O/'07_neighborhood_enrichment/section_pair_results.tsv',sep='\t')
for sid,g in p.groupby('section_id'):
 heat(g.pivot(index='niche_a',columns='niche_b',values='Z'),f'adjacency/sections/{sid}_Z','Enrichment Z','RdBu_r',True)
for balance in ['section','patient','dataset']:
 t=pd.read_csv(O/'07_neighborhood_enrichment'/f'{balance}_balanced_summary.tsv',sep='\t');groups=t.groupby('dataset') if balance=='dataset' else [(balance,t)]
 for name,g in groups:
  for val in ['mean_OE','mean_Z']:heat(g.pivot(index='niche_a',columns='niche_b',values=val),f'adjacency/{name}_{val}',val,'RdBu_r' if val=='mean_Z' else 'viridis',val=='mean_Z')
q=pd.read_csv(O/'03_clustering/clustering_quality_metrics.tsv',sep='\t')
for y in ['n_niches','silhouette']:
 fig,ax=plt.subplots(figsize=(7,5))
 for name,g in q.groupby(['kind','radius','transform','knn']):ax.plot(g.resolution,g[y],'.-',label='_'.join(map(str,name)))
 ax.set_xlabel('Leiden resolution');ax.set_ylabel('Number of niches' if y=='n_niches' else 'Silhouette score');leg=ax.legend(fontsize=7,bbox_to_anchor=(1,1),loc='upper left');export(fig,f'QC/{y}',legend=leg)
a=pd.read_csv(O/'03_clustering/all_models_ARI_NMI.tsv',sep='\t');names=sorted(set(a.model_a)|set(a.model_b))
for val in ['ARI','NMI']:
 mat=pd.DataFrame(np.eye(len(names)),index=names,columns=names)
 for r in a.itertuples():mat.loc[r.model_a,r.model_b]=mat.loc[r.model_b,r.model_a]=getattr(r,val)
 heat(mat,f'QC/{val}',val)
for key in ['dataset','section_id']:
 t=pd.crosstab(d.niche,d[key],normalize='columns');heat(t,f'QC/niche_by_{key}','Spot fraction')
# Plot-ready tables for skill-based R templates (no analysis in R).
d[['niche']+scorefields+['A26_UCell','five_UCell__adjusted','LLPS77_UCell__adjusted']].to_csv(F/'score_plot_data.tsv.gz',sep='\t',index=False)
for role in ['major_best','subtype_best']:
 xy=np.load(O/'03_clustering'/f'{role}_umap.npy');z=d[['spot_id','dataset','section_id','kmeans_niche',role+'_niche']].copy();z['x']=xy[:,0];z['y']=xy[:,1];z.to_csv(F/f'{role}_umap_data.tsv.gz',sep='\t',index=False)
# Section-level scatter/correlations, section is descriptive unit; patient inferential summaries in tables.
s=d.groupby(['section_id','dataset'])[scorefields].mean().reset_index()
for f,h in [('five_UCell','Endothelial_abundance'),('LLPS77_UCell','Endothelial_abundance'),('five_UCell','LLPS77_UCell')]:
 fig,ax=plt.subplots(figsize=(5.5,5));ax.scatter(s[f],s[h],s=20,c='#637BA8');fit=np.polyfit(s[f],s[h],1);xx=np.linspace(s[f].min(),s[f].max(),100);ax.plot(xx,np.polyval(fit,xx),color='black');r=stats.pearsonr(s[f],s[h]).statistic;ax.text(.04,.96,f'Pearson R = {r:.2f}\nUnit: section (n = {len(s)})',transform=ax.transAxes,va='top');ax.set_xlabel({'five_UCell':'C5-associated score','LLPS77_UCell':'LLPS score','Endothelial_abundance':'Endothelial abundance'}[f]);ax.set_ylabel({'five_UCell':'C5-associated score','LLPS77_UCell':'LLPS score','Endothelial_abundance':'Endothelial abundance'}[h]);export(fig,f'correlations/{f}_vs_{h}')
c=pd.read_csv(O/'06_C5_endothelial/section_correlations.tsv',sep='\t');c=c[(c.niche=='ALL')&(c.scale=='center')&(c.method=='pearson')&(c.adjustment=='raw')]
for (f,h),g in c.groupby(['feature_a','feature_b']):
 fig,ax=plt.subplots(figsize=(7,15));ax.scatter(g.r,np.arange(len(g)),s=15);ax.axvline(0,color='grey');ax.set_yticks(np.arange(len(g)),g.section_id,fontsize=7);ax.set_xlabel('Within-section Pearson R (descriptive)');ax.set_xlim(-1,1);export(fig,f'correlations/{f}_vs_{h}_sections')
(F/'FIGURE_NOTES.md').write_text('Spatial maps, heatmaps, QC and correlation plots have no matching routine skill template; custom matplotlib export follows white background, Arial12, vector scatter/PDF,300dpi PNG and separate legend/panel PDFs. H&E is referenced by symlink rather than copied. All63 sections retained. Identical recommended/major/subtype labels share one spatial map set to avoid duplicate figures. Representative selection: maximum niche entropy per dataset, then coverage; no signature/clinical criterion. Continuous scales use global min/max, including section-standardized maps.\n')
# Auxiliary QC: coverage, component structure and scale sensitivity.
coverage=pd.read_csv(O/'01_spatial_graph/neighborhood_coverage_summary.tsv',sep='\t');fig,ax=plt.subplots(figsize=(5,4));ax.boxplot([coverage[coverage.radius==r].mean_coverage for r in [1,2,3]],tick_labels=['R1','R2','R3']);ax.set_ylabel('Mean neighborhood coverage');ax.set_ylim(0,1);export(fig,'QC/radius_coverage')
fragment=pd.read_csv(O/'03_clustering/spatial_fragmentation_metrics.tsv',sep='\t');chosen=['recommended_niche','R2_representative_niche','R3_representative_niche','subtype_best_niche','kmeans_niche'];fig,ax=plt.subplots(figsize=(8,5));ax.boxplot([fragment[fragment.model==m].largest_fraction for m in chosen],tick_labels=chosen);ax.tick_params(axis='x',rotation=60);ax.set_ylabel('Largest component / niche spots');ax.set_ylim(0,1);export(fig,'QC/fragmentation_comparison')
features=pd.read_csv(O/'06_C5_endothelial/niche_feature_summary.tsv',sep='\t',index_col=0);standard=(features-features.mean())/features.std();heat(standard,'correlations/niche_feature_heatmap','Feature z across niches','RdBu_r',True)
# Verified clinical groups only; descriptive section points, not iid spot statistics.
pr=pd.read_csv(O/'08_clinical_association/section_niche_proportions_complete.tsv',sep='\t')
for field,groups in [('tissue',['Normal','Tumor']),('source_sample_type',['NMIBC','MIBC'])]:
 z=pr[pr[field].isin(groups)]
 if not len(z):continue
 for niche,g in z.groupby('niche'):
  available=[v for v in groups if (g[field]==v).any()]
  fig,ax=plt.subplots(figsize=(4.8,4.9));vals=[g.loc[g[field]==v,'proportion'].to_numpy() for v in available];ax.boxplot(vals,tick_labels=available)
  for j,v in enumerate(vals):ax.scatter(np.full(len(v),j+1),v,c='#637BA8',s=20)
  ax.set_ylabel(f'Niche {niche} proportion');ax.set_ylim(0,1);export(fig,f'clinical/{field}_niche{niche}')
# Report adjusted contrast uncertainty directly at the verified-patient unit.
adj=pd.read_csv(O/'05_C5_LLPS_enrichment/fibroblast_adjusted_summary.tsv',sep='\t')
for feature,g in adj[adj.unit=='patient'].groupby('feature'):
 fig,ax=plt.subplots(figsize=(6,5));ax.errorbar(g.effect,np.arange(len(g)),xerr=[g.effect-g.ci_low,g.ci_high-g.effect],fmt='o',color='#637BA8');ax.set_yticks(np.arange(len(g)),['Niche '+str(x) for x in g.niche]);ax.axvline(0,color='grey');ax.set_xlabel('Adjusted standardized contrast (patient mean, 95% CI)');export(fig,f'correlations/{feature}_adjusted_patient_contrasts')
patient=d[d.patient_uid.notna()].groupby(['section_id','patient_uid'])[scorefields].mean().groupby('patient_uid').mean().reset_index()
for f,h in [('five_UCell','Endothelial_abundance'),('LLPS77_UCell','Endothelial_abundance'),('five_UCell','LLPS77_UCell')]:
 fig,ax=plt.subplots(figsize=(5.5,5));ax.scatter(patient[f],patient[h],s=20,c='#637BA8');fit=np.polyfit(patient[f],patient[h],1);xx=np.linspace(patient[f].min(),patient[f].max(),100);ax.plot(xx,np.polyval(fit,xx),color='black');res=stats.pearsonr(patient[f],patient[h]);ax.text(.04,.96,f'Pearson R = {res.statistic:.2f}\nP = {res.pvalue:.3g}\nUnit: patient (n = {len(patient)})',transform=ax.transAxes,va='top');ax.set_xlabel({'five_UCell':'C5-associated score','LLPS77_UCell':'LLPS score','Endothelial_abundance':'Endothelial abundance'}[f]);ax.set_ylabel({'five_UCell':'C5-associated score','LLPS77_UCell':'LLPS score','Endothelial_abundance':'Endothelial abundance'}[h]);export(fig,f'correlations/{f}_vs_{h}_patient')
