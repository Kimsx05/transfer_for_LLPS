from pathlib import Path
import numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import font_manager
import subprocess
O=Path(__file__).resolve().parents[1];F=O/'09_figures/control_spatial_QC';F.mkdir(exist_ok=True);font_manager.fontManager.addfont(subprocess.check_output(['fc-match','-f','%{file}','Arial']).decode());plt.rcParams.update({'font.family':'Arial','font.size':12,'pdf.fonttype':42,'figure.facecolor':'white','savefig.facecolor':'white'})
d=pd.read_parquet(O/'04_niche_characterization/spot_metadata.parquet');meta=pd.read_csv(O/'00_input_audit/section_metadata.tsv',sep='\t').set_index('section_id');rep=pd.read_csv(O/'09_figures/representative_sections.tsv',sep='\t');models=['kmeans_niche','R2_representative_niche','R3_representative_niche','native_R1','native_R2','native_R3']
for sid in rep.section_id:
 g=d[d.section_id==sid]
 for model in models:
  labels=g[model].astype(str);levels=sorted(labels.unique(),key=int);native=model.startswith('native');levels=labels.value_counts().head(5).index.tolist() if native else levels;K=d[model].nunique();colors={v:matplotlib.colors.to_hex(plt.get_cmap('hsv')(int(v)/K)) if native else matplotlib.colors.to_hex(plt.get_cmap('tab20')(int(v)%20)) for v in levels};pointcolors=[colors.get(v,'#DDDDDD') for v in labels];fig,ax=plt.subplots(figsize=(7,5.4));ax.scatter(g.x,g.y,c=pointcolors,s=5,linewidths=0);ax.set_xlim(0,meta.loc[sid,'image_width']);ax.set_ylim(meta.loc[sid,'image_height'],0);ax.set_aspect('equal');ax.axis('off');handles=[Line2D([],[],ls='',marker='o',color=colors[v],label=f'Cluster {v}' if native else f'Niche {v}') for v in levels]
  if native:handles.append(Line2D([],[],ls='',marker='o',color='#DDDDDD',label='Other native clusters'))
  leg=ax.legend(handles=handles,frameon=False,bbox_to_anchor=(1,.5),loc='center left');path=F/(sid+'_'+model);fig.savefig(str(path)+'.pdf',bbox_inches='tight');fig.savefig(str(path)+'.png',dpi=300,bbox_inches='tight');fig.canvas.draw();bbox=leg.get_window_extent().transformed(fig.dpi_scale_trans.inverted());fig.savefig(str(path)+'_legend_only.pdf',bbox_inches=bbox.expanded(1.05,1.05));leg.remove();fig.savefig(str(path)+'_panel_only.pdf',bbox_inches='tight');plt.close(fig)
 print(sid,flush=True)
(F/'README.md').write_text('Score-blind representatives reused from primary maps. Coordinates and field of view match primary maps. Kmeans and R2/R3 continuous candidates display all niches. Native controls have hundreds of clusters: only the five largest within each section have individual colors; all remaining clusters are explicitly pooled in grey for QC readability, with complete original labels retained in spot metadata. Native colors are fixed by global cluster ID within each radius. Numeric IDs across methods are not biologically aligned.\n')
