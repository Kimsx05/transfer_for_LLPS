#!/usr/bin/env python3
import os
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
os.environ['MPLBACKEND']='Agg'
import argparse,gc,json,re,resource,subprocess,time,traceback
from datetime import datetime,timezone
from pathlib import Path
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np,pandas as pd,scvi,torch
from scipy import sparse
from scipy.stats import spearmanr
from cell2location.models import Cell2location

ap=argparse.ArgumentParser();ap.add_argument('--section',required=True);args=ap.parse_args()
P=Path('/home/data/t070721/codex_workspace/Bladder Metabolism');R=P/'cell2location_full_run/v2_30type';O=R/'02_sections'/args.section;O.mkdir(parents=True,exist_ok=True)
torch.set_num_threads(1);torch.set_num_interop_threads(1)
SIG=P/'subtype/v2/reference_training_v2/11_final_spatial_reference/30type_regression_model/learned_signatures_30type.csv.gz'
MAN=R/'01_manifest/01_spatial_manifest_frozen.csv';HE=P/'cell2location_full_run/v1/04_figures/HE/HE_image_manifest.tsv'
STATUS=O/'status.json';seed=20260910
def now():return datetime.now(timezone.utc).isoformat()
def atomjson(p,x):t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(x,indent=2,default=str)+'\n');os.replace(t,p)
def say(x):print(f'[{now()}] {x}',flush=True)
def sync_roll():
    p=R/'04_incremental_tables/rolling_section_status.csv'
    if not p.exists():return
    x=pd.read_csv(p,dtype=object).fillna('');i=x.section.eq(args.section)
    x.loc[i,'status']=status.get('status','');x.loc[i,'start_time']=status.get('start_time','');x.loc[i,'finish_time']=status.get('finish_time','');x.loc[i,'runtime']=status.get('runtime_sec','');x.loc[i,'n_spots']=status.get('n_spots','');x.loc[i,'n_genes']=status.get('n_genes','');x.loc[i,'posterior_ok']=str(status.get('posterior_ok',False));x.loc[i,'plots_ok']=str(status.get('plots_ok',False));x.loc[i,'notes']=status.get('error','');tmp=p.with_suffix('.csv.tmp');x.to_csv(tmp,index=False);os.replace(tmp,p)
def update(**kw):status.update(kw);atomjson(STATUS,status);sync_roll()
def extract(a,key,prefix,labels):
    x=pd.DataFrame(a.obsm[key],index=a.obs_names);x.columns=[str(c).removeprefix(prefix) for c in x.columns]
    if set(x.columns)!=set(labels):raise RuntimeError(f'{key} label mismatch')
    return x[labels]
def safe(x):return re.sub(r'[^A-Za-z0-9._-]+','_',x)
def plot_one(xy,val,path,transparent=False,img=None):
    val=np.asarray(val,float);vmax=float(np.nanquantile(val,.992));vmax=vmax if vmax>0 else 1.0
    fig,ax=plt.subplots(figsize=(6,6));fig.patch.set_alpha(0 if transparent else 1);ax.set_facecolor('none' if transparent else 'white')
    if img is not None:ax.imshow(img,origin='upper')
    q=ax.scatter(xy[:,0],xy[:,1],c=val,s=13,cmap='magma',vmin=0,vmax=vmax,edgecolors='none',alpha=.88)
    ax.set_aspect('equal');ax.invert_yaxis() if img is None else None;ax.axis('off');fig.colorbar(q,ax=ax,fraction=.045,pad=.01);fig.tight_layout(pad=.1);fig.savefig(path,dpi=220,transparent=transparent);plt.close(fig);return vmax

row=pd.read_csv(MAN).query('section_id == @args.section')
if len(row)!=1:raise RuntimeError('section missing/nonunique in frozen manifest')
row=row.iloc[0];status={'section':args.section,'dataset':row.dataset,'status':'RUNNING','start_time':now(),'finish_time':'','runtime_sec':'','n_spots':'','n_genes':'','n_common_genes':'','posterior_ok':False,'plots_ok':False,'error':''};atomjson(STATUS,status);sync_roll();t0=time.time()
try:
    say('THREADS OMP=1 MKL=1 OPENBLAS=1 NUMEXPR=1 affinity<=8')
    say('GPU START '+subprocess.run(['nvidia-smi','--query-gpu=name,memory.total,memory.free,utilization.gpu','--format=csv,noheader'],capture_output=True,text=True).stdout.strip())
    sig=pd.read_csv(SIG,index_col=0);labels=list(sig.columns)
    if sig.shape!=(17987,30):raise RuntimeError(f'signature shape {sig.shape}')
    a=ad.read_h5ad(row.raw_input_path)
    if not sparse.issparse(a.X) or a.X.dtype.kind not in 'iu':raise RuntimeError('spatial X is not sparse integer raw counts')
    genes=[g for g in sig.index if g in set(a.var_names)]
    if len(genes)<12000:raise RuntimeError(f'exact common genes too low {len(genes)}')
    pd.DataFrame([{'section':args.section,'n_ref_genes':len(sig),'n_spatial_genes':a.n_vars,'n_common_genes':len(genes),'common_fraction_reference':len(genes)/len(sig),'common_fraction_spatial':len(genes)/a.n_vars}]).to_csv(O/'common_gene_audit.csv',index=False)
    v=a[:,genes].copy();del a
    if not sparse.isspmatrix_csr(v.X):v.X=sparse.csr_matrix(v.X)
    v.obs['section_batch']=pd.Categorical([args.section]*v.n_obs);Cell2location.setup_anndata(v,batch_key='section_batch')
    scvi.settings.seed=seed;np.random.seed(seed);torch.manual_seed(seed);torch.cuda.manual_seed_all(seed)
    cfg={'N_cells_per_location':20,'detection_alpha':20,'max_epochs':30000,'batch_size':None,'train_size':1.0,'lr':0.002,'seed':seed,'posterior_num_samples':1000,'posterior_batch_size_max':2500,'labels':labels}
    atomjson(O/'model_config.json',cfg);torch.cuda.reset_peak_memory_stats()
    reusable=(O/'model/model.pt').exists() and (O/'training_history.csv.gz').exists()
    if reusable:
        say(f'RESUME SAVED MODEL section={args.section}; training will not be repeated')
        m=Cell2location.load(O/'model',adata=v,accelerator='gpu',device=1)
        hd=pd.read_csv(O/'training_history.csv.gz')
        if hd.value.size==0 or not np.isfinite(hd.value).all():raise RuntimeError('saved training history is non-finite/empty')
        train_sec=0.0
    else:
        m=Cell2location(v,cell_state_df=sig.loc[genes],N_cells_per_location=20,detection_alpha=20)
        say(f'TRAIN START section={args.section} shape={v.shape} common={len(genes)}');tt=time.time()
        m.train(max_epochs=30000,batch_size=None,train_size=1,lr=.002,accelerator='gpu',device=1,enable_checkpointing=False,enable_progress_bar=False)
        train_sec=time.time()-tt;hist=[]
        for key,val in m.history_.items():
            for i,x in enumerate(np.asarray(val).reshape(-1)):hist.append({'metric':key,'epoch':i,'value':float(x)})
        hd=pd.DataFrame(hist);hd.to_csv(O/'training_history.csv.gz',index=False,compression='gzip')
        if hd.value.size==0 or not np.isfinite(hd.value).all():raise RuntimeError('non-finite/empty training history')
        m.save(O/'model',overwrite=False,save_anndata=False)
    update(status='POSTERIOR_EXPORTING',training_seconds=train_sec,reused_trained_model=reusable)
    batch=min(2500,v.n_obs);pt=time.time()
    try:v=m.export_posterior(v,sample_kwargs={'num_samples':1000,'batch_size':batch,'accelerator':'gpu','device':1})
    except torch.cuda.OutOfMemoryError:
        torch.cuda.empty_cache();gc.collect();batch=min(1024,v.n_obs);v=m.export_posterior(v,sample_kwargs={'num_samples':1000,'batch_size':batch,'accelerator':'gpu','device':1})
    post_sec=time.time()-pt
    mean=extract(v,'means_cell_abundance_w_sf','meanscell_abundance_w_sf_',labels);q05=extract(v,'q05_cell_abundance_w_sf','q05cell_abundance_w_sf_',labels);q95=extract(v,'q95_cell_abundance_w_sf','q95cell_abundance_w_sf_',labels)
    total=mean.sum(1);prop=mean.div(total.replace(0,np.nan),axis=0).fillna(0)
    for name,x in [('mean',mean),('q05',q05),('q95',q95),('proportion',prop)]:x.to_csv(O/f'abundance_{name}.csv.gz',compression='gzip')
    arrays=[mean.to_numpy(),q05.to_numpy(),q95.to_numpy(),prop.to_numpy()]
    posterior_ok=(v.n_obs==len(mean) and mean.shape[1]==30 and all(np.isfinite(x).all() for x in arrays) and all((x>=0).all() for x in arrays[:3]) and float(mean.to_numpy().sum())>0)
    if not posterior_ok:raise RuntimeError('posterior integrity failure')
    # HDF5 keys cannot contain '/'. Keep exact biological labels in exported
    # abundance CSVs and reversibly sanitize only DataFrame column keys in H5AD.
    column_map={}
    for key,val in list(v.obsm.items()):
        if isinstance(val,pd.DataFrame):
            old=[str(x) for x in val.columns];new=[x.replace('/','__SLASH__') for x in old]
            if new!=old:
                val=val.copy();val.columns=new;v.obsm[key]=val
                column_map[key]=dict(zip(new,old))
    v.uns['celltype_label_column_mapping_json']=json.dumps(column_map,sort_keys=True)
    v.uns['v2_30type_full_run_config']=cfg;tmp=O/'spatial_posterior.h5ad.incomplete'
    if tmp.exists():tmp.unlink()
    v.write_h5ad(tmp,compression='gzip');os.replace(tmp,O/'spatial_posterior.h5ad')
    # Immediate section plots.
    xy=np.asarray(v.obsm['spatial'],float);plotdir=R/'03_incremental_plots/by_section'/args.section;plotdir.mkdir(parents=True,exist_ok=True)
    dsdir=R/'03_incremental_plots/by_dataset'/safe(str(row.dataset));dsdir.mkdir(parents=True,exist_ok=True);link=dsdir/args.section
    if not link.exists():link.symlink_to(plotdir,target_is_directory=True)
    hm=pd.read_csv(HE,sep='\t');hp=hm.loc[hm.section_id.eq(args.section),'image_path'];img=plt.imread(hp.iloc[0]) if len(hp) and Path(hp.iloc[0]).exists() else None
    scales=[]
    for metric,tab in [('mean',mean),('q05',q05),('proportion',prop)]:
        for ct in labels:
            base=f'{args.section}__{safe(ct)}__{metric}';scales.append({'section':args.section,'cell_type':ct,'metric':metric,'p99_2':plot_one(xy,tab[ct],plotdir/f'{base}__white.png')})
            plot_one(xy,tab[ct],plotdir/f'{base}__transparent.png',transparent=True)
            if img is not None:plot_one(xy,tab[ct],plotdir/f'{base}__HE_overlay.png',img=img)
    pd.DataFrame(scales).to_csv(O/'plot_scale_audit.csv',index=False)
    families={'CD4':['Treg','CD4 Helper/Conventional','CD4 Naive/Memory-like','Tfh/CXCL13-like'],'CD8_NK':['CD8 GZMK/Effector-memory-like','CD8 Exhausted','NK'],'Endothelial':['Capillary EC','Venous EC','Arterial EC'],'TAM':['C1QC TAM','SPP1 TAM'],'DC':['cDC1','cDC2'],'Fibro':['iCAF','Matrix_CAF']}
    corr=[];watch=[];wdir=R/'03_incremental_plots/identifiability_watch'/args.section;wdir.mkdir(parents=True,exist_ok=True)
    for fam,cts in families.items():
        cts=[x for x in cts if x in labels]
        famrho={}
        for metric,tab in [('mean',mean),('q05',q05)]:
            vals=[]
            for x,y in __import__('itertools').combinations(cts,2):vals.append((x,y,float(spearmanr(tab[x],tab[y]).statistic)))
            for x,y,rho in vals:corr.append({'section':args.section,'subtype_A':x,'subtype_B':y,'spearman':rho,'family':fam,'metric':metric})
            famrho[metric]=max([x[2] for x in vals],default=np.nan)
        watch.append({'section':args.section,'dataset':row.dataset,'family':fam,'max_sibling_rho_mean':famrho['mean'],'max_sibling_rho_q05':famrho['q05']})
        for metric,tab in [('mean',mean),('q05',q05),('proportion',prop)]:plot_one(xy,tab[cts].sum(1),wdir/f'{args.section}__{fam}__family_total__{metric}__HE_overlay.png',img=img)
    pd.DataFrame(corr).to_csv(O/'section_sibling_correlation.csv',index=False);pd.DataFrame(watch).to_csv(O/'identifiability_watch_summary.csv',index=False)
    (O/'COMPLETED_AND_VALID').write_text(now()+'\n');update(status='COMPLETED',finish_time=now(),runtime_sec=time.time()-t0,n_spots=v.n_obs,n_genes=v.n_vars,n_common_genes=len(genes),posterior_ok=True,plots_ok=True,training_seconds=train_sec,posterior_seconds=post_sec,peak_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024/1024,peak_vram_mb=torch.cuda.max_memory_allocated()/1024**2)
    say('GPU END '+subprocess.run(['nvidia-smi','--query-gpu=name,memory.total,memory.free,utilization.gpu','--format=csv,noheader'],capture_output=True,text=True).stdout.strip());say('COMPLETED_AND_VALID')
except Exception as e:
    update(status='FAILED',finish_time=now(),runtime_sec=time.time()-t0,error=repr(e),traceback=traceback.format_exc());raise
