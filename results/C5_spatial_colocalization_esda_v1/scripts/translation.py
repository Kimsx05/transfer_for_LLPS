"""Approximate hex-lattice translation sensitivity; all sections, no visual selection.
Changes at irregular tissue boundaries prevent an exact randomization test.
Statistics ALWAYS use esda.Moran_BV, permutations=0, on each overlap's induced W.
"""
from analyze import *
from scipy.signal import fftconvolve

def run(sid):
 dest=O/'05_distance_sensitivity'/f'{sid}_translation.tsv'
 if dest.exists():return
 d=pd.read_parquet(O/'03_local_bivariate_moran'/f'{sid}.parquet');d=d[d.analysis_included].reset_index(drop=True);rc=d[['array_row','array_col']].to_numpy(int);rc-=rc.min(0);lookup={tuple(v):i for i,v in enumerate(rc)};nr,nc=rc.max(0)+1;mask=np.zeros((nr,nc));mask[rc[:,0],rc[:,1]]=1;over=np.rint(fftconvolve(mask,mask[::-1,::-1],mode='full'));rr,cc=np.where(over>=.5*len(d));offsets=[(r-(nr-1),c-(nc-1)) for r,c in zip(rr,cc) if (r-(nr-1),c-(nc-1))!=(0,0) and ((r-(nr-1))-(c-(nc-1)))%2==0];rng=np.random.default_rng(123);rng.shuffle(offsets);offsets=offsets[:99];original=graph(d);rows=[]
 for dr,dc in offsets:
  pair=[(i,lookup[(r+dr,c+dc)]) for i,(r,c) in enumerate(rc) if (r+dr,c+dc) in lookup];baseids=np.array([a for a,b in pair]);shiftids=np.array([b for a,b in pair]);active=set(baseids);keep=np.array([any(j in active for j in original[i]) for i in baseids]);baseids=baseids[keep];shiftids=shiftids[keep];remap={v:i for i,v in enumerate(baseids)};nbs={remap[i]:[remap[j] for j in original[i] if j in remap] for i in baseids};w=weights(nbs)
  for p in PARTNERS:
   for method in ['raw','adjusted']:
    x=d.five_UCell.to_numpy() if method=='raw' else d[p+'__adjusted_five'].to_numpy();y=d['major_mean__'+p].to_numpy() if method=='raw' else d[p+'__adjusted_abundance'].to_numpy()
    if np.std(x[baseids])<1e-12 or np.std(y[shiftids])<1e-12:continue
    observed=esda.Moran_BV(x[baseids],y[baseids],w,permutations=0).I;shift=esda.Moran_BV(x[baseids],y[shiftids],w,permutations=0).I;rows.append(dict(section_id=sid,partner=p,method=method,dr=dr,dc=dc,n_overlap=len(baseids),observed_overlap_I=observed,shifted_I=shift,difference=shift-observed))
 save(rows,dest);print('TRANSLATE',sid,len(offsets),flush=True)
if __name__=='__main__':
 for s in sorted(META.index):run(s)
 d=pd.concat([pd.read_csv(f,sep='\t') for f in (O/'05_distance_sensitivity').glob('*_translation.tsv')]);summ=[]
 for (sid,p,m),z in d.groupby(['section_id','partner','method']):
  a=z.difference.to_numpy();pval=min(1,2*min((np.sum(a>=0)+1)/(len(a)+1),(np.sum(a<=0)+1)/(len(a)+1)));summ.append(dict(section_id=sid,partner=p,method=m,n_translations=len(a),overlap_min=z.n_overlap.min(),median_observed_overlap_I=z.observed_overlap_I.median(),median_shifted_I=z.shifted_I.median(),p_approx=pval))
 z=pd.DataFrame(summ);z['FDR_approx']=z.groupby('method').p_approx.transform(lambda x:bh(x));save(z,O/'tables/translation_sensitivity.tsv')
