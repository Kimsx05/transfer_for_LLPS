from env import *
from scipy.sparse import csc_matrix
import subprocess,time,concurrent.futures,shutil
os.environ.update(R_PROFILE='/dev/null',R_PROFILE_USER='/dev/null',OMP_THREAD_LIMIT='1',RCPP_PARALLEL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',BLIS_NUM_THREADS='1')
mode='official' if '--official' in sys.argv else 'cycle';C=P/'neighborhood_mp10_calibration_v5/01_score_fidelity/counts_cache';genes=(C/'genes.txt').read_text().splitlines();cells=(C/'cells.txt').read_text().splitlines();f=pd.read_csv(B/'02_epi/cell_average_tie_scores.tsv.gz',sep='\t',usecols=['cell_id','sample']);assert f.cell_id.tolist()==cells
ind=np.memmap(C/'indices.bin',dtype='<i4',mode='r');dat=np.memmap(C/'data.bin',dtype='<f8',mode='r');ptr=np.memmap(C/'indptr.bin',dtype='<i4',mode='r');x=csc_matrix((dat,ind,ptr),shape=(len(genes),len(cells)),copy=False);root=B/'02_epi/specimen_checkpoints';root.mkdir(exist_ok=True);log=B/'logs';jobs=list(f.groupby('sample',sort=True).indices.items())
def run(item):
 sample,idx=item;label=sample.replace('|','__');d=root/label;d.mkdir(exist_ok=True);target=d/('official_scores.tsv.gz' if mode=='official' else 'cycle_scores.tsv')
 if target.exists():return dict(sample=sample,status='COMPLETE_REUSED',n=len(idx))
 t=time.time();xs=x[:,idx];xs.data.tofile(d/'data.bin');xs.indices.astype('<i4').tofile(d/'indices.bin');xs.indptr.astype('<i4').tofile(d/'indptr.bin');(d/'genes.txt').write_text('\n'.join(genes)+'\n');(d/'cells.txt').write_text('\n'.join(np.array(cells)[idx])+'\n');del xs
 with (log/(mode+'_'+label+'.log')).open('w') as out:
  cmd=['/usr/bin/time','-v','Rscript',str(B/'scripts/17_epi_official.R'),str(d)]+(['--cycle-only'] if mode=='cycle' else []);r=subprocess.run(cmd,stdout=out,stderr=subprocess.STDOUT)
 # Only this task's transient exact sample slices are removed after successful processing.
 if r.returncode==0:
  for name in ['data.bin','indices.bin','indptr.bin']:(d/name).unlink()
 result=dict(sample=sample,n=len(idx),seconds=time.time()-t,returncode=r.returncode,status='COMPLETE' if r.returncode==0 else 'FAILED');print(result,flush=True);return result
workers=int(os.environ.get('MP10_EPI_WORKERS','2'));assert 1<=workers<=3, 'reserve one CPU for orchestrator within four-thread budget'
with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:results=list(ex.map(run,jobs))
save(pd.DataFrame(results),f'02_epi/{mode}_execution.tsv')
files=list(root.glob('*/'+('official_scores.tsv.gz' if mode=='official' else 'cycle_scores.tsv')));all_=pd.concat([pd.read_csv(p,sep='\t') for p in files],ignore_index=True);save(all_,'02_epi/'+('cell_official_scores.tsv.gz' if mode=='official' else 'cell_cycle_scores.tsv.gz'))
assert all(r['status'].startswith('COMPLETE') for r in results),'some specimens failed; completed checkpoints preserved'
print('ALL_COMPLETE',mode,len(all_),flush=True)
