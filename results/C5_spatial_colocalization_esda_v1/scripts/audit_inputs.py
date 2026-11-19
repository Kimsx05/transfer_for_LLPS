from analyze import *
if __name__=='__main__':
 rows=[];seen=set();q=[]
 for sid in META.index:
  d=pd.read_parquet(S/'01_spatial_scores'/f'{sid}.parquet');ids=set(d.spot_id);assert not ids.intersection(seen);seen|=ids
  assert len(ids)==len(d) and d.dataset.eq(META.loc[sid,'dataset']).all();assert d.five_covered_genes.eq(5).all();assert np.isfinite(d[['array_row','array_col','five_UCell','five_AUCell','nCount','nFeature']]).all().all()
  for p in PARTNERS+['Fibroblast','Mural']:
   y=d['major_mean__'+p];v=d['major_q05__'+p];assert (y>=0).all() and (v>=0).all();q.append(dict(section_id=sid,partner=p,mean_min=y.min(),mean_median=y.median(),mean_max=y.max(),q05_median=v.median(),q05_positive_fraction=np.mean(v>0),q05_to_mean_median=(v/y.replace(0,np.nan)).median(),q05_gt_mean_count=int((v>y).sum()),same_spot_Five_spearman=stats.spearmanr(d.five_UCell,y).statistic,score_constant=np.std(d.five_UCell)==0,target_constant=np.std(y)==0))
  rows.append(dict(section_id=sid,dataset=META.loc[sid,'dataset'],n_spots=len(d),unique_spot_ids=len(ids),covered_five=5,mean_columns=sum(c.startswith('major_mean__') for c in d),q05_columns=sum(c.startswith('major_q05__') for c in d),sub_mean_columns=sum(c.startswith('sub_mean__') for c in d),sub_q05_columns=sum(c.startswith('sub_q05__') for c in d),patient_verified=pd.notna(META.loc[sid,'patient_uid']),Five_UCell_zero_fraction=np.mean(d.five_UCell==0)))
 assert len(seen)==265002 and len(rows)==63;save(rows,O/'00_input_audit/input_QC.tsv');save(q,O/'00_input_audit/abundance_q05_audit.tsv')
 hashes=[]
 for f in [Path('results/C5_fivegene_LLPS77_spatial_v1/COMPREHENSIVE_REPORT.md'),S/'COMPREHENSIVE_REPORT.md',Path('results/llps_C5_cellchat_v2/README.md'),Path('results/llps_C5_cellchat_v2/05_threshold_comparison/C5_partner_network_all_thresholds.csv')]+[Path('results/llps_C5_cellchat_v2')/f/'C5_vs_fibro_common_specimen_strength.csv' for f in ['01_threshold10','02_threshold20','03_threshold30']]:hashes.append(dict(path=str(f.resolve()),sha256=hashlib.sha256(f.read_bytes()).hexdigest(),bytes=f.stat().st_size))
 save(hashes,O/'00_input_audit/context_source_hashes.tsv');print('63 sections / 265002 IDs / 5 frozen genes: PASS')
