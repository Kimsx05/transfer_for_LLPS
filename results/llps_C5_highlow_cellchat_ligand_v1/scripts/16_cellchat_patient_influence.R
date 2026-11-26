source('/home/jinshengxi/bladder_metabolism/results/llps_C5_highlow_cellchat_ligand_v1/scripts/common.R')
p<-fread(file.path(OUT,'03_patient_level_comparison/within_median_t10_patient_pairs.csv.gz'))[analysis_condition=='Tumor'&kind=='network']
r<-p[,{vals<-difference;med<-median(vals);lopo<-vapply(seq_along(vals),function(i)median(vals[-i]),numeric(1));list(n_patients=.N,paired_median=med,LOPO_median_min=min(lopo),LOPO_median_max=max(lopo),LOPO_same_sign_fraction=mean(sign(lopo)==sign(med)),largest_absolute_patient=patient_id[which.max(abs(vals))],largest_absolute_difference=max(abs(vals)),nonzero_patients=sum(vals!=0))},by=.(direction,partner)]
writecsv(r,file.path(OUT,'07_sensitivity_analysis/CellChat_network_leave_patient_out.csv'))
cat('PATIENT INFLUENCE COMPLETE\n')
