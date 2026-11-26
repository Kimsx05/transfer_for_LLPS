source('/home/jinshengxi/bladder_metabolism/results/llps_C5_highlow_cellchat_ligand_v1/scripts/common.R')
f<-qs::qread(FIBRO,nthreads=1);m<-as.data.table(f[[]],keep.rownames='cell');rm(f);gc()
stopifnot(nrow(m[subcluster=='C5'])==716)
m[,analysis_group:=ifelse(subcluster=='C5','C5',ifelse(major_celltype=='Pericyte','Pericyte','Other_fibroblasts'))]
stopifnot(m[analysis_group=='Pericyte',.N]==5895,m[analysis_group=='Other_fibroblasts',.N]==13535)
manifest<-fread(file.path(OUT,'00_input_audit/manifest.csv'));maps<-list()
for(i in seq_len(nrow(manifest))){inp<-qs2::qs_read(manifest$input_path[i],nthreads=1);md<-as.data.table(inp$metadata);maps[[i]]<-data.table(cell=md$cell,specimen_id=manifest$specimen_id[i]);rm(inp);gc(FALSE)}
map<-rbindlist(maps);stopifnot(!anyDuplicated(map$cell));m<-merge(m[,.(cell,subcluster,major_celltype,analysis_group)],map,by='cell',all.x=TRUE);stopifnot(!anyNA(m$specimen_id))
writecsv(m,file.path(OUT,'00_input_audit/stromal_identity_final.csv.gz'))
counts<-dcast(m[analysis_group!='C5',.N,by=.(specimen_id,analysis_group)],specimen_id~analysis_group,value.var='N',fill=0)
writecsv(counts,file.path(OUT,'00_input_audit/stromal_specimen_counts.csv'));writecsv(m[,.N,by=.(subcluster,major_celltype,analysis_group)],file.path(OUT,'00_input_audit/stromal_identity_crosstab.csv'))
cat('IDENTITY AUDIT COMPLETE\n')
