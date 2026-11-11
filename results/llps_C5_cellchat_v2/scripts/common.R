set.seed(123)
BASE<-'/home/jinshengxi/bladder_metabolism';OLD<-file.path(BASE,'cellchat_nichechat_fibro7');OUT<-file.path(BASE,'results/llps_C5_cellchat_v2');AUD<-file.path(BASE,'results/llps_C5_distribution_audit_20261110_v1/v2')
.libPaths(unique(c(file.path(OLD,'Rlib'),'/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_atlas_stage1_v1/Rlib','/home/jinshengxi/Rlibrary/LLPS_pipeline_v1_20260914_150459','/home/jinshengxi/Rlibrary/packages',.libPaths())))
suppressPackageStartupMessages({library(Seurat);library(data.table);library(Matrix);library(CellChat);library(qs2);library(qs)})
data.table::setDTthreads(1);future::plan(future::sequential)
if(requireNamespace('RhpcBLASctl',quietly=TRUE)){RhpcBLASctl::blas_set_num_threads(1);RhpcBLASctl::omp_set_num_threads(1)}
stopifnot(packageVersion('CellChat')=='2.1.2')
GROUPS<-c(paste0('Fibro_C',0:6),'Epithelial','Endothelial','Myeloid','T_NK','B','Plasma','Mast')
MAJOR<-'/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes/results/objects/major_annotated_with_LLPS77_AUC.qs'
FIBRO<-'/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes/fibro_pericyte_subcluster_clean_v4_global_noC5/global_resolution_0.12_numeric_order/global_res012_object.qs'
TD<-function(t)file.path(OUT,c('01_threshold10','02_threshold20','03_threshold30')[match(t,c(10,20,30))])
writecsv<-function(x,p)fwrite(x,p,na='NA')
atomicRDS<-function(x,p){saveRDS(x,paste0(p,'.tmp'),compress='gzip');stopifnot(file.rename(paste0(p,'.tmp'),p))}
if(dir.exists(file.path(OUT,'00_audit'))){proc<-readLines('/proc/self/status');obs<-data.table(timestamp=as.character(Sys.time()),pid=Sys.getpid(),cpu_affinity=trimws(sub('Cpus_allowed_list:','',grep('^Cpus_allowed_list:',proc,value=TRUE))),OS_threads=as.integer(trimws(sub('Threads:','',grep('^Threads:',proc,value=TRUE)))),data_table_threads=data.table::getDTthreads(),future_workers=future::nbrOfWorkers(),OMP=Sys.getenv('OMP_NUM_THREADS'),OPENBLAS=Sys.getenv('OPENBLAS_NUM_THREADS'),MKL=Sys.getenv('MKL_NUM_THREADS'),BLIS=Sys.getenv('BLIS_NUM_THREADS'));fwrite(obs,file.path(OUT,'00_audit/resource_observations.csv'),append=file.exists(file.path(OUT,'00_audit/resource_observations.csv')))}

if(dir.exists(file.path(OUT,'00_audit'))){tdirs<-list.dirs('/proc/self/task',recursive=FALSE);details<-rbindlist(lapply(tdirs,function(d){st<-readLines(file.path(d,'status'));data.table(timestamp=as.character(Sys.time()),pid=Sys.getpid(),thread=basename(d),name=trimws(sub('Name:','',grep('^Name:',st,value=TRUE))),state=trimws(sub('State:','',grep('^State:',st,value=TRUE))),affinity=trimws(sub('Cpus_allowed_list:','',grep('^Cpus_allowed_list:',st,value=TRUE))),wait_channel=readLines(file.path(d,'wchan'),warn=FALSE))}));fwrite(details,file.path(OUT,'00_audit/resource_thread_details.csv'),append=file.exists(file.path(OUT,'00_audit/resource_thread_details.csv')))}
