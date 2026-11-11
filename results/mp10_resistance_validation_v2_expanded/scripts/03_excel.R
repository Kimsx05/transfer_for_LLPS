.libPaths(c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/refdir/Rlib',.libPaths()))
library(readxl);library(data.table);setDTthreads(1)
B='/home/data/t070721/codex_workspace/Bladder Metabolism/results/mp10_resistance_validation_v2_expanded';S=file.path(B,'00_plan_audit/sources');D=file.path(S,'excel_tsv');dir.create(D,showWarnings=FALSE)
for(f in list.files(S,pattern='\\.(xlsx|xls)$',full.names=TRUE)){
 tryCatch({ss=excel_sheets(f);cat('\nFILE',basename(f),'sheets',ss,'\n');for(i in seq_along(ss)){x=read_excel(f,sheet=i,col_names=FALSE,.name_repair='unique_quiet');fwrite(x,file.path(D,paste0(basename(f),'__',i,'.tsv')),sep='\t',na='NA');cat('SHEET',i,'nrow',nrow(x),'ncol',ncol(x),'\n');print(as.data.table(x)[1:min(5,nrow(x)),1:min(12,ncol(x)),with=FALSE])}},error=function(e)cat('ERROR',f,conditionMessage(e),'\n'))
}
