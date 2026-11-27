set.seed(123);options(mc.cores=1)
BASE<-'/home/jinshengxi/bladder_metabolism';O<-file.path(BASE,'results/llps_C5_cellchat_visualization_v3');V<-file.path(BASE,'results/llps_C5_cellchat_v2');H<-file.path(BASE,'results/llps_C5_highlow_cellchat_ligand_v1');N1<-file.path(BASE,'results/C5_Endothelial_NicheNet_spatial_v1');N2<-file.path(BASE,'results/C5_Endothelial_GSEA_nominalNicheNet_v1')
.libPaths(unique(c(file.path(BASE,'cellchat_nichechat_fibro7/Rlib'),'/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_atlas_stage1_v1/Rlib','/home/jinshengxi/Rlibrary/LLPS_pipeline_v1_20260914_150459','/home/jinshengxi/Rlibrary/packages',.libPaths())))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(Matrix);library(qs2)})
setDTthreads(1);if(requireNamespace('RhpcBLASctl',quietly=TRUE)){RhpcBLASctl::blas_set_num_threads(1);RhpcBLASctl::omp_set_num_threads(1)}
skill<-'/home/jinshengxi/.agents/skills/jins-routine-plotting';for(f in c('theme.R','palettes.R','export.R'))source(file.path(skill,'common',f));jin_setup()
FIB<-paste0('Fibro_C',0:6);TME<-c('Epithelial','Endothelial','Myeloid','T_NK','B','Plasma','Mast');G<-c(FIB,TME);C5<-'Fibro_C5';FOCUS<-c('POSTN','FN1','COL1A1','COL1A2','COL3A1','COL6A2','COL6A3','HSPG2','MDK','ADM','SEMA3C','WNT2','VEGFA','INHBA','ANGPTL4','TGFB1','IGF2','TGFB2','WNT5A','TNFSF12')
cols<-c(setNames(JIN_CLUSTER_COLORS,FIB),Epithelial='#E59662',Endothelial='#637BA8',Myeloid='#B86F83',T_NK='#805B9B',B='#619DAE',Plasma='#A3998C',Mast='#AF8660');options(jin.plotting.colors=cols)
nice<-function(x)gsub('T_NK','T/NK',sub('Fibro_','',x),fixed=TRUE)
w<-function(x,p)fwrite(x,file.path(O,p),na='NA');r<-function(p)fread(file.path(O,p));q<-function(x,p=.5)if(all(is.na(x)))NA_real_ else as.numeric(quantile(x,p,na.rm=TRUE,names=FALSE))
TD<-function(t)file.path(V,c('01_threshold10','02_threshold20','03_threshold30')[match(t,c(10,20,30))])
clean<-function(s)gsub('[^A-Za-z0-9_]+','_',s)
cb<-function()guide_colourbar(display='rectangles')
expfig<-function(p,dir,name,data,width=9,height=6,role='Main figure candidates',source='llps_C5_cellchat_v2',threshold='20',unit='patient x clinical condition',purpose='',limitation='Edge-specific eligible patients; descriptive, not causal',script=STAGE){
 stem<-file.path(O,dir,name);fwrite(data,paste0(stem,'_data.csv'),na='NA');jin_export(p,stem,width,height)
 idx<-data.table(figure=name,path=paste0(dir,'/',name,'.png'),pdf_path=paste0(dir,'/',name,'.pdf'),data_path=paste0(dir,'/',name,'_data.csv'),script=paste0('scripts/',script),data_source=source,threshold=threshold,statistical_unit=unit,purpose=purpose,limitations=limitation,figure_class=role)
 ip<-file.path(O,'00_input_audit',paste0('index_',sub('[.]R$','',STAGE),'.csv'));fwrite(idx,ip,append=file.exists(ip));cat('FIGURE ',name,'\n');flush.console()
}
startstage<-function(){ip<-file.path(O,'00_input_audit',paste0('index_',sub('[.]R$','',STAGE),'.csv'));if(file.exists(ip))unlink(ip)}
net<-function(t=20)fread(file.path(V,'04_patient_level_summary',paste0('network_summary_t',t,'.csv')))[stratification=='tumor_only']
features<-function(kind,t=20)fread(file.path(V,'04_patient_level_summary',paste0('C5_',kind,'_summary_t',t,'.csv.gz')))
annotate_dir<-function(z){z[,`:=`(direction=ifelse(sender==C5,'outgoing','incoming'),partner=ifelse(sender==C5,receiver,sender))];z}
heat<-function(z,x,y,value,xorder=NULL,yorder=NULL,limit=NULL,label=FALSE,fillname='Patient median strength'){
 xcol<-x;ycol<-y;vcol<-value;z<-copy(z);z[,xx:=factor(get(xcol),levels=if(is.null(xorder))unique(get(xcol))else xorder)];z[,yy:=factor(get(ycol),levels=rev(if(is.null(yorder))unique(get(ycol))else yorder))];z[,v:=get(vcol)];p<-ggplot(z,aes(xx,yy,fill=v))+geom_tile(color='white',linewidth=.4)+scale_fill_viridis_c(limits=limit,na.value='grey80',guide=cb(),name=fillname)+jin_theme()+theme(axis.text.x=element_text(angle=40,hjust=1),axis.line=element_blank(),axis.ticks=element_blank())+labs(x=NULL,y=NULL)
 if(label)p<-p+geom_text(aes(label=ifelse(is.na(v),'NA',formatC(v,digits=2,format='g'))),size=3,family='Arial');p
}
