source('/home/jinshengxi/bladder_metabolism/cellchat_nichechat_fibro7/scripts/common.R')
suppressPackageStartupMessages(library(ggplot2))
cons<-readRDS(file.path(ROOT,'03_cellchat_consensus/consensus_network.rds'))
net<-cons$filtered
stopifnot(identical(dimnames(net),list(GROUPS,GROUPS)),all(is.finite(net)),all(net>=0))
colors<-setNames(c('#247A6D','#48A999','#77BD9F','#9AC78D','#5D7E45','#A9AE54','#D0B564','#E59662','#637BA8','#B86F83','#805B9B','#619DAE','#A3998C','#AF8660'),GROUPS)
circle<-function(m,title) {
  # The complete 14-node frame is retained even for isolated nodes.
  layout(matrix(c(1,2),nrow=1),widths=c(4.4,1.6));par(mar=c(1,1,3,1),family='sans')
  plot(0,0,type='n',xlim=c(-1.48,1.48),ylim=c(-1.28,1.28),asp=1,axes=FALSE,xlab='',ylab='',main=title)
  theta<-pi/2-2*pi*(0:13)/14;xy<-cbind(cos(theta),sin(theta));maximum<-max(m)
  edges<-which(m>0,arr.ind=TRUE)
  if(nrow(edges))edges<-edges[order(m[edges]),,drop=FALSE]
  t<-seq(0,1,length.out=120)
  for(e in seq_len(nrow(edges))) {
    i<-edges[e,1];j<-edges[e,2];w<-m[i,j];width<-8*w/maximum;col<-adjustcolor(colors[i],alpha.f=.4)
    if(i==j) {
      radial<-xy[i,];perp<-c(-radial[2],radial[1]);start<-radial+.035*perp;end<-radial-.035*perp
      c1<-radial*1.19+.11*perp;c2<-radial*1.19-.11*perp
      curve<-outer((1-t)^3,start)+outer(3*(1-t)^2*t,c1)+outer(3*(1-t)*t^2,c2)+outer(t^3,end)
    } else {
      start<-xy[i,];end<-xy[j,];delta<-end-start;len<-sqrt(sum(delta^2));start<-start+.036*delta/len;end<-end-.036*delta/len
      perp<-c(-delta[2],delta[1]);control<-(start+end)*.28+.12*perp
      curve<-outer((1-t)^2,start)+outer(2*(1-t)*t,control)+outer(t^2,end)
    }
    lines(curve,lwd=width,col=col)
    k<-nrow(curve);arrows(curve[k-5,1],curve[k-5,2],curve[k,1],curve[k,2],length=.045,angle=22,lwd=width,col=col)
  }
  points(xy,pch=21,bg=colors,col='white',cex=2.6,lwd=1)
  for(i in 1:14)text(xy[i,1]*1.07,xy[i,2]*1.07,GROUPS[i],cex=.87,adj=c(if(abs(xy[i,1])<.01).5 else if(xy[i,1]>0)0 else 1,.5))
  par(mar=c(1,0,3,0));plot.new();plot.window(xlim=c(0,1),ylim=c(0,1))
  legend('topleft',legend=GROUPS,pch=21,pt.bg=colors,col='white',pt.cex=1.6,cex=.82,bty='n',y.intersp=1.13,title='Cell population')
  if(maximum>0) {
    values<-maximum*c(.25,.5,1)
    legend(x=0,y=.24,legend=format(values,digits=3),lwd=8*values/maximum,col='#555555',bty='n',cex=.78,title='Consensus strength')
  } else text(.03,.20,'No positive median edges',adj=0,cex=.75)
  text(.03,.045,'Equal node sizes\nArrows: source to target',adj=c(0,0),cex=.7)
  layout(1)
}
export_base<-function(stem,draw,width=12,height=9) {
 cairo_pdf(paste0(stem,'.pdf'),width=width,height=height,onefile=TRUE);draw();dev.off()
 png(paste0(stem,'.png'),width=width,height=height,units='in',res=450,type='cairo');draw();dev.off()
 svg(paste0(stem,'.svg'),width=width,height=height,onefile=TRUE);draw();dev.off()
}
export_base(file.path(ROOT,'04_cellchat_figures/Figure_A_cohort_consensus_circle'),function()circle(net,'Cohort CellChat consensus'))
fibro<-net;mask<-outer(grepl('^Fibro_C',GROUPS),grepl('^Fibro_C',GROUPS),'|');fibro[!mask]<-0
export_base(file.path(ROOT,'04_cellchat_figures/Figure_A_fibro_centered_circle'),function()circle(fibro,'Fibro-centered CellChat consensus'))
target<-GROUPS[8:14]
out<-cons$statistics[source%in%GROUPS[1:7]&target%in%GROUPS[8:14],.(fibro=source,major=target,strength=display_strength,direction='Outgoing: fibro to major')]
inc<-cons$statistics[target%in%GROUPS[1:7]&source%in%GROUPS[8:14],.(fibro=target,major=source,strength=display_strength,direction='Incoming: major to fibro')]
h<-rbind(out,inc);h[,fibro:=factor(fibro,levels=rev(GROUPS[1:7]))];h[,major:=factor(major,levels=GROUPS[8:14])]
p<-ggplot(h,aes(major,fibro,fill=strength))+geom_tile(width=.94,height=.94)+facet_wrap(~direction,nrow=1)+scale_fill_gradient(low='#F7F7F7',high='#A73D55',name='Consensus\nstrength')+scale_x_discrete(drop=FALSE)+scale_y_discrete(drop=FALSE)+labs(x=NULL,y=NULL,caption='Zero-inclusive sample median; n eligible >=3, prevalence >=0.30. Blank/zero display tiles can also reflect insufficient coverage; see audit table.')+theme_minimal(base_size=12)+theme(panel.grid=element_blank(),axis.text.x=element_text(angle=45,hjust=1),strip.text=element_text(face='bold'),plot.caption=element_text(size=8,hjust=0),legend.position='right')
stem<-file.path(ROOT,'09_figures/Figure_B_fibro_major_consensus_heatmap')
ggsave(paste0(stem,'.pdf'),p,width=12,height=5.7,device=cairo_pdf,bg='white')
ggsave(paste0(stem,'.png'),p,width=12,height=5.7,dpi=450,device=ragg::agg_png,bg='white')
ggsave(paste0(stem,'.svg'),p,width=12,height=5.7,device=grDevices::svg,bg='white')
writeLines(c('Figure A: full cohort circle is primary; fibro-centered circle is supplemental.','Both circle plots contain 14 equally sized nodes, fixed order, no edge labels, source-colored directed edges, and legend.','Edge line widths are exactly proportional to the displayed median communication strength. All retained positive edges are drawn; no aesthetic cutoff.','Figure B: outgoing and incoming fibro x seven major-population heatmaps; consistent linear strength color scale; no gridlines.','PDF and SVG are vector outputs. PNG is 450 dpi.','NicheNet ligand-activity and mechanism figures are deferred until a receiver geneset_oi is explicitly designed. NICHECHAT_STATUS=NOT_APPLICABLE.'),file.path(ROOT,'09_figures/figure_definition_audit.txt'))
cat('FIGURES COMPLETE\n')
