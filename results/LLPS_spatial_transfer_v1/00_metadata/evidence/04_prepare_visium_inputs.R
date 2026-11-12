#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,warn=1)
Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",NUMEXPR_NUM_THREADS="1")
suppressPackageStartupMessages({library(data.table);library(Matrix);library(qs);library(SeuratObject);library(reticulate)})

project <- "/home/data/t070721/codex_workspace/Bladder Metabolism"
root <- file.path(project,"cell2location_preflight_v1"); out <- file.path(root,"02_visium"); adir <- file.path(out,"anndata")
dir.create(adir,recursive=TRUE,showWarnings=FALSE);dir.create(file.path(root,"logs"),recursive=TRUE,showWarnings=FALSE)
source_qs <- file.path(project,"ST/all_ST_combined_list_SeuratObject_updated.qs")
qc_path <- file.path(project,"ST/audit/spatial_final_QC_manifest.csv")
inv_path <- file.path(project,"ST/audit/spatial_dataset_inventory.csv")
python_path <- file.path(project,"cell2location_setup_v1/cell2location_env/bin/python")
ref_gene_path <- file.path(root,"01_reference/source_per_gene_totals.tsv.gz")
logcon<-file(file.path(root,"logs/04_prepare_visium_inputs.log"),"wt");sink(logcon,type="output",split=TRUE);sink(logcon,type="message")
on.exit({sink(type="message");sink(type="output");close(logcon)},add=TRUE)
stamp<-function(...)cat(sprintf("[%s] %s\n",format(Sys.time(),"%FT%TZ",tz="UTC"),paste(...,collapse=" ")))
safe<-function(x)gsub("[^A-Za-z0-9._-]+","_",x)

qc<-fread(qc_path);inv<-fread(inv_path)
if(nrow(qc)!=73L || sum(qc$eligible_for_RCTD)!=63L)stop("Frozen 73/63 manifest invariant failed")
if(anyDuplicated(qc$sample))stop("Duplicate section IDs")
ref_genes<-fread(ref_gene_path,select="gene_id")$gene_id

stamp("LOAD_73_SECTION_SOURCE",source_qs)
sp<-qread(source_qs,nthreads=10)
if(!is.list(sp)||length(sp)!=73L)stop("Combined source is not a 73-section list")

get_obj<-function(sm){
  hit<-inv[inv$sample==sm,object_name]
  if(length(hit)!=1L||!hit%in%names(sp))stop("Cannot uniquely map object: ",sm)
  sp[[hit]]
}
get_coords<-function(o,cells){
  for(im in Images(o)){
    z<-tryCatch(GetTissueCoordinates(o[[im]]),error=function(e)NULL)
    if(!is.null(z)&&all(cells%in%rownames(z))){
      pref<-list(c("imagecol","imagerow"),c("x","y"),c("col","row"))
      for(p in pref)if(all(p%in%colnames(z)))return(list(xy=as.matrix(z[cells,p,drop=FALSE]),image=im,columns=p))
      nums<-colnames(z)[vapply(z,is.numeric,logical(1))]
      if(length(nums)>=2L)return(list(xy=as.matrix(z[cells,nums[1:2],drop=FALSE]),image=im,columns=nums[1:2]))
    }
  }
  NULL
}
artifact_keep<-function(md){
  hit<-colnames(md)[grepl("artifact|fold|acellular|off.?tissue",colnames(md),ignore.case=TRUE)]
  keep<-rep(TRUE,nrow(md))
  for(f in hit){x<-tolower(as.character(md[[f]]));bad<-x%in%c("1","true","yes","artifact","fold","acellular","off_tissue","off-tissue");keep<-keep&!(!is.na(bad)&bad)}
  keep
}

# HARD GATE: deterministically reconstruct the already-defined hard-QC spot set.
recon<-vector("list",nrow(qc))
for(i in seq_len(nrow(qc))){
  d<-qc[i];o<-get_obj(d$sample);counts<-LayerData(o[["Spatial"]],layer="counts");md<-o[[]][colnames(counts),,drop=FALSE]
  keep<-Matrix::colSums(counts>0)>=500 & artifact_keep(md)
  # Frozen audit states all stored spots are tissue-associated; no new tissue inference is made.
  recon[[i]]<-data.table(section_id=d$sample,reconstructed=sum(keep),expected=d$n_hardQC_spots,identical=sum(keep)==d$n_hardQC_spots)
  rm(o,counts,md);gc(FALSE)
}
recon<-rbindlist(recon);fwrite(recon,file.path(out,"hardQC_spot_reconstruction_check.tsv"),sep="\t")
if(any(!recon$identical))stop("HARD GATE FAIL: reconstructed hard-QC spot counts differ from frozen manifest")

use_python(python_path,required=TRUE);ad<-import("anndata",convert=FALSE)
audits<-list();overlaps<-list();manifest<-vector("list",nrow(qc))
for(i in seq_len(nrow(qc))){
  d<-qc[i];eligible<-isTRUE(d$eligible_for_RCTD);o<-get_obj(d$sample)
  counts0<-LayerData(o[["Spatial"]],layer="counts");keep<-Matrix::colSums(counts0>0)>=500 & artifact_keep(o[[]][colnames(counts0),,drop=FALSE])
  counts<-counts0[,keep,drop=FALSE];cells<-colnames(counts);ci<-get_coords(o,cells)
  coord_ok<-!is.null(ci)&&nrow(ci$xy)==length(cells)&&identical(rownames(ci$xy),cells)
  images<-Images(o);image_ok<-length(images)>0L&&all(vapply(images,function(z)length(o[[z]]@image)>0L,logical(1)))
  sf_ok<-length(images)>0L&&all(vapply(images,function(z)length(o[[z]]@scale.factors)>0L,logical(1)))
  vals<-counts@x; totals<-Matrix::colSums(counts);features<-Matrix::colSums(counts>0)
  numeric_ok<-!anyNA(vals)&&all(is.finite(vals))&&!any(vals<0)&&!any(abs(vals-round(vals))>1e-8)
  status<-if(eligible&&coord_ok&&numeric_ok&&!anyDuplicated(rownames(counts))&&!anyDuplicated(cells))"PASS" else if(!eligible)"NOT_ELIGIBLE_FROZEN_QC" else "FAIL"
  h5<-if(eligible)file.path(adir,paste0(safe(d$sample),".h5ad")) else ""
  manifest[[i]]<-data.table(section_id=d$sample,dataset=d$dataset,specimen=d$sample,patient_uid="",clinical_group="",
    input_path=h5,counts_path=paste0(source_qs,"::",d$sample,"::Spatial/counts"),
    spatial_path=paste0(source_qs,"::",d$sample,"::",if(is.null(ci))"NO_COORDINATES" else ci$image),
    image_path=if(image_ok)paste0(source_qs,"::",d$sample,"::EMBEDDED_VISIUM_IMAGE") else "",
    eligible_for_cell2location=eligible,exclusion_reason=if(eligible)"" else paste0("FROZEN_SPATIAL_QC: ",d$final_QC_status,"; ",d$final_QC_reason),
    source_qc_status=d$final_QC_status,source_checkpoint_recorded=d$final_checkpoint_path,
    source_checkpoint_currently_exists=file.exists(d$final_checkpoint_path))
  audits[[i]]<-data.table(section_id=d$sample,dataset=d$dataset,status=status,n_spots=ncol(counts),n_genes=nrow(counts),
    n_in_tissue_spots=ncol(counts),total_counts=as.numeric(sum(vals)),median_counts_per_spot=median(totals),
    median_genes_per_spot=median(features),zero_count_spots=sum(totals==0),duplicate_genes=anyDuplicated(rownames(counts)),
    duplicate_barcodes=anyDuplicated(cells),coordinate_availability=coord_ok,image_availability=image_ok,
    scale_factor_availability=sf_ok,counts_nonnegative=!any(vals<0),counts_integer_like=!any(abs(vals-round(vals))>1e-8),
    no_NA=!anyNA(vals),no_Inf=all(is.finite(vals)),counts_source="Spatial/counts from combined Seurat; existing hard-QC rule reconstructed",
    input_matrix_scope="hard-QC retained tissue spots (not Space Ranger raw grid)")
  sg<-rownames(counts);ov<-intersect(ref_genes,sg)
  overlaps[[i]]<-data.table(section_id=d$sample,dataset=d$dataset,n_reference_genes=length(ref_genes),n_spatial_genes=length(sg),
    n_overlap_genes=length(ov),reference_overlap_fraction=length(ov)/length(ref_genes),spatial_overlap_fraction=length(ov)/length(sg),
    reference_ensembl_fraction=mean(grepl("^ENS[A-Z]*G[0-9]+",ref_genes)),spatial_ensembl_fraction=mean(grepl("^ENS[A-Z]*G[0-9]+",sg)),
    exact_case_overlap=length(ov),uppercase_overlap=length(intersect(toupper(ref_genes),toupper(sg))),
    gene_name_harmonization="NONE_EXACT_NAMES_USED",status=if(length(ov)/length(sg)>=0.5)"PASS" else "REVIEW")
  if(eligible){
    if(status!="PASS")stop("Eligible section input gate failed: ",d$sample)
    if(file.exists(h5)||file.exists(paste0(h5,".incomplete")))stop("Refusing overwrite: ",h5)
    obs<-data.frame(section_id=rep(d$sample,length(cells)),specimen=rep(d$sample,length(cells)),
                    dataset=rep(d$dataset,length(cells)),in_tissue=rep(TRUE,length(cells)),
                    source_qc_status=rep(d$final_QC_status,length(cells)),row.names=cells,check.names=FALSE)
    var<-data.frame(feature_name=rownames(counts),gene_symbol=rownames(counts),row.names=rownames(counts),check.names=FALSE)
    pyx<-r_to_py(as(Matrix::t(counts),"RsparseMatrix"),convert=FALSE)$astype("int32")
    aa<-ad$AnnData(X=pyx,obs=r_to_py(obs,convert=FALSE),var=r_to_py(var,convert=FALSE))
    aa$obsm$`__setitem__`("spatial",r_to_py(unname(ci$xy),convert=FALSE))
    sf<-as.list(o[[ci$image]]@scale.factors)
    aa$uns$`__setitem__`("input_provenance",r_to_py(list(counts_location="X",source_file=source_qs,
      source_assay="Spatial",source_layer="counts",matrix_scope="existing hard-QC retained spots",
      section_id=d$sample,coordinate_columns=as.list(ci$columns),embedded_image_available=image_ok,
      scale_factors=sf),convert=FALSE))
    tmp<-paste0(h5,".incomplete");aa$write_h5ad(tmp,compression="lzf");if(!file.rename(tmp,h5))stop("Atomic rename failed: ",h5)
    stamp("WROTE",d$sample,ncol(counts),nrow(counts),h5)
  }
  rm(o,counts0,counts,aa);gc(FALSE)
}
manifest<-rbindlist(manifest);audits<-rbindlist(audits);overlaps<-rbindlist(overlaps)
fwrite(manifest,file.path(out,"visium_section_manifest.tsv"),sep="\t")
fwrite(audits,file.path(out,"visium_input_audit.tsv"),sep="\t")
fwrite(overlaps,file.path(out,"reference_visium_gene_overlap.tsv"),sep="\t")
writeLines(c("# Standard 10x Visium input audit","",
  paste0("- Frozen QC sections: ",nrow(qc),"; eligible: ",sum(qc$eligible_for_RCTD),"."),
  "- Technology: standard 10x Visium (Seurat image class VisiumV1); no highres_ST/Visium HD source was used.",
  paste0("- Existing combined source: `",source_qs,"`."),
  "- Space Ranger matrix/spatial directories were not present; raw counts were taken from the preserved Seurat `Spatial/counts` layer.",
  "- Deleted final-QC checkpoints were reconstructed deterministically with the frozen hard-QC rule and required to match every manifest spot count.",
  "- Each eligible section was saved as a separate h5ad with raw counts in X and coordinates in obsm['spatial']; embedded H&E pixels were not duplicated.",
  paste0("- PASS: ",sum(audits$status=="PASS"),"; failed eligible: ",sum(audits$status=="FAIL"),".")),
  file.path(out,"visium_input_audit.md"))
