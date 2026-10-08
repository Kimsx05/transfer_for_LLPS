# This file defines a gated function only. Sourcing it never runs ligand activity.
source('/home/jinshengxi/bladder_metabolism/cellchat_nichechat_fibro7/scripts/nichenet_common.R')

run_frozen_receiver_program <- function(sample, source, receiver, geneset_oi,
                                       background_expressed_genes,
                                       potential_ligands, program_definition,
                                       target_links_per_ligand, output_directory) {
  set.seed(123)
  required_definition<-c('program_label','biological_rationale','source_reference',
                         'selection_rule','sample_replication_design','frozen')
  if(!all(required_definition%in%names(program_definition)) ||
     !isTRUE(program_definition$frozen))
    stop('A documented, explicitly frozen receiver program is required; no geneset_oi is inferred here.')
  for(field in setdiff(required_definition,'frozen')) {
    value<-program_definition[[field]]
    if(!is.character(value) || length(value)!=1L || is.na(value) || !nzchar(trimws(value)))
      stop('Receiver-program definition must explicitly document: ',field)
  }
  if(!is.character(geneset_oi) || !length(geneset_oi) || anyNA(geneset_oi))
    stop('Supply an explicit receiver geneset_oi; there is no automatic default.')
  if(!is.character(background_expressed_genes) || !length(background_expressed_genes))
    stop('Supply a receiver-specific expressed background for this biological sample.')
  if(!source%in%GROUPS[1:7] || !receiver%in%GROUPS[8:14])
    stop('Invalid frozen sender or receiver identity.')
  cfg<-readRDS(file.path(ROOT,'01_metadata/frozen_analysis_inputs.rds'))
  sample_key<-sample;source_key<-source;receiver_key<-receiver
  eligibility<-cfg$eligibility[sample==sample_key & source==source_key & target==receiver_key]
  if(nrow(eligibility)!=1L || !isTRUE(eligibility$eligible))
    stop('This biological sample/source/receiver pair is not evaluable (both groups need >=30 cells).')
  if(!is.numeric(target_links_per_ligand) || length(target_links_per_ligand)!=1L || target_links_per_ligand<1)
    stop('Explicitly specify the target-link display/reporting limit.')
  nn<-readRDS(file.path(NN_OUT,'frozen_nichenet_preflight_inputs.rds'))
  lt<-readRDS(nn$resource_paths[['ligand_target']])
  if(any(!geneset_oi%in%background_expressed_genes) || any(!background_expressed_genes%in%rownames(lt)))
    stop('geneset_oi must lie within the expressed receiver background and verified model target universe; audit exclusions explicitly first.')
  if(length(setdiff(background_expressed_genes,geneset_oi))==0L)
    stop('Background must include genes outside the program so activity can be assessed.')
  if(!length(potential_ligands) || any(!potential_ligands%in%colnames(lt)))
    stop('Supply expressed candidate sender ligands with expressed receiver receptors and verified ligand-target model support.')
  if(dir.exists(output_directory) && length(list.files(output_directory)))
    stop('Refusing to overwrite an existing activity result directory.')
  dir.create(output_directory,recursive=TRUE,showWarnings=FALSE)
  suppressPackageStartupMessages(library(nichenetr))
  activities<-nichenetr::predict_ligand_activities(
    geneset=unique(geneset_oi),background_expressed_genes=unique(background_expressed_genes),
    ligand_target_matrix=lt,potential_ligands=unique(potential_ligands),single=TRUE)
  activities<-as.data.table(activities)
  activities[,`:=`(sample=sample,source=source,receiver=receiver,program=program_definition$program_label)]
  setorder(activities,-pearson,-aupr_corrected)
  # Preserve every ligand and all four activity metrics; no guessed significance cutoff.
  fwrite(activities,file.path(output_directory,'ligand_activity.tsv'),sep='\t',na='NA')
  links<-rbindlist(lapply(unique(potential_ligands),function(l)nichenetr::get_weighted_ligand_target_links(
    ligand=l,geneset=unique(geneset_oi),ligand_target_matrix=lt,n=target_links_per_ligand)),fill=TRUE)
  fwrite(links,file.path(output_directory,'ligand_target_links.tsv'),sep='\t',na='NA')
  saveRDS(list(sample=sample,source=source,receiver=receiver,program_definition=program_definition,
               geneset_oi=geneset_oi,background_expressed_genes=background_expressed_genes,
               potential_ligands=potential_ligands,target_links_per_ligand=target_links_per_ligand,
               seed=123,package_version=as.character(packageVersion('nichenetr'))),
          file.path(output_directory,'frozen_analysis_definition.rds'))
  invisible(list(activities=activities,ligand_target_links=links))
}
