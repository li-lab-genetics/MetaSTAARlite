coding_MetaSTAARlite_worker <- function(chr,gene_name,genofile,obj_nullmodel,genes,known_loci=NULL,
                                        cov_maf_cutoff=0.05,signif.digits=NULL,
                                        QC_label="annotation/filter",check_qc_label=FALSE,variant_type=c("SNV","Indel","variant"),
                                        Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                        Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                        silent=FALSE){

  ## evaluate choices
  variant_type <- match.arg(variant_type)

  phenotype.id <- as.character(obj_nullmodel$id_include)

  filter <- seqGetData(genofile, QC_label)
  if(variant_type=="variant")
  {
    if(check_qc_label)
    {
      SNVlist <- TRUE
    }else
    {
      SNVlist <- filter == "PASS"
    }
  }

  if(variant_type=="SNV")
  {
    if(check_qc_label)
    {
      SNVlist <- isSNV(genofile)
    }else
    {
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }
  }

  if(variant_type=="Indel")
  {
    if(check_qc_label)
    {
      SNVlist <- !isSNV(genofile)
    }else
    {
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }
  }

  position <- as.numeric(seqGetData(genofile, "position"))
  variant.id <- seqGetData(genofile, "variant.id")

  if(!is.null(known_loci))
  {
    allele <- seqGetData(genofile, "allele")

    loc_SNV <- c()
    for (i in 1:dim(known_loci)[1]){
      loc_SNV <- c(loc_SNV,which((position==known_loci$POS[i])&(allele==paste0(known_loci$REF[i],",",known_loci$ALT[i]))))
    }

    seqSetFilter(genofile,variant.id=variant.id[loc_SNV],sample.id=phenotype.id)

    G_SNV <- seqGetData(genofile, "$dosage")

    if (!is.null(G_SNV)){
      id.SNV <- seqGetData(genofile,"sample.id")
      id.SNV.match <- rep(0,length(id.SNV))

      for(i in 1:length(id.SNV))
      {
        id.SNV.match[i] <- which.max(id.SNV==phenotype.id[i])
      }

      G_SNV <- G_SNV[id.SNV.match,,drop=FALSE]
      G_SNV_MAF <- apply(G_SNV, 2, function(x){sum(x[!is.na(x)])/2/sum(!is.na(x))})
      for (i in 1:length(G_SNV_MAF)){
        if (G_SNV_MAF[i] <= 0.5){
          G_SNV[is.na(G_SNV[,i]),i] <- 0
        }
        else{
          G_SNV[is.na(G_SNV[,i]),i] <- 2
        }
      }

      pos_adj <- as.integer(seqGetData(genofile, "position"))
      ref_adj <- as.character(seqGetData(genofile, "$ref"))
      alt_adj <- as.character(seqGetData(genofile, "$alt"))
      variant_adj_info <- data.frame(chr,pos_adj,ref_adj,alt_adj)
      colnames(variant_adj_info) <- c("chr","pos","ref","alt")
      variant_adj_info

      seqResetFilter(genofile)
    }
  }

  rm(filter)
  gc()

  ### Gene
  kk <- which(genes[,1]==gene_name)

  sub_start_loc <- genes[kk,3]
  sub_end_loc <- genes[kk,4]

  is.in <- (SNVlist)&(position>=sub_start_loc)&(position<=sub_end_loc)
  variant.id.gene <- variant.id[is.in]

  rm(position)
  gc()

  seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

  ## Gencode_Exonic
  GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
  ## Gencode
  GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
  ## Meta.SVM.Pred
  MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))

  ################################################
  #           Coding
  ################################################
  variant.id.gene <- seqGetData(genofile, "variant.id")
  lof.in.coding <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|(GENCODE.EXONIC.Category=="nonsynonymous SNV")|(GENCODE.EXONIC.Category=="synonymous SNV")
  variant.id.gene <- variant.id.gene[lof.in.coding]

  seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

  ## Gencode_Exonic
  GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
  ## Gencode
  GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
  ## Meta.SVM.Pred
  MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))

  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL

  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))

          if(Annotation_name[k]=="CADD")
          {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }

      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }

  ################################################
  #                  plof_ds
  ################################################
  variant.id.gene <- seqGetData(genofile, "variant.id")
  lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
  variant.id.gene.category <- variant.id.gene[lof.in.plof]

  seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

  pos <- as.integer(seqGetData(genofile, "position"))
  ref <- as.character(seqGetData(genofile, "$ref"))
  alt <- as.character(seqGetData(genofile, "$alt"))
  if(check_qc_label){
    qc_label <- as.character(seqGetData(genofile, QC_label))
  }else{
    qc_label <- NULL
  }

  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))

  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index

  ## Genotype
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match,,drop=FALSE]

  summary_stat_list <- list()
  GTSinvG_rare_list <- list()
  cov_cond_list <- list()

  summary_stat <- NULL
  GTSinvG_rare <- NULL
  cov_cond <- NULL

  if(!is.null(Geno))
  {
    # Summary statistics
    genotype <- matrix_impute(Geno)
    variant_info <- data.frame(chr,pos,ref,alt,row.names=NULL)

    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.plof,]

    try(summary_stat <- MetaSTAARlite_worker_sumstat(genotype,obj_nullmodel,variant_info,qc_label,
                                                     Anno.Int.PHRED.sub.category),silent=silent)

    # Covariance matrices
    genotype <- matrix_flip(Geno)$Geno
    genotype <- as(genotype,"dgCMatrix")

    rm(Geno)
    gc()

    try(GTSinvG_rare <- MetaSTAARlite_worker_cov(genotype,obj_nullmodel,cov_maf_cutoff,
                                                 qc_label,signif.digits),silent=silent)

    # Covariance matrices for conditional analysis
    if(!is.null(known_loci))
    {
      try(cov_cond <- MetaSTAARlite_worker_cov_cond(genotype,G_SNV,obj_nullmodel,variant_info,variant_adj_info),silent=silent)
    }
  }

  summary_stat_list[["plof_ds"]] <- summary_stat
  GTSinvG_rare_list[["plof_ds"]] <- GTSinvG_rare
  cov_cond_list[["plof_ds"]] <- cov_cond

  #####################################################
  #                      plof
  #####################################################
  lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")
  variant.id.gene.category <- variant.id.gene[lof.in.plof]

  seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

  pos <- as.integer(seqGetData(genofile, "position"))
  ref <- as.character(seqGetData(genofile, "$ref"))
  alt <- as.character(seqGetData(genofile, "$alt"))
  if(check_qc_label){
    qc_label <- as.character(seqGetData(genofile, QC_label))
  }else{
    qc_label <- NULL
  }

  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))

  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index

  ## Genotype
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match,,drop=FALSE]

  summary_stat <- NULL
  GTSinvG_rare <- NULL
  cov_cond <- NULL

  if(!is.null(Geno))
  {
    # Summary statistics
    genotype <- matrix_impute(Geno)
    variant_info <- data.frame(chr,pos,ref,alt,row.names=NULL)

    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.plof,]

    try(summary_stat <- MetaSTAARlite_worker_sumstat(genotype,obj_nullmodel,variant_info,qc_label,
                                                     Anno.Int.PHRED.sub.category),silent=silent)

    # Covariance matrices
    genotype <- matrix_flip(Geno)$Geno
    genotype <- as(genotype,"dgCMatrix")

    rm(Geno)
    gc()

    try(GTSinvG_rare <- MetaSTAARlite_worker_cov(genotype,obj_nullmodel,cov_maf_cutoff,
                                                 qc_label,signif.digits),silent=silent)

    # Covariance matrices for conditional analysis
    if(!is.null(known_loci))
    {
      try(cov_cond <- MetaSTAARlite_worker_cov_cond(genotype,G_SNV,obj_nullmodel,variant_info,variant_adj_info),silent=silent)
    }
  }

  summary_stat_list[["plof"]] <- summary_stat
  GTSinvG_rare_list[["plof"]] <- GTSinvG_rare
  cov_cond_list[["plof"]] <- cov_cond

  #############################################
  #             synonymous
  #############################################
  lof.in.synonymous <- (GENCODE.EXONIC.Category=="synonymous SNV")
  variant.id.gene.category <- variant.id.gene[lof.in.synonymous]

  seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

  pos <- as.integer(seqGetData(genofile, "position"))
  ref <- as.character(seqGetData(genofile, "$ref"))
  alt <- as.character(seqGetData(genofile, "$alt"))
  if(check_qc_label){
    qc_label <- as.character(seqGetData(genofile, QC_label))
  }else{
    qc_label <- NULL
  }

  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))

  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index

  ## Genotype
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match,,drop=FALSE]

  summary_stat <- NULL
  GTSinvG_rare <- NULL
  cov_cond <- NULL

  if(!is.null(Geno))
  {
    # Summary statistics
    genotype <- matrix_impute(Geno)
    variant_info <- data.frame(chr,pos,ref,alt,row.names=NULL)

    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.synonymous,]

    try(summary_stat <- MetaSTAARlite_worker_sumstat(genotype,obj_nullmodel,variant_info,qc_label,
                                                     Anno.Int.PHRED.sub.category),silent=silent)

    # Covariance matrices
    genotype <- matrix_flip(Geno)$Geno
    genotype <- as(genotype,"dgCMatrix")

    rm(Geno)
    gc()

    try(GTSinvG_rare <- MetaSTAARlite_worker_cov(genotype,obj_nullmodel,cov_maf_cutoff,
                                                 qc_label,signif.digits),silent=silent)

    # Covariance matrices for conditional analysis
    if(!is.null(known_loci))
    {
      try(cov_cond <- MetaSTAARlite_worker_cov_cond(genotype,G_SNV,obj_nullmodel,variant_info,variant_adj_info),silent=silent)
    }
  }

  summary_stat_list[["synonymous"]] <- summary_stat
  GTSinvG_rare_list[["synonymous"]] <- GTSinvG_rare
  cov_cond_list[["synonymous"]] <- cov_cond

  #################################################
  #        missense
  #################################################
  lof.in.missense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")
  variant.id.gene.category <- variant.id.gene[lof.in.missense]

  seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

  pos <- as.integer(seqGetData(genofile, "position"))
  ref <- as.character(seqGetData(genofile, "$ref"))
  alt <- as.character(seqGetData(genofile, "$alt"))
  if(check_qc_label){
    qc_label <- as.character(seqGetData(genofile, QC_label))
  }else{
    qc_label <- NULL
  }

  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))

  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index

  ## Genotype
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match,,drop=FALSE]

  summary_stat <- NULL
  GTSinvG_rare <- NULL
  cov_cond <- NULL

  if(!is.null(Geno))
  {
    # Summary statistics
    genotype <- matrix_impute(Geno)
    variant_info <- data.frame(chr,pos,ref,alt,row.names=NULL)

    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.missense,]

    try(summary_stat <- MetaSTAARlite_worker_sumstat(genotype,obj_nullmodel,variant_info,qc_label,
                                                     Anno.Int.PHRED.sub.category),silent=silent)

    # Covariance matrices
    genotype <- matrix_flip(Geno)$Geno
    genotype <- as(genotype,"dgCMatrix")

    rm(Geno)
    gc()

    try(GTSinvG_rare <- MetaSTAARlite_worker_cov(genotype,obj_nullmodel,cov_maf_cutoff,
                                                 qc_label,signif.digits),silent=silent)

    # Covariance matrices for conditional analysis
    if(!is.null(known_loci))
    {
      try(cov_cond <- MetaSTAARlite_worker_cov_cond(genotype,G_SNV,obj_nullmodel,variant_info,variant_adj_info),silent=silent)
    }
  }

  summary_stat_list[["missense"]] <- summary_stat
  GTSinvG_rare_list[["missense"]] <- GTSinvG_rare
  cov_cond_list[["missense"]] <- cov_cond

  #################################################
  #         disruptive missense
  #################################################
  lof.in.dmissense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D")
  variant.id.gene.category <- variant.id.gene[lof.in.dmissense]

  seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

  pos <- as.integer(seqGetData(genofile, "position"))
  ref <- as.character(seqGetData(genofile, "$ref"))
  alt <- as.character(seqGetData(genofile, "$alt"))
  if(check_qc_label){
    qc_label <- as.character(seqGetData(genofile, QC_label))
  }else{
    qc_label <- NULL
  }

  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))

  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index

  ## Genotype
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match,,drop=FALSE]

  summary_stat <- NULL
  GTSinvG_rare <- NULL
  cov_cond <- NULL

  if(!is.null(Geno))
  {
    # Summary statistics
    genotype <- matrix_impute(Geno)
    variant_info <- data.frame(chr,pos,ref,alt,row.names=NULL)

    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.dmissense,]

    try(summary_stat <- MetaSTAARlite_worker_sumstat(genotype,obj_nullmodel,variant_info,qc_label,
                                                     Anno.Int.PHRED.sub.category),silent=silent)

    # Covariance matrices
    genotype <- matrix_flip(Geno)$Geno
    genotype <- as(genotype,"dgCMatrix")

    rm(Geno)
    gc()

    try(GTSinvG_rare <- MetaSTAARlite_worker_cov(genotype,obj_nullmodel,cov_maf_cutoff,
                                                 qc_label,signif.digits),silent=silent)

    # Covariance matrices for conditional analysis
    if(!is.null(known_loci))
    {
      try(cov_cond <- MetaSTAARlite_worker_cov_cond(genotype,G_SNV,obj_nullmodel,variant_info,variant_adj_info),silent=silent)
    }
  }

  summary_stat_list[["disruptive_missense"]] <- summary_stat
  GTSinvG_rare_list[["disruptive_missense"]] <- GTSinvG_rare
  cov_cond_list[["disruptive_missense"]] <- cov_cond

  seqResetFilter(genofile)

  if(!is.null(known_loci))
  {
    return(list(summary_stat_list=summary_stat_list,
                GTSinvG_rare_list=GTSinvG_rare_list,
                cov_cond_list=cov_cond_list))
  }else
  {
    return(list(summary_stat_list=summary_stat_list,
                GTSinvG_rare_list=GTSinvG_rare_list))
  }
}
