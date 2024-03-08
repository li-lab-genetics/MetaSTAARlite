noncoding_MetaSTAARlite_worker <- function(chr,gene_name,genofile,obj_nullmodel,known_loci=NULL,
                                           cov_maf_cutoff=0.05,signif.digits=NULL,
                                           QC_label="annotation/filter",check_qc_label=FALSE,variant_type=c("SNV","Indel","variant"),
                                           Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                           Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                           silent=FALSE){

  ## evaluate choices
  variant_type <- match.arg(variant_type)

  phenotype.id <- as.character(obj_nullmodel$id_include)

  #####################################
  #   Gene Info
  ## get SNV id
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

  variant.id <- seqGetData(genofile, "variant.id")

  if(!is.null(known_loci))
  {
    position <- as.integer(seqGetData(genofile, "position"))
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

  ########################################
  #   Downstream

  GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
  is.in <- (GENCODE.Category=="downstream")&(SNVlist)
  variant.id.downstream <- variant.id[is.in]

  seqSetFilter(genofile,variant.id=variant.id.downstream)

  rm(variant.id.downstream)
  gc()

  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
  variant_gene_num <- sapply(GENCODE.Info.split,function(z) length(z))

  variant.id.SNV <- seqGetData(genofile, "variant.id")
  variant.id.SNV <- rep(variant.id.SNV,variant_gene_num)

  rm(GENCODE.Info)
  gc()

  rm(variant_gene_num)
  gc()

  Gene <- as.character(unlist(GENCODE.Info.split))

  rm(GENCODE.Info.split)
  gc()

  seqResetFilter(genofile)

  ### Gene
  is.in <- which(Gene==gene_name)
  variant.is.in <- variant.id.SNV[is.in]

  seqSetFilter(genofile,variant.id=variant.is.in)

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

    try(summary_stat <- MetaSTAARlite_worker_sumstat(genotype,obj_nullmodel,variant_info,qc_label,
                                                     Anno.Int.PHRED.sub),silent=silent)

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

  seqResetFilter(genofile)

  summary_stat_list[["downstream"]] <- summary_stat
  GTSinvG_rare_list[["downstream"]] <- GTSinvG_rare
  cov_cond_list[["downstream"]] <- cov_cond

  ########################################
  #   Upstream

  is.in <- (GENCODE.Category=="upstream")&(SNVlist)
  variant.id.upstream <- variant.id[is.in]

  seqSetFilter(genofile,variant.id=variant.id.upstream)

  rm(variant.id.upstream)
  gc()

  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
  variant_gene_num <- sapply(GENCODE.Info.split,function(z) length(z))

  variant.id.SNV <- seqGetData(genofile, "variant.id")
  variant.id.SNV <- rep(variant.id.SNV,variant_gene_num)

  rm(GENCODE.Info)
  gc()

  rm(variant_gene_num)
  gc()

  Gene <- as.character(unlist(GENCODE.Info.split))

  rm(GENCODE.Info.split)
  gc()

  seqResetFilter(genofile)

  ### Gene
  is.in <- which(Gene==gene_name)
  variant.is.in <- variant.id.SNV[is.in]

  seqSetFilter(genofile,variant.id=variant.is.in)

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

    try(summary_stat <- MetaSTAARlite_worker_sumstat(genotype,obj_nullmodel,variant_info,qc_label,
                                                     Anno.Int.PHRED.sub),silent=silent)

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

  seqResetFilter(genofile)

  summary_stat_list[["upstream"]] <- summary_stat
  GTSinvG_rare_list[["upstream"]] <- GTSinvG_rare
  cov_cond_list[["upstream"]] <- cov_cond

  ########################################################
  #                UTR

  is.in <- ((GENCODE.Category=="UTR3")|(GENCODE.Category=="UTR5")|(GENCODE.Category=="UTR5;UTR3"))&(SNVlist)
  variant.id.UTR <- variant.id[is.in]

  rm(GENCODE.Category)
  gc()

  seqSetFilter(genofile,variant.id=variant.id.UTR)

  rm(variant.id.UTR)
  gc()

  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[(]")

  rm(GENCODE.Info)
  gc()

  Gene <- as.character(sapply(GENCODE.Info.split,function(z) z[1]))

  rm(GENCODE.Info.split)
  gc()

  variant.id.SNV <- seqGetData(genofile, "variant.id")

  seqResetFilter(genofile)

  ### Gene
  is.in <- which(Gene==gene_name)
  variant.is.in <- variant.id.SNV[is.in]

  seqSetFilter(genofile,variant.id=variant.is.in)

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

    try(summary_stat <- MetaSTAARlite_worker_sumstat(genotype,obj_nullmodel,variant_info,qc_label,
                                                     Anno.Int.PHRED.sub),silent=silent)

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

  seqResetFilter(genofile)

  summary_stat_list[["UTR"]] <- summary_stat
  GTSinvG_rare_list[["UTR"]] <- GTSinvG_rare
  cov_cond_list[["UTR"]] <- cov_cond

  #############################################
  #   Promoter-CAGE

  ## Promoter
  varid <- seqGetData(genofile, "variant.id")
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  promGobj <- promoters(genes(txdb), upstream = 3000, downstream = 3000)

  # Subsetting Promoters that within +/-3kb of TSS and have CAGE signals
  CAGEAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="CAGE")]))
  CAGEBvt <- CAGEAnno!=""
  CAGEidx <- which(CAGEBvt,useNames=TRUE)
  seqSetFilter(genofile,variant.id=varid[CAGEidx])
  seqSetFilter(genofile,promGobj,intersect=TRUE)
  CAGEpromgene <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  CAGEGene <- unlist(lapply(strsplit(CAGEpromgene,"\\(|\\,|;|-"),`[[`,1))
  ##obtain variants info
  CAGEvchr <- as.numeric(seqGetData(genofile,"chromosome"))
  CAGEvpos <- as.numeric(seqGetData(genofile,"position"))
  CAGEvref <- as.character(seqGetData(genofile,"$ref"))
  CAGEvalt <- as.character(seqGetData(genofile,"$alt"))
  dfPromCAGEVarGene <- data.frame(CAGEvchr,CAGEvpos,CAGEvref,CAGEvalt,CAGEGene)

  ## get SNV id
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

  variant.id <- seqGetData(genofile, "variant.id")
  variant.id.SNV <- variant.id[SNVlist]

  dfPromCAGEVarGene.SNV <- dfPromCAGEVarGene[SNVlist,]
  dfPromCAGEVarGene.SNV$CAGEvpos <- as.character(dfPromCAGEVarGene.SNV$CAGEvpos)
  dfPromCAGEVarGene.SNV$CAGEvref <- as.character(dfPromCAGEVarGene.SNV$CAGEvref)
  dfPromCAGEVarGene.SNV$CAGEvalt <- as.character(dfPromCAGEVarGene.SNV$CAGEvalt)

  seqResetFilter(genofile)

  rm(dfPromCAGEVarGene)
  gc()

  ### Gene
  is.in <- which(dfPromCAGEVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV[is.in]

  seqSetFilter(genofile,variant.id=variant.is.in)

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

    try(summary_stat <- MetaSTAARlite_worker_sumstat(genotype,obj_nullmodel,variant_info,qc_label,
                                                     Anno.Int.PHRED.sub),silent=silent)

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

  seqResetFilter(genofile)

  summary_stat_list[["promoter_CAGE"]] <- summary_stat
  GTSinvG_rare_list[["promoter_CAGE"]] <- GTSinvG_rare
  cov_cond_list[["promoter_CAGE"]] <- cov_cond

  ##################################################
  #       Promoter-DHS

  # Subsetting Promoters that within +/-3kb of TSS and have rOCRs signals
  rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
  rOCRsBvt <- rOCRsAnno!=""
  rOCRsidx <- which(rOCRsBvt,useNames=TRUE)
  seqSetFilter(genofile,variant.id=varid[rOCRsidx])

  seqSetFilter(genofile,promGobj,intersect=TRUE)
  rOCRspromgene <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  rOCRsGene <- unlist(lapply(strsplit(rOCRspromgene,"\\(|\\,|;|-"),`[[`,1))
  ## obtain variants info
  rOCRsvchr <- as.numeric(seqGetData(genofile,"chromosome"))
  rOCRsvpos <- as.numeric(seqGetData(genofile,"position"))
  rOCRsvref <- as.character(seqGetData(genofile,"$ref"))
  rOCRsvalt <- as.character(seqGetData(genofile,"$alt"))
  dfPromrOCRsVarGene <- data.frame(rOCRsvchr,rOCRsvpos,rOCRsvref,rOCRsvalt,rOCRsGene)

  ## get SNV id
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

  variant.id <- seqGetData(genofile, "variant.id")
  variant.id.SNV <- variant.id[SNVlist]

  dfPromrOCRsVarGene.SNV <- dfPromrOCRsVarGene[SNVlist,]
  dfPromrOCRsVarGene.SNV$rOCRsvpos <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvpos)
  dfPromrOCRsVarGene.SNV$rOCRsvref <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvref)
  dfPromrOCRsVarGene.SNV$rOCRsvalt <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvalt)

  seqResetFilter(genofile)

  rm(dfPromrOCRsVarGene)
  gc()

  ### Gene
  is.in <- which(dfPromrOCRsVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV[is.in]

  seqSetFilter(genofile,variant.id=variant.is.in)

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

    try(summary_stat <- MetaSTAARlite_worker_sumstat(genotype,obj_nullmodel,variant_info,qc_label,
                                                     Anno.Int.PHRED.sub),silent=silent)

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

  seqResetFilter(genofile)

  summary_stat_list[["promoter_DHS"]] <- summary_stat
  GTSinvG_rare_list[["promoter_DHS"]] <- GTSinvG_rare
  cov_cond_list[["promoter_DHS"]] <- cov_cond

  ###########################################
  #        Enhancer-CAGE

  #Now extract the GeneHancer with CAGE Signal Overlay
  genehancerAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
  genehancer <- genehancerAnno!=""

  CAGEAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="CAGE")]))
  CAGE <- CAGEAnno!=""
  CAGEGeneHancervt <- CAGEAnno!=""&genehancerAnno!=""
  CAGEGeneHanceridx <- which(CAGEGeneHancervt,useNames=TRUE)
  seqSetFilter(genofile,variant.id=varid[CAGEGeneHanceridx])

  # variants that covered by whole GeneHancer without CAGE overlap.
  genehancerSet <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
  enhancerGene <- unlist(lapply(strsplit(genehancerSet,"="),`[[`,4))
  enhancer2GENE <- unlist(lapply(strsplit(enhancerGene,";"),`[[`,1))
  enhancervchr <- as.numeric(seqGetData(genofile,"chromosome"))
  enhancervpos <- as.numeric(seqGetData(genofile,"position"))
  enhancervref <- as.character(seqGetData(genofile,"$ref"))
  enhancervalt <- as.character(seqGetData(genofile,"$alt"))
  dfHancerCAGEVarGene <- data.frame(enhancervchr,enhancervpos,enhancervref,enhancervalt,enhancer2GENE)

  ## get SNV id
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

  variant.id <- seqGetData(genofile, "variant.id")
  variant.id.SNV <- variant.id[SNVlist]

  dfHancerCAGEVarGene.SNV <- dfHancerCAGEVarGene[SNVlist,]
  dfHancerCAGEVarGene.SNV$enhancervpos <- as.character(dfHancerCAGEVarGene.SNV$enhancervpos)
  dfHancerCAGEVarGene.SNV$enhancervref <- as.character(dfHancerCAGEVarGene.SNV$enhancervref)
  dfHancerCAGEVarGene.SNV$enhancervalt <- as.character(dfHancerCAGEVarGene.SNV$enhancervalt)

  seqResetFilter(genofile)

  rm(dfHancerCAGEVarGene)
  gc()

  ### Gene
  is.in <- which(dfHancerCAGEVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV[is.in]

  seqSetFilter(genofile,variant.id=variant.is.in)

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

    try(summary_stat <- MetaSTAARlite_worker_sumstat(genotype,obj_nullmodel,variant_info,qc_label,
                                                     Anno.Int.PHRED.sub),silent=silent)

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

  seqResetFilter(genofile)

  summary_stat_list[["enhancer_CAGE"]] <- summary_stat
  GTSinvG_rare_list[["enhancer_CAGE"]] <- GTSinvG_rare
  cov_cond_list[["enhancer_CAGE"]] <- cov_cond

  ##################################################
  #       Enhancer-DHS

  rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
  rOCRs <- rOCRsAnno!=""
  rOCRsGeneHancervt <- rOCRsAnno!=""&genehancerAnno!=""
  rOCRsGeneHanceridx <- which(rOCRsGeneHancervt,useNames=TRUE)
  seqSetFilter(genofile,variant.id=varid[rOCRsGeneHanceridx])
  # variants that covered by whole GeneHancer without rOCRs overlap.

  genehancerSet <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
  enhancerGene <- unlist(lapply(strsplit(genehancerSet,"="),`[[`,4))
  enhancer2GENE <- unlist(lapply(strsplit(enhancerGene,";"),`[[`,1))
  enhancervchr <- as.numeric(seqGetData(genofile,"chromosome"))
  enhancervpos <- as.numeric(seqGetData(genofile,"position"))
  enhancervref <- as.character(seqGetData(genofile,"$ref"))
  enhancervalt <- as.character(seqGetData(genofile,"$alt"))
  dfHancerrOCRsVarGene <- data.frame(enhancervchr,enhancervpos,enhancervref,enhancervalt,enhancer2GENE)

  rm(varid)
  gc()

  ## get SNV id
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

  variant.id <- seqGetData(genofile, "variant.id")
  variant.id.SNV <- variant.id[SNVlist]

  dfHancerrOCRsVarGene.SNV <- dfHancerrOCRsVarGene[SNVlist,]
  dfHancerrOCRsVarGene.SNV$enhancervpos <- as.character(dfHancerrOCRsVarGene.SNV$enhancervpos)
  dfHancerrOCRsVarGene.SNV$enhancervref <- as.character(dfHancerrOCRsVarGene.SNV$enhancervref)
  dfHancerrOCRsVarGene.SNV$enhancervalt <- as.character(dfHancerrOCRsVarGene.SNV$enhancervalt)

  seqResetFilter(genofile)

  rm(dfHancerrOCRsVarGene)
  gc()

  ### Gene
  is.in <- which(dfHancerrOCRsVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV[is.in]

  seqSetFilter(genofile,variant.id=variant.is.in)

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

    try(summary_stat <- MetaSTAARlite_worker_sumstat(genotype,obj_nullmodel,variant_info,qc_label,
                                                     Anno.Int.PHRED.sub),silent=silent)

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

  seqResetFilter(genofile)

  summary_stat_list[["enhancer_DHS"]] <- summary_stat
  GTSinvG_rare_list[["enhancer_DHS"]] <- GTSinvG_rare
  cov_cond_list[["enhancer_DHS"]] <- cov_cond

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
