#' Perform the worker step of MetaSTAARlite for ncRNA masks
#'
#' This function uses MetaSTAARliteWorker to generate variant summary statistics
#' and sparse LD matrices for gene-centric coding analysis.
#'
#' @param chr an integer which specifies the chromosome number.
#' @param gene_name a character which specifies the name of the gene to be meta-analyzed using
#'the MetaSTAAR procedure.
#' @param genofile
#' @param obj_nullmodel an object from fitting the null model, which is the
#' output from either \code{\link{fit_null_glm}} function for unrelated samples or
#' \code{\link{fit_null_glmmkin}} function for related samples in the \code{\link{STAAR}} package.
#' @param genes
#' @param known_loci the data frame of variants to be adjusted for in conditional analysis. Should
#' contain four columns in the following order: chromosome (CHR), position (POS), reference allele (REF),
#' and alternative allele (ALT). Default is NULL.
#' @param cov_maf_cutoff a numeric value indicating the maximum minor allele frequency cutoff
#' under which the sparse weighted covariance file between variants is stored.
#' @param signif.digits an integer specifying the number of digits to be included beyond the
#' decimal point.
#' @param QC_label a character specifying the channel name of the QC label in the GDS/aGDS file.
#' Default is "annotation/filter".
#' @param check_qc_label #' @param check_qc_label a logical value indicating whether variants need to be dropped according to \code{qc_label}
#' specified in \code{\link{generate_MetaSTAAR_sumstat}} and \code{\link{generate_MetaSTAAR_cov}}. Default is FALSE.
#' @param variant_type a character vector specifying the types of variants to be considered.
#' Default is c("SNV","Indel","variant").
#' @param Annotation_dir a character specifying the channel name of the annotations in the aGDS file.
#' Default is "annotation/info/FunctionalAnnotation"
#' @param Annotation_name_catalog a data frame containing the annotation name and the corresponding
#' channel name in the aGDS file.
#' @param Use_annotation_weights a logical value which specifies if annotations will be used as weights
#' or not. Default is TRUE.
#' @param Annotation_name a vector of annotation names used in MetaSTAAR. Default is NULL.
#' @param silent a logical value which determines if the report of error messages will be suppressed. Default is FALSE.
#' @return two objects. First, the data frame of all variants in the variant-set (the summary statistics file),
#' including the following information: chromosome (chr), position (pos), reference allele (ref),
#' alternative allele (alt), quality control status (qc_label, optional), alternative allele count (alt_AC), minor allele count (MAC),
#' minor allele frequency (MAF), study sample size (N), score statistic (U), variance (V), and
#' the (low-rank decomposed) dense component of the covariance file. Second, the sparse matrix of all variants in the variant-set
#' whose minor allele frequency is below \code{cov_maf_cutoff} (the sparse weighted
#' covariance file), stored as a rectangle format.
#' @export
ncRNA_MetaSTAARlite_worker <- function(chr,gene_name,genofile,obj_nullmodel,known_loci=NULL,
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

  ## ncRNA SNVs
  GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
  is.in <- ((GENCODE.Category=="ncRNA_exonic")|(GENCODE.Category=="ncRNA_exonic;splicing")|(GENCODE.Category=="ncRNA_splicing"))&(SNVlist)

  variant.id.ncRNA <- variant.id[is.in]

  rm(GENCODE.Category)
  gc()

  seqSetFilter(genofile,variant.id=variant.id.ncRNA,sample.id=phenotype.id)

  rm(variant.id.ncRNA)
  gc()

  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[;]")
  Gene <- as.character(sapply(GENCODE.Info.split,function(z) gsub("\\(.*\\)","",z[1])))

  Gene_list_1 <- as.character(sapply(strsplit(Gene,','),'[',1))
  Gene_list_2 <- as.character(sapply(strsplit(Gene,','),'[',2))
  Gene_list_3 <- as.character(sapply(strsplit(Gene,','),'[',3))

  rm(GENCODE.Info)
  gc()

  rm(GENCODE.Info.split)
  gc()

  variant.id.ncRNA <- seqGetData(genofile, "variant.id")

  seqResetFilter(genofile)

  ### Gene
  is.in <- union(which(Gene_list_1==gene_name),which(Gene_list_2==gene_name))
  is.in <- union(is.in,which(Gene_list_3==gene_name))

  variant.is.in <- variant.id.ncRNA[is.in]

  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

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

  if(!is.null(known_loci))
  {
    return(list(summary_stat=summary_stat,
                GTSinvG_rare=GTSinvG_rare,
                cov_cond=cov_cond))
  }else
  {
    return(list(summary_stat=summary_stat,
                GTSinvG_rare=GTSinvG_rare))
  }
}
