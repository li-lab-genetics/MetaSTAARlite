#' Performs meta-analysis of noncoding functional categories using MetaSTAAR
#'
#' This function performs meta-analysis to detect associations between a
#' quantitative/dichotomous phenotype and coding functional categories of a gene by using the MetaSTAAR procedure.
#' For each coding functional category, the MetaSTAAR-O p-value is a p-value from an omnibus test
#' that aggregated SKAT-MS(1,25), SKAT-MS(1,1), Burden-MS(1,25), Burden-MS(1,1), ACAT-V-MS(1,25),
#' and ACAT-V-MS(1,1) together with p-values of each test weighted by each annotation
#' using Cauchy method.
#' @param chr an integer which specifies the chromosome number.
#' @param gene_name a character which specifies the name of the gene to be meta-analyzed using the MetaSTAAR procedure.
#' @param sample.sizes a numeric vector with the length of \code{study.names}
#' indicating the sample size of each study.
#' @param noncoding_sumstat_gene_list a list containing study-specific summary statistics corresponding to the specified gene.
#' @param noncoding_cov_gene_list a list containing study-specific sparse weighted covariance matrices corresponding to the specified gene.
#' @param noncoding_cov_cond_gene_list a list containing the list of variants to condition on
#' @param cov_maf_cutoff a numeric vector with the length of \code{study.names}
#' indicating the maximum minor allele frequency cutoffs under which the sparse weighted
#' covariance files between variants are stored.
#' @param rare_maf_cutoff the cutoff of maximum minor allele frequency in
#' defining rare variants (default = 0.01).
#' @param rv_num_cutoff the cutoff of minimum number of variants of meta-analyzing
#' a given variant-set (default = 2).
#' @param check_qc_label a logical value indicating whether variants need to be dropped according to \code{qc_label}
#' specified in \code{\link{generate_MetaSTAAR_sumstat}} and \code{\link{generate_MetaSTAAR_cov}}. Default is FALSE.
#' @param variant_type type of variant included in the analysis. Choices include "SNV", "Indel", or "variant".
#' Default is c("SNV","Indel","variant").
#' @param Use_annotation_weights use annotations as weights or not. Default is TRUE.
#' @param Annotation_name a vector of annotation names used in MetaSTAAR. Default is NULL.
#' @param silent a logical value which determines if the report of error messages will be suppressed. Default is FALSE.
#' @return a list of data frames containing the MetaSTAAR p-values (including MetaSTAAR-O) corresponding to the coding functional category of the given gene.
#' @references Li, X., et al. (2023). Powerful, scalable and resource-efficient
#' meta-analysis of rare variant associations in large whole genome sequencing studies.
#' \emph{Nature Genetics}, \emph{55}(1), 154-164.
#' (\href{https://doi.org/10.1038/s41588-022-01225-6}{pub})
#' @references Li, Z., Li, X., et al. (2022). A framework for detecting noncoding
#' rare-variant associations of large-scale whole-genome sequencing studies.
#' \emph{Nature Methods}.
#' (\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
#' @references Li, X., Li, Z., et al. (2020). Dynamic incorporation of multiple
#' in silico functional annotations empowers rare variant association analysis of
#' large whole-genome sequencing studies at scale. \emph{Nature Genetics}, \emph{52}(9), 969-983.
#' (\href{https://doi.org/10.1038/s41588-020-0676-4}{pub})
#' @export


noncoding_MetaSTAARlite_cond <- function(chr,gene_name,
                                         sample.sizes,noncoding_sumstat_gene_list,noncoding_cov_gene_list,noncoding_cov_cond_gene_list,
                                         cov_maf_cutoff,rare_maf_cutoff=0.01,rv_num_cutoff=2,effect.cond = c("homogeneous","heterogeneous"),
                                         check_qc_label=FALSE,variant_type=c("SNV","Indel","variant"),
                                         Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,silent=FALSE){

  ## evaluate choices
  effect.cond <- match.arg(effect.cond)
  variant_type <- match.arg(variant_type)

  ########################################
  #   Downstream

  results_downstream <- c()

  sumstat.list <- lapply(noncoding_sumstat_gene_list, function(x) {
    x[["downstream"]]
  })
  cov.list <- lapply(noncoding_cov_gene_list, function(x) {
    x[["downstream"]]
  })
  covcond.list <- lapply(noncoding_cov_cond_gene_list, function(x) {
    x[["downstream"]]
  })

  gene_test_merge_cond <- MetaSTAARlite_merge_cond(chr,sample.sizes,sumstat.list,cov.list,covcond.list,
                                                   cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,effect.cond=effect.cond,
                                                   check_qc_label=check_qc_label,variant_type=variant_type,
                                                   Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge_cond$U_cond)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge_cond$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR_cond(gene_test_merge_cond,annotation_phred,rv_num_cutoff),silent=silent)

    if(inherits(pvalues,"list"))
    {
      results_temp <- rep(NA,4)
      results_temp[3] <- "downstream_cond"
      results_temp[2] <- chr
      results_temp[1] <- as.character(gene_name)
      results_temp[4] <- pvalues$num_variant


      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)

      results_downstream <- rbind(results_downstream,results_temp)
    }
  }

  if(!is.null(results_downstream))
  {
    colnames(results_downstream) <- colnames(results_downstream, do.NULL = FALSE, prefix = "col")
    colnames(results_downstream)[1:4] <- c("Gene name","Chr","Category","#SNV")
    colnames(results_downstream)[(dim(results_downstream)[2]-1):dim(results_downstream)[2]] <- c("ACAT-O-MS","MetaSTAAR-O")
  }

  ########################################
  #   Upstream

  results_upstream <- c()

  sumstat.list <- lapply(noncoding_sumstat_gene_list, function(x) {
    x[["upstream"]]
  })
  cov.list <- lapply(noncoding_cov_gene_list, function(x) {
    x[["upstream"]]
  })
  covcond.list <- lapply(noncoding_cov_cond_gene_list, function(x) {
    x[["upstream"]]
  })

  gene_test_merge_cond <- MetaSTAARlite_merge_cond(chr,sample.sizes,sumstat.list,cov.list,covcond.list,
                                                   cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,effect.cond=effect.cond,
                                                   check_qc_label=check_qc_label,variant_type=variant_type,
                                                   Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge_cond$U_cond)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge_cond$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR_cond(gene_test_merge_cond,annotation_phred,rv_num_cutoff),silent=silent)

    if(inherits(pvalues,"list"))
    {
      results_temp <- rep(NA,4)
      results_temp[3] <- "upstream_cond"
      results_temp[2] <- chr
      results_temp[1] <- as.character(gene_name)
      results_temp[4] <- pvalues$num_variant


      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)

      results_upstream <- rbind(results_upstream,results_temp)
    }
  }

  if(!is.null(results_upstream))
  {
    colnames(results_upstream) <- colnames(results_upstream, do.NULL = FALSE, prefix = "col")
    colnames(results_upstream)[1:4] <- c("Gene name","Chr","Category","#SNV")
    colnames(results_upstream)[(dim(results_upstream)[2]-1):dim(results_upstream)[2]] <- c("ACAT-O-MS","MetaSTAAR-O")
  }

  ########################################################
  #                UTR

  results_UTR <- c()

  sumstat.list <- lapply(noncoding_sumstat_gene_list, function(x) {
    x[["UTR"]]
  })
  cov.list <- lapply(noncoding_cov_gene_list, function(x) {
    x[["UTR"]]
  })
  covcond.list <- lapply(noncoding_cov_cond_gene_list, function(x) {
    x[["UTR"]]
  })

  gene_test_merge_cond <- MetaSTAARlite_merge_cond(chr,sample.sizes,sumstat.list,cov.list,covcond.list,
                                                   cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,effect.cond=effect.cond,
                                                   check_qc_label=check_qc_label,variant_type=variant_type,
                                                   Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge_cond$U_cond)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge_cond$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR_cond(gene_test_merge_cond,annotation_phred,rv_num_cutoff),silent=silent)

    if(inherits(pvalues,"list"))
    {
      results_temp <- rep(NA,4)
      results_temp[3] <- "UTR_cond"
      results_temp[2] <- chr
      results_temp[1] <- as.character(gene_name)
      results_temp[4] <- pvalues$num_variant


      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)

      results_UTR <- rbind(results_UTR,results_temp)
    }
  }

  if(!is.null(results_UTR))
  {
    colnames(results_UTR) <- colnames(results_UTR, do.NULL = FALSE, prefix = "col")
    colnames(results_UTR)[1:4] <- c("Gene name","Chr","Category","#SNV")
    colnames(results_UTR)[(dim(results_UTR)[2]-1):dim(results_UTR)[2]] <- c("ACAT-O-MS","MetaSTAAR-O")
  }

  #############################################
  #   Promoter-CAGE

  results_promoter_CAGE <- c()

  sumstat.list <- lapply(noncoding_sumstat_gene_list, function(x) {
    x[["promoter_CAGE"]]
  })
  cov.list <- lapply(noncoding_cov_gene_list, function(x) {
    x[["promoter_CAGE"]]
  })
  covcond.list <- lapply(noncoding_cov_cond_gene_list, function(x) {
    x[["promoter_CAGE"]]
  })

  gene_test_merge_cond <- MetaSTAARlite_merge_cond(chr,sample.sizes,sumstat.list,cov.list,covcond.list,
                                                   cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,effect.cond=effect.cond,
                                                   check_qc_label=check_qc_label,variant_type=variant_type,
                                                   Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge_cond$U_cond)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge_cond$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR_cond(gene_test_merge_cond,annotation_phred,rv_num_cutoff),silent=silent)

    if(inherits(pvalues,"list"))
    {
      results_temp <- rep(NA,4)
      results_temp[3] <- "promoter_CAGE_cond"
      results_temp[2] <- chr
      results_temp[1] <- as.character(gene_name)
      results_temp[4] <- pvalues$num_variant


      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)

      results_promoter_CAGE <- rbind(results_promoter_CAGE,results_temp)
    }
  }

  if(!is.null(results_promoter_CAGE))
  {
    colnames(results_promoter_CAGE) <- colnames(results_promoter_CAGE, do.NULL = FALSE, prefix = "col")
    colnames(results_promoter_CAGE)[1:4] <- c("Gene name","Chr","Category","#SNV")
    colnames(results_promoter_CAGE)[(dim(results_promoter_CAGE)[2]-1):dim(results_promoter_CAGE)[2]] <- c("ACAT-O-MS","MetaSTAAR-O")
  }

  ##################################################
  #       Promoter-DHS

  results_promoter_DHS <- c()

  sumstat.list <- lapply(noncoding_sumstat_gene_list, function(x) {
    x[["promoter_DHS"]]
  })
  cov.list <- lapply(noncoding_cov_gene_list, function(x) {
    x[["promoter_DHS"]]
  })
  covcond.list <- lapply(noncoding_cov_cond_gene_list, function(x) {
    x[["promoter_DHS"]]
  })

  gene_test_merge_cond <- MetaSTAARlite_merge_cond(chr,sample.sizes,sumstat.list,cov.list,covcond.list,
                                                   cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,effect.cond=effect.cond,
                                                   check_qc_label=check_qc_label,variant_type=variant_type,
                                                   Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge_cond$U_cond)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge_cond$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR_cond(gene_test_merge_cond,annotation_phred,rv_num_cutoff),silent=silent)

    if(inherits(pvalues,"list"))
    {
      results_temp <- rep(NA,4)
      results_temp[3] <- "promoter_DHS_cond"
      results_temp[2] <- chr
      results_temp[1] <- as.character(gene_name)
      results_temp[4] <- pvalues$num_variant


      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)

      results_promoter_DHS <- rbind(results_promoter_DHS,results_temp)
    }
  }

  if(!is.null(results_promoter_DHS))
  {
    colnames(results_promoter_DHS) <- colnames(results_promoter_DHS, do.NULL = FALSE, prefix = "col")
    colnames(results_promoter_DHS)[1:4] <- c("Gene name","Chr","Category","#SNV")
    colnames(results_promoter_DHS)[(dim(results_promoter_DHS)[2]-1):dim(results_promoter_DHS)[2]] <- c("ACAT-O-MS","MetaSTAAR-O")
  }

  ###########################################
  #        Enhancer-CAGE

  results_enhancer_CAGE <- c()

  sumstat.list <- lapply(noncoding_sumstat_gene_list, function(x) {
    x[["enhancer_CAGE"]]
  })
  cov.list <- lapply(noncoding_cov_gene_list, function(x) {
    x[["enhancer_CAGE"]]
  })
  covcond.list <- lapply(noncoding_cov_cond_gene_list, function(x) {
    x[["enhancer_CAGE"]]
  })

  gene_test_merge_cond <- MetaSTAARlite_merge_cond(chr,sample.sizes,sumstat.list,cov.list,covcond.list,
                                                   cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,effect.cond=effect.cond,
                                                   check_qc_label=check_qc_label,variant_type=variant_type,
                                                   Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge_cond$U_cond)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge_cond$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR_cond(gene_test_merge_cond,annotation_phred,rv_num_cutoff),silent=silent)

    if(inherits(pvalues,"list"))
    {
      results_temp <- rep(NA,4)
      results_temp[3] <- "enhancer_CAGE_cond"
      results_temp[2] <- chr
      results_temp[1] <- as.character(gene_name)
      results_temp[4] <- pvalues$num_variant


      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)

      results_enhancer_CAGE <- rbind(results_enhancer_CAGE,results_temp)
    }
  }

  if(!is.null(results_enhancer_CAGE))
  {
    colnames(results_enhancer_CAGE) <- colnames(results_enhancer_CAGE, do.NULL = FALSE, prefix = "col")
    colnames(results_enhancer_CAGE)[1:4] <- c("Gene name","Chr","Category","#SNV")
    colnames(results_enhancer_CAGE)[(dim(results_enhancer_CAGE)[2]-1):dim(results_enhancer_CAGE)[2]] <- c("ACAT-O-MS","MetaSTAAR-O")
  }

  ##################################################
  #       Enhancer-DHS

  results_enhancer_DHS <- c()

  sumstat.list <- lapply(noncoding_sumstat_gene_list, function(x) {
    x[["enhancer_DHS"]]
  })
  cov.list <- lapply(noncoding_cov_gene_list, function(x) {
    x[["enhancer_DHS"]]
  })
  covcond.list <- lapply(noncoding_cov_cond_gene_list, function(x) {
    x[["enhancer_DHS"]]
  })

  gene_test_merge_cond <- MetaSTAARlite_merge_cond(chr,sample.sizes,sumstat.list,cov.list,covcond.list,
                                                   cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,effect.cond=effect.cond,
                                                   check_qc_label=check_qc_label,variant_type=variant_type,
                                                   Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge_cond$U_cond)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge_cond$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR_cond(gene_test_merge_cond,annotation_phred,rv_num_cutoff),silent=silent)

    if(inherits(pvalues,"list"))
    {
      results_temp <- rep(NA,4)
      results_temp[3] <- "enhancer_DHS_cond"
      results_temp[2] <- chr
      results_temp[1] <- as.character(gene_name)
      results_temp[4] <- pvalues$num_variant


      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)

      results_enhancer_DHS <- rbind(results_enhancer_DHS,results_temp)
    }
  }

  if(!is.null(results_enhancer_DHS))
  {
    colnames(results_enhancer_DHS) <- colnames(results_enhancer_DHS, do.NULL = FALSE, prefix = "col")
    colnames(results_enhancer_DHS)[1:4] <- c("Gene name","Chr","Category","#SNV")
    colnames(results_enhancer_DHS)[(dim(results_enhancer_DHS)[2]-1):dim(results_enhancer_DHS)[2]] <- c("ACAT-O-MS","MetaSTAAR-O")
  }

  ############################################
  #           results

  results_noncoding <- list(upstream=results_upstream,downstream=results_downstream,UTR=results_UTR,
                            promoter_CAGE=results_promoter_CAGE,promoter_DHS=results_promoter_DHS,
                            enhancer_CAGE=results_enhancer_CAGE,enhancer_DHS=results_enhancer_DHS)

  return(results_noncoding)
}
