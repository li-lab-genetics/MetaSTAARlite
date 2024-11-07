#' Performs meta-analysis of coding functional categories using MetaSTAAR
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
#' @param coding_sumstat_gene_list a list containing study-specific summary statistics corresponding to the specified gene.
#' @param coding_cov_gene_list a list containing study-specific sparse weighted covariance matrices corresponding to the specified gene.
#' @param cov_maf_cutoff a numeric vector with the length of \code{study.names}
#' indicating the maximum minor allele frequency cutoffs under which the sparse weighted
#' covariance files between variants are stored.
#' @param rare_maf_cutoff the cutoff of maximum minor allele frequency in
#' defining rare variants (default = 0.01).
#' @param rv_num_cutoff the cutoff of minimum number of variants of meta-analyzing
#' a given variant-set (default = 2).
#' @param check_qc_label a logical value indicating whether variants need to be dropped according to \code{qc_label}
#' specified in \code{\link{generate_MetaSTAAR_sumstat}} and \code{\link{generate_MetaSTAAR_cov}}. Default is FALSE.
#' @param variant_type type of variant included in the analysis. Choices include "SNV", "Indel", or "variant" (default = "SNV").
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

coding_MetaSTAARlite <- function(chr,gene_name,genes,
                                 sample.sizes,coding_sumstat_gene_list,coding_cov_gene_list,
                                 cov_maf_cutoff,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                 check_qc_label=FALSE,variant_type=c("SNV","Indel","variant"),
                                 Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,silent=FALSE){

  ## evaluate choices
  variant_type <- match.arg(variant_type)

  ### Gene
  kk <- which(genes[,1]==gene_name)

  ################################################
  #                  plof_ds
  ################################################
  results_plof_ds <- c()

  sumstat.list <- lapply(coding_sumstat_gene_list, function(x) {
    x[["plof_ds"]]
  })
  cov.list <- lapply(coding_cov_gene_list, function(x) {
    x[["plof_ds"]]
  })

  gene_test_merge <- MetaSTAARlite_merge(chr,sample.sizes,sumstat.list,cov.list,
                                         cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,
                                         check_qc_label=check_qc_label,variant_type=variant_type,
                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge$U)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR(gene_test_merge,annotation_phred,rv_num_cutoff),silent=silent)

    if(inherits(pvalues,"list"))
    {
      results_temp <- as.vector(genes[kk,])
      results_temp[3] <- "plof_ds"
      results_temp[2] <- chr
      results_temp[1] <- as.character(genes[kk,1])
      results_temp[4] <- pvalues$num_variant


      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)

      results_plof_ds <- rbind(results_plof_ds,results_temp)
    }
  }

  if(!is.null(results_plof_ds))
  {
    colnames(results_plof_ds) <- colnames(results_plof_ds, do.NULL = FALSE, prefix = "col")
    colnames(results_plof_ds)[1:4] <- c("Gene name","Chr","Category","#SNV")
    colnames(results_plof_ds)[(dim(results_plof_ds)[2]-1):dim(results_plof_ds)[2]] <- c("ACAT-O-MS","MetaSTAAR-O")
  }

  #####################################################
  #                      plof
  #####################################################
  results_plof <- c()

  sumstat.list <- lapply(coding_sumstat_gene_list, function(x) {
    x[["plof"]]
  })
  cov.list <- lapply(coding_cov_gene_list, function(x) {
    x[["plof"]]
  })

  gene_test_merge <- MetaSTAARlite_merge(chr,sample.sizes,sumstat.list,cov.list,
                                         cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,
                                         check_qc_label=check_qc_label,variant_type=variant_type,
                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge$U)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR(gene_test_merge,annotation_phred,rv_num_cutoff),silent=silent)

    if(inherits(pvalues,"list"))
    {
      results_temp <- as.vector(genes[kk,])
      results_temp[3] <- "plof"
      results_temp[2] <- chr
      results_temp[1] <- as.character(genes[kk,1])
      results_temp[4] <- pvalues$num_variant


      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)

      results_plof <- rbind(results_plof,results_temp)
    }
  }

  if(!is.null(results_plof))
  {
    colnames(results_plof) <- colnames(results_plof, do.NULL = FALSE, prefix = "col")
    colnames(results_plof)[1:4] <- c("Gene name","Chr","Category","#SNV")
    colnames(results_plof)[(dim(results_plof)[2]-1):dim(results_plof)[2]] <- c("ACAT-O-MS","MetaSTAAR-O")
  }

  #############################################
  #             synonymous
  #############################################
  results_synonymous <- c()

  sumstat.list <- lapply(coding_sumstat_gene_list, function(x) {
    x[["synonymous"]]
  })
  cov.list <- lapply(coding_cov_gene_list, function(x) {
    x[["synonymous"]]
  })

  gene_test_merge <- MetaSTAARlite_merge(chr,sample.sizes,sumstat.list,cov.list,
                                         cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,
                                         check_qc_label=check_qc_label,variant_type=variant_type,
                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge$U)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR(gene_test_merge,annotation_phred,rv_num_cutoff),silent=silent)

    if(inherits(pvalues,"list"))
    {
      results_temp <- as.vector(genes[kk,])
      results_temp[3] <- "synonymous"
      results_temp[2] <- chr
      results_temp[1] <- as.character(genes[kk,1])
      results_temp[4] <- pvalues$num_variant


      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)

      results_synonymous <- rbind(results_synonymous,results_temp)
    }
  }

  if(!is.null(results_synonymous))
  {
    colnames(results_synonymous) <- colnames(results_synonymous, do.NULL = FALSE, prefix = "col")
    colnames(results_synonymous)[1:4] <- c("Gene name","Chr","Category","#SNV")
    colnames(results_synonymous)[(dim(results_synonymous)[2]-1):dim(results_synonymous)[2]] <- c("ACAT-O-MS","MetaSTAAR-O")
  }

  #################################################
  #        missense
  #################################################
  results <- c()

  sumstat.list <- lapply(coding_sumstat_gene_list, function(x) {
    x[["missense"]]
  })
  cov.list <- lapply(coding_cov_gene_list, function(x) {
    x[["missense"]]
  })

  gene_test_merge <- MetaSTAARlite_merge(chr,sample.sizes,sumstat.list,cov.list,
                                         cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,
                                         check_qc_label=check_qc_label,variant_type=variant_type,
                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge$U)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR(gene_test_merge,annotation_phred,rv_num_cutoff),silent=silent)

    if(inherits(pvalues,"list"))
    {
      results_temp <- as.vector(genes[kk,])
      results_temp[3] <- "missense"
      results_temp[2] <- chr
      results_temp[1] <- as.character(genes[kk,1])
      results_temp[4] <- pvalues$num_variant


      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)

      results <- rbind(results,results_temp)
    }
  }

  #################################################
  #         disruptive missense
  #################################################
  sumstat.list <- lapply(coding_sumstat_gene_list, function(x) {
    x[["disruptive_missense"]]
  })
  cov.list <- lapply(coding_cov_gene_list, function(x) {
    x[["disruptive_missense"]]
  })

  gene_test_merge <- MetaSTAARlite_merge(chr,sample.sizes,sumstat.list,cov.list,
                                         cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,
                                         check_qc_label=check_qc_label,variant_type=variant_type,
                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(gene_test_merge$U)>=2)
  {
    ## Annotation
    annotation_phred <- gene_test_merge$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR(gene_test_merge,annotation_phred,rv_num_cutoff),silent=silent)

    if(inherits(pvalues,"list"))
    {
      results_temp <- as.vector(genes[kk,])
      results_temp[3] <- "disruptive_missense"
      results_temp[2] <- chr
      results_temp[1] <- as.character(genes[kk,1])
      results_temp[4] <- pvalues$num_variant


      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)

      results <- rbind(results,results_temp)
    }
  }

  if(!is.null(results))
  {
    colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
    colnames(results)[1:4] <- c("Gene name","Chr","Category","#SNV")
    colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O-MS","MetaSTAAR-O")

    if(dim(results)[1]==1)
    {
      if(results[3]!="disruptive_missense")
      {
        results <- cbind(results,matrix(1,1,6))
        colnames(results)[(dim(results)[2]-5):dim(results)[2]] <- c("SKAT-MS(1,25)-Disruptive","SKAT-MS(1,1)-Disruptive","Burden-MS(1,25)-Disruptive","Burden-MS(1,1)-Disruptive","ACAT-V-MS(1,25)-Disruptive","ACAT-V-MS(1,1)-Disruptive")
        results_missense <- results
        results_ds <- c()
      }else
      {
        results_missense <- c()
        results_ds <- results
        results <- c()
      }
    }

    if(!is.null(results))
    {
      if(dim(results)[1]==2)
      {
        results_m <- c(results[1,],rep(0,6))
        names(results_m)[(length(results_m)-5):length(results_m)] <- c("SKAT-MS(1,25)-Disruptive","SKAT-MS(1,1)-Disruptive","Burden-MS(1,25)-Disruptive","Burden-MS(1,1)-Disruptive","ACAT-V-MS(1,25)-Disruptive","ACAT-V-MS(1,1)-Disruptive")
        results_m[(length(results_m)-5):length(results_m)] <- results[2,c("SKAT-MS(1,25)","SKAT-MS(1,1)","Burden-MS(1,25)","Burden-MS(1,1)","ACAT-V-MS(1,25)","ACAT-V-MS(1,1)")]
        apc_num <- (length(results_m)-18)/6
        p_seq <- c(1:apc_num,1:apc_num+(apc_num+1),1:apc_num+2*(apc_num+1),1:apc_num+3*(apc_num+1),1:apc_num+4*(apc_num+1),1:apc_num+5*(apc_num+1),(6*apc_num+9):(6*apc_num+14))
        results_m["MetaSTAAR-O"] <- CCT(as.numeric(results_m[5:length(results_m)][p_seq]))
        results_m["MetaSTAAR-S(1,25)"] <- CCT(as.numeric(results_m[5:length(results_m)][c(1:apc_num,6*apc_num+9)]))
        results_m["MetaSTAAR-S(1,1)"] <- CCT(as.numeric(results_m[5:length(results_m)][c(1:apc_num+(apc_num+1),6*apc_num+10)]))
        results_m["MetaSTAAR-B(1,25)"] <- CCT(as.numeric(results_m[5:length(results_m)][c(1:apc_num+2*(apc_num+1),6*apc_num+11)]))
        results_m["MetaSTAAR-B(1,1)"] <- CCT(as.numeric(results_m[5:length(results_m)][c(1:apc_num+3*(apc_num+1),6*apc_num+12)]))
        results_m["MetaSTAAR-A(1,25)"] <- CCT(as.numeric(results_m[5:length(results_m)][c(1:apc_num+4*(apc_num+1),6*apc_num+13)]))
        results_m["MetaSTAAR-A(1,1)"] <- CCT(as.numeric(results_m[5:length(results_m)][c(1:apc_num+5*(apc_num+1),6*apc_num+14)]))

        results_ds <- c()
        results_ds <- rbind(results_ds,results[2,])

        results <- c()
        results <- rbind(results,results_m)
      }
    }
  }else
  {
    results <- c()
    results_ds <- c()
  }

  results_coding <- list(plof=results_plof,plof_ds=results_plof_ds,missense=results,disruptive_missense=results_ds,synonymous=results_synonymous)

  return(results_coding)
}
