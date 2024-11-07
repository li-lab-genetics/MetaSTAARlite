#' Performs conditional meta-analysis of individual variants using using MetaSTAAR.
#'
#' This function performs meta-analysis to detect associations between a quantitative/dichotomous phenotype
#' and each individual variant in a genetic region by using score test.
#' @param sample.sizes a numeric vector with the length of \code{study.names}
#' indicating the sample size of each study.
#' @param sumstat.dir a character vector containing the directories of the study-specific summary statistics file folders.
#' @param covcond.list a character vector containing the names of variants to be conditioned on.
#' @param mac_cutoff the cutoff of minimum minor allele count in
#' defining individual variants. Default is 20.
#' @param check_qc_label a logical value indicating whether variants need to be dropped according to \code{qc_label}. Default is FALSE.
#' specified in \code{\link{generate_MetaSTAAR_sumstat}} and \code{\link{generate_MetaSTAAR_cov}}.
#' If \code{check_qc_label} is FALSE, it is assumed that no variant will be dropped (default = FALSE).
#' @return a data frame containing the score test p-value and the estimated effect size of the alternative allele for each individual variant in the given genetic region.
#' @export

individual_analysis_MetaSTAARlite_cond <- function(individual_results,sample.sizes,sumstat.list,covcond.list,
                                                   mac_cutoff=20,effect.cond=c("homogeneous","heterogeneous"),
                                                   check_qc_label=FALSE){

  ## evaluate choices
  effect.cond <- match.arg(effect.cond)

  cov_maf_cutoff <- rep(0.5 + 1e-16,length(sample.sizes))
  ### summary statistics
  sumstat.list <- lapply(sumstat.list, function(x) {
    if (!(is.null(x)) && !("qc_label" %in% colnames(x))){
      data.frame(x[,1:4],qc_label="PASS",x[,5:dim(x)[2]],stringsAsFactors = FALSE)
    }else{
      x
    }
  })
  sumstat.varid.list <- mapply(function(x,y) {
    x[(x$MAF<y)&(x$MAF>0)&(x$qc_label=="PASS"),1:4]
  }, x = sumstat.list, y = cov_maf_cutoff, SIMPLIFY = FALSE)
  sumstat.varid.merge <- do.call("rbind",sumstat.varid.list)
  sumstat.varid.nodup <- sumstat.varid.merge[!duplicated(sumstat.varid.merge),]
  if (is.null(sumstat.varid.nodup)) {
    return(NULL)
  }else if (dim(sumstat.varid.nodup)[1] == 0) {
    return(NULL)
  }
  sumstat.merge.list <- lapply(sumstat.list, function(x) {
    if (is.null(x)) {
      cbind(sumstat.varid.nodup,qc_label=NA,alt_AC=NA,MAC=NA,MAF=NA,N=NA,U=NA,V=NA,"1"=NA)
    }else {
      left_join(sumstat.varid.nodup,x,by=c("chr"="chr",
                                           "pos"="pos",
                                           "ref"="ref",
                                           "alt"="alt"))
    }
  })
  sumstat.merge.list <- mapply(function(x,y) {
    x[is.na(x[,"N"]),"N"] <- y
    x[is.na(x[,"qc_label"]),"qc_label"] <- "PASS"
    x[is.na(x)] <- 0
    return(x)
  }, x = sumstat.merge.list, y = sample.sizes, SIMPLIFY = FALSE)
  sumstat.varid.nodup$index <- 1:dim(sumstat.varid.nodup)[1]
  sumstat.index.list <- lapply(sumstat.varid.list, function(x) {
    if (is.null(x)) {
      integer(0)
    }else{
      left_join(x,sumstat.varid.nodup,by=c("chr"="chr",
                                           "pos"="pos",
                                           "ref"="ref",
                                           "alt"="alt"))$index
    }
  })
  rm(sumstat.list,sumstat.varid.list,sumstat.varid.merge)
  gc()

  ### select "common" variant based on the input cutoff
  alt_AC.merge <- as.integer(Reduce("+",lapply(sumstat.merge.list, function(x) {x$alt_AC})))
  N.merge.nonzero <- Reduce("+",lapply(sumstat.merge.list, function(x) {x$N * (x$MAC != 0)}))
  alt_AF.merge.nonzero <- alt_AC.merge / (2 * N.merge.nonzero)
  N.merge.zero <- Reduce("+",lapply(sumstat.merge.list, function(x) {x$N * (x$MAC == 0)}))
  alt_AC.merge <- alt_AC.merge + (alt_AF.merge.nonzero > 0.5) * (2 * N.merge.zero)
  N.merge <- Reduce("+",lapply(sumstat.merge.list, function(x) {x$N}))
  MAC.merge <- pmin(alt_AC.merge, 2 * N.merge - alt_AC.merge)
  MAF.merge <- MAC.merge / (2 * N.merge)
  cv.index <- (MAF.merge>(mac_cutoff - 0.5) / (2 * N.merge)) & Reduce("*",mapply(function(x,y) {
    x$MAF<y
  }, x = sumstat.merge.list, y = cov_maf_cutoff, SIMPLIFY = FALSE))
  if (check_qc_label){
    cv.index <- cv.index & Reduce("*",lapply(sumstat.merge.list, function(x) {
      x$qc_label=="PASS"
    }))
  }

  info <- cbind(sumstat.varid.nodup[,c("chr","pos","ref","alt")],
                alt_AC=alt_AC.merge,MAC=MAC.merge,MAF=MAF.merge,N=N.merge)[cv.index,]

  diagV.merge.list <- lapply(sumstat.merge.list, function(x) {
    return(diag(as.vector(x$V)))
  })

  ### covariance files for conditional analysis
  variant_info_list <- lapply(covcond.list, function(x) {x$variant_info})
  variant_info_merge <- do.call("rbind",variant_info_list)
  variant_info_nodup <- variant_info_merge[!duplicated(variant_info_merge),]
  variant_info_nodup$index <- 1:dim(variant_info_nodup)[1]
  variant_info_index_list <- lapply(variant_info_list, function(x) {
    if (is.null(x)) {
      integer(0)
    }else{
      left_join(x,variant_info_nodup,by=c("chr"="chr",
                                          "pos"="pos",
                                          "ref"="ref",
                                          "alt"="alt"))$index
    }
  })
  variant_adj_info_list <- lapply(covcond.list, function(x) {x$variant_adj_info[,1:4]})
  variant_adj_info_merge <- do.call("rbind",variant_adj_info_list)
  variant_adj_info_nodup <- variant_adj_info_merge[!duplicated(variant_adj_info_merge),]
  variant_adj_info_nodup$index <- 1:dim(variant_adj_info_nodup)[1]
  variant_adj_info_index_list <- lapply(variant_adj_info_list, function(x) {
    if (is.null(x)) {
      integer(0)
    }else{
      left_join(x,variant_adj_info_nodup,by=c("chr"="chr",
                                              "pos"="pos",
                                              "ref"="ref",
                                              "alt"="alt"))$index
    }
  })
  covcond.list <- mapply(function(x,y,z) {
    GTPG_cond_nodup <- matrix(0, nrow = dim(variant_info_nodup)[1], ncol = dim(variant_adj_info_nodup)[1])
    GTPG_cond_nodup[y,z] <- x$GTPG_cond
    variant_adj_info_nodup <- data.frame(variant_adj_info_nodup[,c("chr","pos","ref","alt")], U = 0, V = 0)
    variant_adj_info_nodup[z,] <- x$variant_adj_info
    return(list(GTPG_cond = GTPG_cond_nodup,
                G_condTPG_cond = x$G_condTPG_cond,
                variant_info = variant_info_nodup[,c("chr","pos","ref","alt")],
                variant_adj_info = variant_adj_info_nodup))
  },x = covcond.list, y = variant_info_index_list, z = variant_adj_info_index_list, SIMPLIFY = FALSE)
  rm(variant_info_list,variant_info_merge,variant_info_nodup,variant_info_index_list,
     variant_adj_info_list,variant_adj_info_merge,variant_adj_info_nodup,variant_adj_info_index_list)
  gc()

  variant_info <- covcond.list[[1]]$variant_info
  variant_adj_info <- covcond.list[[1]]$variant_adj_info[,c("chr","pos","ref","alt")]
  U_adj.list <- lapply(covcond.list, function(x) {x$variant_adj_info$U})
  GTPG_cond.list <- lapply(covcond.list, function(x) {x$GTPG_cond})
  G_condTPG_cond.list <- lapply(covcond.list, function(x) {x$G_condTPG_cond})
  info$index <- 1:dim(info)[1]
  ex.index <- left_join(variant_adj_info,info,by=c("chr"="chr",
                                                   "pos"="pos",
                                                   "ref"="ref",
                                                   "alt"="alt"))$index
  ex.index <- ex.index[!is.na(ex.index)]
  if (length(ex.index) > 0) {
    info <- info[-ex.index,!names(info)%in%c("index")]
  }else{
    info <- info[,!names(info)%in%c("index")]
  }

  variant_info$index <- 1:dim(variant_info)[1]
  var.common.index <- left_join(info,variant_info,by=c("chr"="chr",
                                                       "pos"="pos",
                                                       "ref"="ref",
                                                       "alt"="alt"))$index

  if (length(ex.index) > 0){
    if (effect.cond == "homogeneous") {
      U.common.merge <- Reduce("+", lapply(sumstat.merge.list, function(x) {
        x <- x[cv.index,][-ex.index,]
        return(x$U * (2 * (x$alt_AC == x$MAC) - 1))
      }))
      U.adj.merge <- Reduce("+", U_adj.list)
      cov.common.adj.merge <- Reduce("+", mapply(function(x,y) {
        x <- x[cv.index,][-ex.index,]
        return(y[var.common.index,,drop=FALSE] * (2 * (x$alt_AC == x$MAC) - 1))
      }, x = sumstat.merge.list, y = GTPG_cond.list, SIMPLIFY = FALSE))
      cov.adj.merge <- Reduce("+", G_condTPG_cond.list)
      cov.common.merge <- Reduce("+", lapply(diagV.merge.list, function(x) {x[cv.index, cv.index][-ex.index,-ex.index]}))

      U.merge <- U.common.merge-cov.common.adj.merge%*%ginv(cov.adj.merge)%*%U.adj.merge # Do not use solve()
      cov.merge <- cov.common.merge-cov.common.adj.merge%*%ginv(cov.adj.merge)%*%t(cov.common.adj.merge)
    }else{
      U.merge <- Reduce("+", mapply(function(x,y,z,z.adj) {
        x <- x[cv.index,][-ex.index,]
        u <- x$U * (2 * (x$alt_AC == x$MAC) - 1)
        z.common <- z[var.common.index,,drop=FALSE] * (2 * (x$alt_AC == x$MAC) - 1)
        return(u-z.common%*%ginv(z.adj)%*%y) # Do not use solve()
      }, x = sumstat.merge.list, y = U_adj.list, z = GTPG_cond.list, z.adj = G_condTPG_cond.list, SIMPLIFY = FALSE))

      cov.merge <- Reduce("+", mapply(function(x,y,z,z.adj) {
        cov <- x[cv.index, cv.index][-ex.index,-ex.index]
        y <- y[cv.index,][-ex.index,]
        z.common <- z[var.common.index,,drop=FALSE] * (2 * (y$alt_AC == y$MAC) - 1)
        return(cov-z.common%*%ginv(z.adj)%*%t(z.common))
      }, x = diagV.merge.list, y = sumstat.merge.list, z = GTPG_cond.list, z.adj = G_condTPG_cond.list, SIMPLIFY = FALSE))
    }
  }else{
    if (effect.cond == "homogeneous") {
      U.common.merge <- Reduce("+", lapply(sumstat.merge.list, function(x) {
        x <- x[cv.index,]
        return(x$U * (2 * (x$alt_AC == x$MAC) - 1))
      }))
      U.adj.merge <- Reduce("+", U_adj.list)
      cov.common.adj.merge <- Reduce("+", mapply(function(x,y) {
        x <- x[cv.index,]
        return(y[var.common.index,,drop=FALSE] * (2 * (x$alt_AC == x$MAC) - 1))
      }, x = sumstat.merge.list, y = GTPG_cond.list, SIMPLIFY = FALSE))
      cov.adj.merge <- Reduce("+", G_condTPG_cond.list)
      cov.common.merge <- Reduce("+", lapply(diagV.merge.list, function(x) {x[cv.index, cv.index]}))

      U.merge <- U.common.merge-cov.common.adj.merge%*%ginv(cov.adj.merge)%*%U.adj.merge # Do not use solve()
      cov.merge <- cov.common.merge-cov.common.adj.merge%*%ginv(cov.adj.merge)%*%t(cov.common.adj.merge)
    }else{
      U.merge <- Reduce("+", mapply(function(x,y,z,z.adj) {
        x <- x[cv.index,]
        u <- x$U * (2 * (x$alt_AC == x$MAC) - 1)
        z.common <- z[var.common.index,,drop=FALSE] * (2 * (x$alt_AC == x$MAC) - 1)
        return(u-z.common%*%ginv(z.adj)%*%y) # Do not use solve()
      }, x = sumstat.merge.list, y = U_adj.list, z = GTPG_cond.list, z.adj = G_condTPG_cond.list, SIMPLIFY = FALSE))

      cov.merge <- Reduce("+", mapply(function(x,y,z,z.adj) {
        cov <- x[cv.index, cv.index]
        y <- y[cv.index,]
        z.common <- z[var.common.index,,drop=FALSE] * (2 * (y$alt_AC == y$MAC) - 1)
        return(cov-z.common%*%ginv(z.adj)%*%t(z.common))
      }, x = diagV.merge.list, y = sumstat.merge.list, z = GTPG_cond.list, z.adj = G_condTPG_cond.list, SIMPLIFY = FALSE))
    }
  }

  rm(list=setdiff(ls(), c("info","U.merge","cov.merge")))
  gc()

  U.merge <- as.vector(U.merge)
  V.merge <- diag(as.matrix(cov.merge))
  p.merge <- pchisq(U.merge^2/V.merge,df=1,lower.tail=FALSE)
  logp.merge <- -pchisq(U.merge^2/V.merge,df=1,lower.tail=FALSE,log.p=TRUE)

  individual_results_cond <- data.frame(CHR=info$chr,
                                        POS=info$pos,
                                        REF=info$ref,
                                        ALT=info$alt,
                                        pvalue_cond=p.merge,
                                        pvalue_cond_log10=logp.merge/log(10))

  individual_results <- left_join(individual_results,individual_results_cond,by=c("CHR"="CHR",
                                                                                  "POS"="POS",
                                                                                  "REF"="REF",
                                                                                  "ALT"="ALT"))
  individual_results[is.na(individual_results$pvalue_cond),"pvalue_cond"] <- 1
  individual_results[is.na(individual_results$pvalue_cond_log10),"pvalue_cond_log10"] <- 0
  row.names(individual_results) <- NULL

  return(individual_results)
}
