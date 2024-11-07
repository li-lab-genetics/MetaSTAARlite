#' Performs meta-analysis of individual variants using using MetaSTAAR.
#'
#' This function performs meta-analysis to detect associations between a quantitative/dichotomous phenotype
#' and each individual variant in a genetic region by using score test.
#' @param sample.sizes a numeric vector with the length of \code{study.names}
#' indicating the sample size of each study.
#' @param sumstat.dir a character vector containing the directories of the study-specific summary statistics file folders.
#' @param mac_cutoff the cutoff of minimum minor allele count in
#' defining individual variants. Default is 20.
#' @param check_qc_label a logical value indicating whether variants need to be dropped according to \code{qc_label}. Default is FALSE.
#' specified in \code{\link{generate_MetaSTAAR_sumstat}} and \code{\link{generate_MetaSTAAR_cov}}.
#' If \code{check_qc_label} is FALSE, it is assumed that no variant will be dropped (default = FALSE).
#' @return a data frame containing the score test p-value and the estimated effect size of the alternative allele for each individual variant in the given genetic region.
#' @export


individual_analysis_MetaSTAARlite <- function(sample.sizes,sumstat.list,
                                              mac_cutoff=20,check_qc_label=FALSE){

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

  U.merge <- Reduce("+", lapply(sumstat.merge.list, function(x) {
    x <- x[cv.index,]
    return(x$U * (2 * (x$alt_AC == x$MAC) - 1))
  }))

  V.merge <- Reduce("+", lapply(sumstat.merge.list, function(x) {
    x <- x[cv.index,]
    return(x$V)
  }))

  p.merge <- pchisq(U.merge^2/V.merge,df=1,lower.tail=FALSE)
  logp.merge <- -pchisq(U.merge^2/V.merge,df=1,lower.tail=FALSE,log.p=TRUE)

  rm(list=setdiff(ls(), c("info","U.merge","V.merge","p.merge","logp.merge")))
  gc()

  results_temp <- data.frame(chr=info$chr,
                             pos=info$pos,
                             ref=info$ref,
                             alt=info$alt,
                             alt_AC=info$alt_AC,
                             MAC=info$MAC,
                             MAF=info$MAF,
                             N=info$N,
                             p=p.merge,
                             logp=logp.merge,
                             Score=U.merge,
                             Score_se=sqrt(V.merge),
                             Est=U.merge/V.merge,
                             Est_se=1/sqrt(V.merge))

  if(!is.null(results_temp))
  {
    results <- data.frame(CHR=results_temp$chr,POS=results_temp$pos,REF=results_temp$ref,ALT=results_temp$alt,
                          ALT_AF=results_temp$alt_AC/(2*results_temp$N),MAF=results_temp$MAF,N=results_temp$N,
                          pvalue=results_temp$p,pvalue_log10=results_temp$logp/log(10),
                          Score=results_temp$Score,Score_se=results_temp$Score_se,
                          Est=results_temp$Est,Est_se=results_temp$Est_se)
    results <- results[order(results$POS,results$REF,results$ALT),]
    row.names(results) <- NULL
  }else
  {
    results <- NULL
  }

  return(results)
}
