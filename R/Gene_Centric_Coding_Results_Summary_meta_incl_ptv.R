#' Generates summary table and/or visualization for the meta-analysis of coding functional categories that was conducted.
#'
#' This function generates a summary table, Manhattan plot, and QQ plot for the meta-analysis of coding functional categories that was conducted
#' based on the parameters provided by the user.
#'
#' @param gene_centric_coding_jobs_num an integer which specifies the number of jobs done in the meta-analysis.
#' @param input_path a character which specifies the file path to the meta-analysis results files.
#' @param output_path a character which specifies the file path to the desired location of the produced summary table and visualizations.
#' @param gene_centric_results_name a character which specifies the name (excluding the jobs number) of the meta-analysis results files.
#' @param alpha a numeric value which specifies the desired significance threshold. Default is 2.5E-06.
#' @param manhattan_plot a logical value which determines if a Manhattan plot is generated. Default is FALSE.
#' @param QQ_plot a logical value which determines if a QQ plot is generated. Default is FALSE.

Gene_Centric_Coding_Results_Summary_meta_incl_ptv <- function(gene_centric_coding_jobs_num,input_path,output_path,gene_centric_results_name,
                                                              alpha=2.5E-06,manhattan_plot=FALSE,QQ_plot=FALSE){

  #######################################################
  #     summarize unconditional analysis results
  #######################################################

  coding_categories <- c("plof","plof_ds","missense","disruptive_missense","synonymous","ptv","ptv_ds")
  coding_category_labels <- c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous","PTV","PTV+D")
  coding_sig_columns <- c("Gene name","Chr","Category","#SNV","SKAT-MS(1,25)","Burden-MS(1,1)","ACAT-V-MS(1,25)","MetaSTAAR-O")

  results_coding_genome <- c()

  for(kk in 1:gene_centric_coding_jobs_num)
  {
    print(kk)
    results_coding <- get(load(paste0(input_path,gene_centric_results_name,"_",kk,".Rdata")))

    results_coding_genome <- c(results_coding_genome,results_coding)
  }

  results_coding_by_category <- setNames(vector("list",length(coding_categories)),coding_categories)

  for(kk in 1:length(results_coding_genome))
  {
    results <- results_coding_genome[[kk]]

    if(is.null(results)==FALSE)
    {
      category <- as.character(results[3])
      if(category %in% coding_categories)
      {
        results_coding_by_category[[category]] <- rbind(results_coding_by_category[[category]],results)
      }
    }

    if(kk%%1000==0)
    {
      print(kk)
    }
  }

  results_plof_genome <- results_coding_by_category[["plof"]]
  results_ptv_genome <- results_coding_by_category[["ptv"]]
  results_plof_ds_genome <- results_coding_by_category[["plof_ds"]]
  results_ptv_ds_genome <- results_coding_by_category[["ptv_ds"]]
  results_missense_genome <- results_coding_by_category[["missense"]]
  results_disruptive_missense_genome <- results_coding_by_category[["disruptive_missense"]]
  results_synonymous_genome <- results_coding_by_category[["synonymous"]]

  ###### whole-genome results
  # plof
  save(results_plof_genome,file=paste0(output_path,"plof.Rdata"))
  # ptv
  save(results_ptv_genome,file=paste0(output_path,"ptv.Rdata"))
  # plof + disruptive missense
  save(results_plof_ds_genome,file=paste0(output_path,"plof_ds.Rdata"))
  # ptv + disruptive missense
  save(results_ptv_ds_genome,file=paste0(output_path,"ptv_ds.Rdata"))
  # missense
  save(results_missense_genome,file=paste0(output_path,"missense.Rdata"))
  # disruptive missense
  save(results_disruptive_missense_genome,file=paste0(output_path,"disruptive_missense.Rdata"))
  # synonymous
  save(results_synonymous_genome,file=paste0(output_path,"synonymous.Rdata"))

  ###### significant results
  # plof
  plof_sig <- results_plof_genome[results_plof_genome[,"MetaSTAAR-O"]<alpha,,drop=FALSE]
  write.csv(plof_sig,file=paste0(output_path,"plof_sig.csv"))
  # ptv
  ptv_sig <- results_ptv_genome[results_ptv_genome[,"MetaSTAAR-O"]<alpha,,drop=FALSE]
  write.csv(ptv_sig,file=paste0(output_path,"ptv_sig.csv"))
  # plof_ds
  plof_ds_sig <- results_plof_ds_genome[results_plof_ds_genome[,"MetaSTAAR-O"]<alpha,,drop=FALSE]
  write.csv(plof_ds_sig,file=paste0(output_path,"plof_ds_sig.csv"))
  # ptv_ds
  ptv_ds_sig <- results_ptv_ds_genome[results_ptv_ds_genome[,"MetaSTAAR-O"]<alpha,,drop=FALSE]
  write.csv(ptv_ds_sig,file=paste0(output_path,"ptv_ds_sig.csv"))
  # missense
  missense_sig <- results_missense_genome[results_missense_genome[,"MetaSTAAR-O"]<alpha,,drop=FALSE]
  write.csv(missense_sig,file=paste0(output_path,"missense_sig.csv"))
  # synonymous
  synonymous_sig <- results_synonymous_genome[results_synonymous_genome[,"MetaSTAAR-O"]<alpha,,drop=FALSE]
  write.csv(synonymous_sig,file=paste0(output_path,"synonymous_sig.csv"))
  # disruptive_missense
  disruptive_missense_sig <- results_disruptive_missense_genome[results_disruptive_missense_genome[,"MetaSTAAR-O"]<alpha,,drop=FALSE]
  write.csv(disruptive_missense_sig,file=paste0(output_path,"disruptive_missense_sig.csv"))
  # coding results
  coding_sig_by_category <- list(plof=plof_sig,
                                 plof_ds=plof_ds_sig,
                                 missense=missense_sig,
                                 disruptive_missense=disruptive_missense_sig,
                                 synonymous=synonymous_sig,
                                 ptv=ptv_sig,
                                 ptv_ds=ptv_ds_sig)
  coding_sig <- do.call(rbind,lapply(coding_categories,function(category) {
    coding_sig_by_category[[category]][,coding_sig_columns]
  }))
  write.csv(coding_sig,file=paste0(output_path,"coding_sig.csv"))

  prepare_MetaSTAAR_results <- function(results_genome,category)
  {
    if(category=="missense")
    {
      results_MetaSTAAR <- results_genome[,c(1,2,dim(results_genome)[2]-6)]
    }else
    {
      results_MetaSTAAR <- results_genome[,c(1,2,dim(results_genome)[2])]
    }

    results_m <- c()
    for(i in 1:dim(results_MetaSTAAR)[2])
    {
      results_m <- cbind(results_m,unlist(results_MetaSTAAR[,i]))
    }

    colnames(results_m) <- colnames(results_MetaSTAAR)
    results_m <- data.frame(results_m,stringsAsFactors = FALSE)
    results_m[,2] <- as.numeric(results_m[,2])
    results_m[,3] <- as.numeric(results_m[,3])

    return(results_m)
  }

  if(manhattan_plot || QQ_plot)
  {
    genes_info_manhattan <- genes_info
    for(category in coding_categories)
    {
      results_m <- prepare_MetaSTAAR_results(results_coding_by_category[[category]],category)
      genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
      genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
      colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- category
    }

    ## ylim
    coding_minp <- min(genes_info_manhattan[,(dim(genes_info_manhattan)[2]-6):dim(genes_info_manhattan)[2]])
    min_y <- ceiling(-log10(coding_minp)) + 1
  }

  ## manhattan plot
  if(manhattan_plot)
  {
    pch <- 0:6
    figures <- vector("list",length(coding_categories))

    for(i in 1:length(coding_categories))
    {
      figures[[i]] <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan[,coding_categories[i]],sig.level=alpha,pch=pch[i],col = c("blue4", "orange3"),ylim=c(0,min_y),
                                     auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(coding_category_labels)))
    }

    print("Manhattan plot")

    png(paste0(output_path,"gene_centric_coding_manhattan.png"), width = 9, height = 6, units = 'in', res = 600)

    print(figures[[1]])
    for(i in 2:length(figures))
    {
      print(figures[[i]],newpage = FALSE)
    }

    dev.off()
  }

  ## Q-Q plot
  if(QQ_plot)
  {
    print("Q-Q plot")
    cex_point <- 1

    png(paste0(output_path,"gene_centric_coding_qqplot.png"), width = 8, height = 8, units = 'in', res = 600)

    for(i in 1:length(coding_categories))
    {
      ## remove unconverged p-values
      observed <- sort(genes_info_manhattan[,coding_categories[i]][genes_info_manhattan[,coding_categories[i]] < 1])
      lobs <- -(log10(observed))

      expected <- c(1:length(observed))
      lexp <- -(log10(expected / (length(expected)+1)))

      if(i>1)
      {
        par(new=T)
      }
      par(mar=c(5,6,4,4))
      plot(lexp,lobs,pch=i-1, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
           xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
           font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)

      abline(0, 1, col="red",lwd=2)
    }

    legend("topleft",legend=coding_category_labels,ncol=1,bty="o",box.lwd=1,pch=0:6,cex=1.5,text.font=2)

    dev.off()
  }

}

