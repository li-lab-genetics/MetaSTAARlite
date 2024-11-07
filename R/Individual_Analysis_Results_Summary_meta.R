#' Generates summary table and/or visualization for the meta-analysis of individual variants that was conducted.
#'
#' This function generates a summary table, Manhattan plot, and QQ plot for the meta-analysis of individual variants that was conducted
#' based on the parameters provided by the user.
#'
#' @param jobs_num an integer which specifies the number of jobs done in the meta-analysis.
#' @param input_path a character which specifies the file path to the meta-analysis results files.
#' @param output_path a character which specifies the file path to the desired location of the produced summary table and visualizations.
#' @param indivdual_results_name a character which specifies the name (excluding the jobs number) of the meta-analysis results files.
#' @param alpha a numeric value which specifies the desired significance threshold. Default is 5E-09.
#' @param manhattan_plot a logical value which determines if a Manhattan plot is generated. Default is FALSE.
#' @param QQ_plot a logical value which determines if a QQ plot is generated. Default is FALSE.

Individual_Analysis_Results_Summary_meta <- function(jobs_num,input_path,output_path,individual_results_name,
                                                     alpha=5E-09,manhattan_plot=FALSE,QQ_plot=FALSE){

  use_SPA <- FALSE
  ## Summarize Individual Analysis Results
  results_individual_analysis_genome <- c()
  num <- 0
  for(chr in 1:22)
  {
    print(chr)

    if(chr > 1)
    {
      num <- num + jobs_num$individual_analysis_num[chr-1]
    }

    results_individual_analysis_chr <- c()
    for(kk in 1:jobs_num$individual_analysis_num[chr])
    {
      print(kk)
      job_id <- kk + num
      results_individual_analysis <- get(load(paste0(input_path,individual_results_name,"_",job_id,".Rdata")))

      results_individual_analysis_chr <- rbind(results_individual_analysis_chr,results_individual_analysis)
    }
    results_individual_analysis_genome <- rbind(results_individual_analysis_genome,results_individual_analysis_chr)

    rm(results_individual_analysis_chr)
  }

  # save results
  save(results_individual_analysis_genome,file=paste0(output_path,"results_individual_analysis_genome.Rdata"))
  ## Significant findings
  results_sig <- results_individual_analysis_genome[results_individual_analysis_genome$pvalue<alpha,]

  # save significant results
  save(results_sig,file=paste0(output_path,"results_individual_analysis_sig.Rdata"))
  write.csv(results_sig,paste0(output_path,"results_individual_analysis_sig.csv"))

  ## manhattan plot
  if(manhattan_plot)
  {
    png(paste0(output_path,"manhattan_MAC_20.png"), width = 9, height = 6, units = 'in', res = 600)

    pvalue <- results_individual_analysis_genome$pvalue

    if(min(pvalue)==0)
    {
      if(!use_SPA)
      {
        print(manhattan_plot(results_individual_analysis_genome$CHR, results_individual_analysis_genome$POS, results_individual_analysis_genome$pvalue_log10, use_logp=TRUE, col = c("blue4", "orange3"),sig.level=alpha))
      }else
      {
        pvalue_log10 <- -log10(pvalue)
        pvalue_log10[!is.finite(pvalue_log10)] <- 308

        print(manhattan_plot(results_individual_analysis_genome$CHR, results_individual_analysis_genome$POS, pvalue_log10, use_logp=TRUE, col = c("blue4", "orange3"),sig.level=alpha))

        rm(pvalue_log10)
      }

    }else
    {
      print(manhattan_plot(results_individual_analysis_genome$CHR, results_individual_analysis_genome$POS, results_individual_analysis_genome$pvalue, col = c("blue4", "orange3"),sig.level=alpha))
    }

    rm(pvalue)
    gc()

    dev.off()
  }

  ## Q-Q plot
  if(QQ_plot)
  {
    observed <- sort(results_individual_analysis_genome$pvalue)

    if(min(observed)==0)
    {
      if(!use_SPA)
      {
        lobs <- sort(results_individual_analysis_genome$pvalue_log10,decreasing = TRUE)
      }else
      {
        lobs <- -(log10(observed))
        lobs[!is.finite(lobs)] <- 308
      }
    }else
    {
      lobs <- -(log10(observed))
    }

    expected <- c(1:length(observed))
    lexp <- -(log10(expected / (length(expected)+1)))

    rm(results_individual_analysis_genome)
    gc()

    png(paste0(output_path,"qqplot_MAC_20.png"), width = 9, height = 9, units = 'in', res = 600)

    par(mar=c(5,6,4,4))
    plot(lexp,lobs,pch=20, cex=1, xlim = c(0, max(lexp)), ylim = c(0, max(lobs)),
         xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
         font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)

    abline(0, 1, col="red",lwd=2)

    dev.off()
  }

}
