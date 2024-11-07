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

Gene_Centric_Coding_Results_Summary_meta <- function(gene_centric_coding_jobs_num,input_path,output_path,gene_centric_results_name,
                                                     alpha=2.5E-06,manhattan_plot=FALSE,QQ_plot=FALSE){

  #######################################################
  #     summarize unconditional analysis results
  #######################################################

  results_coding_genome <- c()

  for(kk in 1:gene_centric_coding_jobs_num)
  {
    print(kk)
    results_coding <- get(load(paste0(input_path,gene_centric_results_name,"_",kk,".Rdata")))

    results_coding_genome <- c(results_coding_genome,results_coding)
  }

  results_plof_genome <- c()
  results_plof_ds_genome <- c()
  results_missense_genome <- c()
  results_disruptive_missense_genome <- c()
  results_synonymous_genome <- c()

  for(kk in 1:length(results_coding_genome))
  {
    results <- results_coding_genome[[kk]]

    if(is.null(results)==FALSE)
    {
      ### plof
      if(results[3]=="plof")
      {
        results_plof_genome <- rbind(results_plof_genome,results)
      }
      ### plof_ds
      if(results[3]=="plof_ds")
      {
        results_plof_ds_genome <- rbind(results_plof_ds_genome,results)
      }
      ### missense
      if(results[3]=="missense")
      {
        results_missense_genome <- rbind(results_missense_genome,results)
      }
      ### disruptive_missense
      if(results[3]=="disruptive_missense")
      {
        results_disruptive_missense_genome <- rbind(results_disruptive_missense_genome,results)
      }
      ### synonymous
      if(results[3]=="synonymous")
      {
        results_synonymous_genome <- rbind(results_synonymous_genome,results)
      }
    }

    if(kk%%1000==0)
    {
      print(kk)
    }
  }

  ###### whole-genome results
  # plof
  save(results_plof_genome,file=paste0(output_path,"plof.Rdata"))
  # plof + disruptive missense
  save(results_plof_ds_genome,file=paste0(output_path,"plof_ds.Rdata"))
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
  # missense
  missense_sig <- results_missense_genome[results_missense_genome[,"MetaSTAAR-O"]<alpha,,drop=FALSE]
  write.csv(missense_sig,file=paste0(output_path,"missense_sig.csv"))
  # synonymous
  synonymous_sig <- results_synonymous_genome[results_synonymous_genome[,"MetaSTAAR-O"]<alpha,,drop=FALSE]
  write.csv(synonymous_sig,file=paste0(output_path,"synonymous_sig.csv"))
  # plof_ds
  plof_ds_sig <- results_plof_ds_genome[results_plof_ds_genome[,"MetaSTAAR-O"]<alpha,,drop=FALSE]
  write.csv(plof_ds_sig,file=paste0(output_path,"plof_ds_sig.csv"))
  # disruptive_missense
  disruptive_missense_sig <- results_disruptive_missense_genome[results_disruptive_missense_genome[,"MetaSTAAR-O"]<alpha,,drop=FALSE]
  write.csv(disruptive_missense_sig,file=paste0(output_path,"disruptive_missense_sig.csv"))
  # coding results
  coding_sig <- rbind(plof_sig[,c("Gene name","Chr","Category","#SNV","SKAT-MS(1,25)","Burden-MS(1,1)","ACAT-V-MS(1,25)","MetaSTAAR-O")],
                      missense_sig[,c("Gene name","Chr","Category","#SNV","SKAT-MS(1,25)","Burden-MS(1,1)","ACAT-V-MS(1,25)","MetaSTAAR-O")])
  coding_sig <- rbind(coding_sig,synonymous_sig[,c("Gene name","Chr","Category","#SNV","SKAT-MS(1,25)","Burden-MS(1,1)","ACAT-V-MS(1,25)","MetaSTAAR-O")])
  coding_sig <- rbind(coding_sig,plof_ds_sig[,c("Gene name","Chr","Category","#SNV","SKAT-MS(1,25)","Burden-MS(1,1)","ACAT-V-MS(1,25)","MetaSTAAR-O")])
  coding_sig <- rbind(coding_sig,disruptive_missense_sig[,c("Gene name","Chr","Category","#SNV","SKAT-MS(1,25)","Burden-MS(1,1)","ACAT-V-MS(1,25)","MetaSTAAR-O")])
  write.csv(coding_sig,file=paste0(output_path,"coding_sig.csv"))

  ## manhattan plot
  if(manhattan_plot)
  {
    ### plof
    results_MetaSTAAR <- results_plof_genome[,c(1,2,dim(results_plof_genome)[2])]

    results_m <- c()
    for(i in 1:dim(results_MetaSTAAR)[2])
    {
      results_m <- cbind(results_m,unlist(results_MetaSTAAR[,i]))
    }

    colnames(results_m) <- colnames(results_MetaSTAAR)
    results_m <- data.frame(results_m,stringsAsFactors = FALSE)
    results_m[,2] <- as.numeric(results_m[,2])
    results_m[,3] <- as.numeric(results_m[,3])

    genes_info_manhattan <- dplyr::left_join(genes_info,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))

    genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
    colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "plof"

    ### plof_ds
    results_MetaSTAAR <- results_plof_ds_genome[,c(1,2,dim(results_plof_ds_genome)[2])]

    results_m <- c()
    for(i in 1:dim(results_MetaSTAAR)[2])
    {
      results_m <- cbind(results_m,unlist(results_MetaSTAAR[,i]))
    }

    colnames(results_m) <- colnames(results_MetaSTAAR)
    results_m <- data.frame(results_m,stringsAsFactors = FALSE)
    results_m[,2] <- as.numeric(results_m[,2])
    results_m[,3] <- as.numeric(results_m[,3])

    genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
    genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
    colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "plof_ds"

    ### missense
    results_MetaSTAAR <- results_missense_genome[,c(1,2,dim(results_missense_genome)[2]-6)]

    results_m <- c()
    for(i in 1:dim(results_MetaSTAAR)[2])
    {
      results_m <- cbind(results_m,unlist(results_MetaSTAAR[,i]))
    }

    colnames(results_m) <- colnames(results_MetaSTAAR)
    results_m <- data.frame(results_m,stringsAsFactors = FALSE)
    results_m[,2] <- as.numeric(results_m[,2])
    results_m[,3] <- as.numeric(results_m[,3])

    genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
    genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
    colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "missense"

    ### disruptive_missense
    results_MetaSTAAR <- results_disruptive_missense_genome[,c(1,2,dim(results_disruptive_missense_genome)[2])]

    results_m <- c()
    for(i in 1:dim(results_MetaSTAAR)[2])
    {
      results_m <- cbind(results_m,unlist(results_MetaSTAAR[,i]))
    }

    colnames(results_m) <- colnames(results_MetaSTAAR)
    results_m <- data.frame(results_m,stringsAsFactors = FALSE)
    results_m[,2] <- as.numeric(results_m[,2])
    results_m[,3] <- as.numeric(results_m[,3])

    genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
    genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
    colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "disruptive_missense"

    ### synonymous
    results_MetaSTAAR <- results_synonymous_genome[,c(1,2,dim(results_synonymous_genome)[2])]

    results_m <- c()
    for(i in 1:dim(results_MetaSTAAR)[2])
    {
      results_m <- cbind(results_m,unlist(results_MetaSTAAR[,i]))
    }

    colnames(results_m) <- colnames(results_MetaSTAAR)
    results_m <- data.frame(results_m,stringsAsFactors = FALSE)
    results_m[,2] <- as.numeric(results_m[,2])
    results_m[,3] <- as.numeric(results_m[,3])

    genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
    genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
    colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "synonymous"

    ## ylim
    coding_minp <- min(genes_info_manhattan[,(dim(genes_info_manhattan)[2]-4):dim(genes_info_manhattan)[2]])
    min_y <- ceiling(-log10(coding_minp)) + 1

    pch <- c(0,1,2,3,4)
    figure1 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$plof,sig.level=alpha,pch=pch[1],col = c("blue4", "orange3"),ylim=c(0,min_y),
                              auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous"))))

    figure2 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$plof_ds,sig.level=alpha,pch=pch[2],col = c("blue4", "orange3"),ylim=c(0,min_y),
                              auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous"))))

    figure3 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$missense,sig.level=alpha,pch=pch[3],col = c("blue4", "orange3"),ylim=c(0,min_y),
                              auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous"))))

    figure4 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$disruptive_missense,sig.level=alpha,pch=pch[4],col = c("blue4", "orange3"),ylim=c(0,min_y),
                              auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous"))))

    figure5 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$synonymous,sig.level=alpha,pch=pch[5],col = c("blue4", "orange3"),ylim=c(0,min_y),
                              auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous"))))

    print("Manhattan plot")

    png(paste0(output_path,"gene_centric_coding_manhattan.png"), width = 9, height = 6, units = 'in', res = 600)

    print(figure1)
    print(figure2,newpage = FALSE)
    print(figure3,newpage = FALSE)
    print(figure4,newpage = FALSE)
    print(figure5,newpage = FALSE)

    dev.off()
  }

  ## Q-Q plot
  if(QQ_plot)
  {
    if(!manhattan_plot)
    {
      ### plof
      results_MetaSTAAR <- results_plof_genome[,c(1,2,dim(results_plof_genome)[2])]

      results_m <- c()
      for(i in 1:dim(results_MetaSTAAR)[2])
      {
        results_m <- cbind(results_m,unlist(results_MetaSTAAR[,i]))
      }

      colnames(results_m) <- colnames(results_MetaSTAAR)
      results_m <- data.frame(results_m,stringsAsFactors = FALSE)
      results_m[,2] <- as.numeric(results_m[,2])
      results_m[,3] <- as.numeric(results_m[,3])

      genes_info_manhattan <- dplyr::left_join(genes_info,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))

      genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
      colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "plof"

      ### plof_ds
      results_MetaSTAAR <- results_plof_ds_genome[,c(1,2,dim(results_plof_ds_genome)[2])]

      results_m <- c()
      for(i in 1:dim(results_MetaSTAAR)[2])
      {
        results_m <- cbind(results_m,unlist(results_MetaSTAAR[,i]))
      }

      colnames(results_m) <- colnames(results_MetaSTAAR)
      results_m <- data.frame(results_m,stringsAsFactors = FALSE)
      results_m[,2] <- as.numeric(results_m[,2])
      results_m[,3] <- as.numeric(results_m[,3])

      genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
      genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
      colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "plof_ds"

      ### missense
      results_MetaSTAAR <- results_missense_genome[,c(1,2,dim(results_missense_genome)[2]-6)]

      results_m <- c()
      for(i in 1:dim(results_MetaSTAAR)[2])
      {
        results_m <- cbind(results_m,unlist(results_MetaSTAAR[,i]))
      }

      colnames(results_m) <- colnames(results_MetaSTAAR)
      results_m <- data.frame(results_m,stringsAsFactors = FALSE)
      results_m[,2] <- as.numeric(results_m[,2])
      results_m[,3] <- as.numeric(results_m[,3])

      genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
      genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
      colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "missense"

      ### disruptive_missense
      results_MetaSTAAR <- results_disruptive_missense_genome[,c(1,2,dim(results_disruptive_missense_genome)[2])]

      results_m <- c()
      for(i in 1:dim(results_MetaSTAAR)[2])
      {
        results_m <- cbind(results_m,unlist(results_MetaSTAAR[,i]))
      }

      colnames(results_m) <- colnames(results_MetaSTAAR)
      results_m <- data.frame(results_m,stringsAsFactors = FALSE)
      results_m[,2] <- as.numeric(results_m[,2])
      results_m[,3] <- as.numeric(results_m[,3])

      genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
      genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
      colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "disruptive_missense"

      ### synonymous
      results_MetaSTAAR <- results_synonymous_genome[,c(1,2,dim(results_synonymous_genome)[2])]

      results_m <- c()
      for(i in 1:dim(results_MetaSTAAR)[2])
      {
        results_m <- cbind(results_m,unlist(results_MetaSTAAR[,i]))
      }

      colnames(results_m) <- colnames(results_MetaSTAAR)
      results_m <- data.frame(results_m,stringsAsFactors = FALSE)
      results_m[,2] <- as.numeric(results_m[,2])
      results_m[,3] <- as.numeric(results_m[,3])

      genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
      genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
      colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "synonymous"

      ## ylim
      coding_minp <- min(genes_info_manhattan[,(dim(genes_info_manhattan)[2]-4):dim(genes_info_manhattan)[2]])
      min_y <- ceiling(-log10(coding_minp)) + 1
    }

    print("Q-Q plot")
    cex_point <- 1

    png(paste0(output_path,"gene_centric_coding_qqplot.png"), width = 9, height = 9, units = 'in', res = 600)

    ### plof
    ## remove unconverged p-values
    observed <- sort(genes_info_manhattan$plof[genes_info_manhattan$plof < 1])
    lobs <- -(log10(observed))

    expected <- c(1:length(observed))
    lexp <- -(log10(expected / (length(expected)+1)))

    # par(mar=c(5,6,4,4))
    plot(lexp,lobs,pch=0, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
         xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
         font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)

    abline(0, 1, col="red",lwd=1)

    ### plof_ds
    ## remove unconverged p-values
    observed <- sort(genes_info_manhattan$plof_ds[genes_info_manhattan$plof_ds < 1])

    lobs <- -(log10(observed))

    expected <- c(1:length(observed))
    lexp <- -(log10(expected / (length(expected)+1)))

    par(new=T)
    # par(mar=c(5,6,4,4))
    plot(lexp,lobs,pch=1, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
         xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
         font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)

    abline(0, 1, col="red",lwd=1)

    ### missense
    ## remove unconverged p-values
    observed <- sort(genes_info_manhattan$missense[genes_info_manhattan$missense < 1])

    lobs <- -(log10(observed))

    expected <- c(1:length(observed))
    lexp <- -(log10(expected / (length(expected)+1)))

    par(new=T)
    # par(mar=c(5,6,4,4))
    plot(lexp,lobs,pch=2, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
         xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
         font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)

    abline(0, 1, col="red",lwd=1)

    ### disruptive_missense
    ## remove unconverged p-values
    observed <- sort(genes_info_manhattan$disruptive_missense[genes_info_manhattan$disruptive_missense < 1])


    lobs <- -(log10(observed))

    expected <- c(1:length(observed))
    lexp <- -(log10(expected / (length(expected)+1)))

    par(new=T)
    # par(mar=c(5,6,4,4))
    plot(lexp,lobs,pch=3, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
         xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
         font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)

    abline(0, 1, col="red",lwd=1)

    ### synonymous
    ## remove unconverged p-values
    observed <- sort(genes_info_manhattan$synonymous[genes_info_manhattan$synonymous < 1])

    lobs <- -(log10(observed))

    expected <- c(1:length(observed))
    lexp <- -(log10(expected / (length(expected)+1)))

    par(new=T)
    # par(mar=c(5,6,4,4))
    plot(lexp,lobs,pch=4, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
         xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
         font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)

    abline(0, 1, col="red",lwd=1)

    legend("topleft",legend=c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous"),ncol=1,bty="o",box.lwd=1,pch=0:4,cex=1,text.font=2)

    dev.off()
  }

}

