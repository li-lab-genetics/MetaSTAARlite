#' Generates summary table and/or visualization for the meta-analysis of noncoding functional categories that was conducted.
#'
#' This function generates a summary table, Manhattan plot, and QQ plot for the meta-analysis of noncoding functional categories that was conducted
#' based on the parameters provided by the user.
#'
#' @param gene_centric_noncoding_jobs_num an integer which specifies the number of jobs done in the gene-centric noncoding meta-analysis.
#' @param input_path a character which specifies the file path to the gene-centric noncoding meta-analysis results files.
#' @param output_path a character which specifies the file path to the desired location of the produced summary table and visualizations for gene-centric noncoding meta-analysis.
#' @param gene_centric_results_name a character which specifies the name (excluding the jobs number) of the gene-centric noncoding meta-analysis results files.
#' @param ncRNA_jobs_num an integer which specifies the number of jobs done in the ncRNA meta-analysis.
#' @param ncRNA_input_path a character which specifies the file path to the ncRNA meta-analysis results files.
#' @param ncRNA_output_path a character which specifies the file path to the desired location of the produced summary table and visualizations for the ncRNA meta-analysis.
#' @param ncRNA_results_name a character which specifies the name (excluding the jobs number) of the ncRNA meta-analysis results files.
#' @param alpha a numeric value which specifies the desired significance threshold for the gene-centric noncoding meta-analysis. Default is 2.5E-06.
#' @param alpha_ncRNA a numeric value which specifies the desired significance threshold for the ncRNA meta-analysis. Default is 2.5E-06.
#' @param manhattan_plot a logical value which determines if a Manhattan plot is generated. Default is FALSE.
#' @param QQ_plot a logical value which determines if a QQ plot is generated. Default is FALSE.

Gene_Centric_Noncoding_Results_Summary_meta <- function(gene_centric_noncoding_jobs_num,input_path,output_path,gene_centric_results_name,
                                                        ncRNA_jobs_num,ncRNA_input_path,ncRNA_output_path,ncRNA_results_name,
                                                        alpha=2.5E-06,alpha_ncRNA=2.5E-06,
                                                        ncRNA_pos=NULL,manhattan_plot=FALSE,QQ_plot=FALSE){

  #######################################################
  #     summarize unconditional analysis results
  #######################################################

  results_noncoding_genome <- c()

  for(kk in 1:gene_centric_noncoding_jobs_num)
  {
    print(kk)
    results_noncoding <- get(load(paste0(input_path,gene_centric_results_name,"_",kk,".Rdata")))

    results_noncoding_genome <- c(results_noncoding_genome,results_noncoding)
  }

  results_UTR_genome <- c()
  results_upstream_genome <- c()
  results_downstream_genome <- c()
  results_promoter_CAGE_genome <- c()
  results_promoter_DHS_genome <- c()
  results_enhancer_CAGE_genome <- c()
  results_enhancer_DHS_genome <- c()

  for(kk in 1:length(results_noncoding_genome))
  {
    results <- results_noncoding_genome[[kk]]
    if(is.null(results)==FALSE)
    {
      ### UTR
      if(results[3]=="UTR")
      {
        results_UTR_genome <- rbind(results_UTR_genome,results)
      }
      ### upstream
      if(results[3]=="upstream")
      {
        results_upstream_genome <- rbind(results_upstream_genome,results)
      }
      ### downstream
      if(results[3]=="downstream")
      {
        results_downstream_genome <- rbind(results_downstream_genome,results)
      }
      ### promoter_CAGE
      if(results[3]=="promoter_CAGE")
      {
        results_promoter_CAGE_genome <- rbind(results_promoter_CAGE_genome,results)
      }
      ### promoter_DHS
      if(results[3]=="promoter_DHS")
      {
        results_promoter_DHS_genome <- rbind(results_promoter_DHS_genome,results)
      }
      ### enhancer_CAGE
      if(results[3]=="enhancer_CAGE")
      {
        results_enhancer_CAGE_genome <- rbind(results_enhancer_CAGE_genome,results)
      }
      ### enhancer_DHS
      if(results[3]=="enhancer_DHS")
      {
        results_enhancer_DHS_genome <- rbind(results_enhancer_DHS_genome,results)
      }
    }
    if(kk%%1000==0)
    {
      print(kk)
    }
  }

  ###### whole-genome results
  # UTR
  save(results_UTR_genome,file=paste0(output_path,"UTR.Rdata"))
  # upstream
  save(results_upstream_genome,file=paste0(output_path,"upstream.Rdata"))
  # downstream
  save(results_downstream_genome,file=paste0(output_path,"downstream.Rdata"))
  # promoter CAGE
  save(results_promoter_CAGE_genome,file=paste0(output_path,"promoter_CAGE.Rdata"))
  # promoter DHS
  save(results_promoter_DHS_genome,file=paste0(output_path,"promoter_DHS.Rdata"))
  # enhancer CAGE
  save(results_enhancer_CAGE_genome,file=paste0(output_path,"enhancer_CAGE.Rdata"))
  # enahncer DHS
  save(results_enhancer_DHS_genome,file=paste0(output_path,"enhancer_DHS.Rdata"))

  ###### ncRNA
  results_ncRNA_genome <- c()

  for(kk in 1:ncRNA_jobs_num)
  {
    print(kk)
    results_ncRNA <- get(load(paste0(ncRNA_input_path,ncRNA_results_name,"_",kk,".Rdata")))
    results_ncRNA_genome <- rbind(results_ncRNA_genome,results_ncRNA)
  }

  ###### whole-genome results
  save(results_ncRNA_genome,file=paste0(ncRNA_output_path,"results_ncRNA_genome.Rdata"))

  ###### significant results
  ### UTR
  UTR_sig <- results_UTR_genome[results_UTR_genome[,"MetaSTAAR-O"]<alpha,,drop=FALSE]
  write.csv(UTR_sig,file=paste0(output_path,"UTR_sig.csv"))

  noncoding_sig <- c()
  noncoding_sig <- rbind(noncoding_sig,UTR_sig)

  ### upstream
  upstream_sig <- results_upstream_genome[results_upstream_genome[,"MetaSTAAR-O"]<alpha,,drop=FALSE]
  write.csv(upstream_sig,file=paste0(output_path,"upstream_sig.csv"))

  noncoding_sig <- rbind(noncoding_sig,upstream_sig)

  ### downstream
  downstream_sig <- results_downstream_genome[results_downstream_genome[,"MetaSTAAR-O"]<alpha,,drop=FALSE]
  write.csv(downstream_sig,file=paste0(output_path,"downstream_sig.csv"))

  noncoding_sig <- rbind(noncoding_sig,downstream_sig)

  ### promoter_CAGE
  promoter_CAGE_sig <- results_promoter_CAGE_genome[results_promoter_CAGE_genome[,"MetaSTAAR-O"]<alpha,,drop=FALSE]
  write.csv(promoter_CAGE_sig,file=paste0(output_path,"promoter_CAGE_sig.csv"))

  noncoding_sig <- rbind(noncoding_sig,promoter_CAGE_sig)

  ### promoter_DHS
  promoter_DHS_sig <- results_promoter_DHS_genome[results_promoter_DHS_genome[,"MetaSTAAR-O"]<alpha,,drop=FALSE]
  write.csv(promoter_DHS_sig,file=paste0(output_path,"promoter_DHS_sig.csv"))

  noncoding_sig <- rbind(noncoding_sig,promoter_DHS_sig)

  ### enhancer_CAGE
  enhancer_CAGE_sig <- results_enhancer_CAGE_genome[results_enhancer_CAGE_genome[,"MetaSTAAR-O"]<alpha,,drop=FALSE]
  write.csv(enhancer_CAGE_sig,file=paste0(output_path,"enhancer_CAGE_sig.csv"))

  noncoding_sig <- rbind(noncoding_sig,enhancer_CAGE_sig)

  ### enhancer_DHS
  enhancer_DHS_sig <- results_enhancer_DHS_genome[results_enhancer_DHS_genome[,"MetaSTAAR-O"]<alpha,,drop=FALSE]
  write.csv(enhancer_DHS_sig,file=paste0(output_path,"enhancer_DHS_sig.csv"))

  noncoding_sig <- rbind(noncoding_sig,enhancer_DHS_sig)

  ### ncRNA
  ncRNA_sig <- results_ncRNA_genome[results_ncRNA_genome[,"MetaSTAAR-O"]<alpha_ncRNA,,drop=FALSE]
  write.csv(ncRNA_sig,file=paste0(ncRNA_output_path,"ncRNA_sig.csv"))

  noncoding_sig <- rbind(noncoding_sig,ncRNA_sig)

  write.csv(noncoding_sig,file=paste0(output_path,"noncoding_sig.csv"))

  ## manhattan plot for protein coding genes
  if(manhattan_plot)
  {
    ############## noncoding
    ### UTR
    results_MetaSTAAR <- results_UTR_genome[,c(1,2,dim(results_UTR_genome)[2])]

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
    colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "UTR"

    ### upstream
    results_MetaSTAAR <- results_upstream_genome[,c(1,2,dim(results_upstream_genome)[2])]

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
    colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "upstream"

    ### downstream
    results_MetaSTAAR <- results_downstream_genome[,c(1,2,dim(results_downstream_genome)[2])]

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
    colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "downstream"

    ### promoter_CAGE
    results_MetaSTAAR <- results_promoter_CAGE_genome[,c(1,2,dim(results_promoter_CAGE_genome)[2])]

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
    colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "promoter_CAGE"

    ### promoter_DHS
    results_MetaSTAAR <- results_promoter_DHS_genome[,c(1,2,dim(results_promoter_DHS_genome)[2])]

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
    colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "promoter_DHS"

    ### enhancer_CAGE
    results_MetaSTAAR <- results_enhancer_CAGE_genome[,c(1,2,dim(results_enhancer_CAGE_genome)[2])]

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
    colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "enhancer_CAGE"

    ### enhancer_DHS
    results_MetaSTAAR <- results_enhancer_DHS_genome[,c(1,2,dim(results_enhancer_DHS_genome)[2])]

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
    colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "enhancer_DHS"

    ## ylim
    noncoding_minp <- min(genes_info_manhattan[,(dim(genes_info_manhattan)[2]-6):dim(genes_info_manhattan)[2]])
    min_y <- ceiling(-log10(noncoding_minp)) + 1

    pch <- c(0,1,2,3,4,5,6)

    figure1 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$UTR,sig.level=alpha,pch=pch[1],col = c("blue4", "orange3"),ylim=c(0,min_y),
                              auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("UTR","upstream","downstream","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"))))

    figure2 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$upstream,sig.level=alpha,pch=pch[2],col = c("blue4", "orange3"),ylim=c(0,min_y),
                              auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("UTR","upstream","downstream","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"))))

    figure3 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$downstream,sig.level=alpha,pch=pch[3],col = c("blue4", "orange3"),ylim=c(0,min_y),
                              auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("UTR","upstream","downstream","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"))))

    figure4 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$promoter_CAGE,sig.level=alpha,pch=pch[4],col = c("blue4", "orange3"),ylim=c(0,min_y),
                              auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("UTR","upstream","downstream","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"))))

    figure5 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$promoter_DHS,sig.level=alpha,pch=pch[5],col = c("blue4", "orange3"),ylim=c(0,min_y),
                              auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("UTR","upstream","downstream","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"))))

    figure6 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$enhancer_CAGE,sig.level=alpha,pch=pch[6],col = c("blue4", "orange3"),ylim=c(0,min_y),
                              auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("UTR","upstream","downstream","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"))))

    figure7 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$enhancer_DHS,sig.level=alpha,pch=pch[7],col = c("blue4", "orange3"),ylim=c(0,min_y),
                              auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("UTR","upstream","downstream","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"))))

    print("Manhattan plot")

    png(paste0(output_path,"gene_centric_noncoding_manhattan.png"), width = 9, height = 6, units = 'in', res = 600)

    print(figure1)
    print(figure2,newpage = FALSE)
    print(figure3,newpage = FALSE)
    print(figure4,newpage = FALSE)
    print(figure5,newpage = FALSE)
    print(figure6,newpage = FALSE)
    print(figure7,newpage = FALSE)

    dev.off()
  }

  ## manhattan plot for ncRNA genes
  if(manhattan_plot)
  {
    if(!is.null(ncRNA_pos))
    {
      results_ncRNA_genome_temp <- data.frame(results_ncRNA_genome[,c(1,2,dim(results_ncRNA_genome)[2])],stringsAsFactors = FALSE)
      results_ncRNA_genome_temp[,2] <- as.numeric(results_ncRNA_genome_temp[,2])
      results_ncRNA_genome_temp[,1] <- as.character(results_ncRNA_genome_temp[,1])
      results_ncRNA_genome_temp[,3] <- as.numeric(results_ncRNA_genome_temp[,3])

      ncRNA_gene_pos_results <- dplyr::left_join(ncRNA_pos,results_ncRNA_genome_temp,by=c("chr"="Chr","ncRNA"="Gene.name"))
      ncRNA_gene_pos_results[is.na(ncRNA_gene_pos_results[,5]),5] <- 1

      print("ncRNA Manhattan plot")

      png(paste0(ncRNA_output_path,"gene_centric_ncRNA_manhattan.png"), width = 9, height = 6, units = 'in', res = 600)

      print(manhattan_plot(as.numeric(ncRNA_gene_pos_results[,1]), (as.numeric(ncRNA_gene_pos_results[,3])+as.numeric(ncRNA_gene_pos_results[,4]))/2, as.numeric(ncRNA_gene_pos_results[,5]), col = c("blue4", "orange3"),sig.level=alpha_ncRNA))

      dev.off()
    }
  }

  ## Q-Q plot for protein coding genes
  if(QQ_plot)
  {
    if(!manhattan_plot)
    {
      ############## noncoding
      ### UTR
      results_MetaSTAAR <- results_UTR_genome[,c(1,2,dim(results_UTR_genome)[2])]

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
      colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "UTR"

      ### upstream
      results_MetaSTAAR <- results_upstream_genome[,c(1,2,dim(results_upstream_genome)[2])]

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
      colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "upstream"

      ### downstream
      results_MetaSTAAR <- results_downstream_genome[,c(1,2,dim(results_downstream_genome)[2])]

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
      colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "downstream"

      ### promoter_CAGE
      results_MetaSTAAR <- results_promoter_CAGE_genome[,c(1,2,dim(results_promoter_CAGE_genome)[2])]

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
      colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "promoter_CAGE"

      ### promoter_DHS
      results_MetaSTAAR <- results_promoter_DHS_genome[,c(1,2,dim(results_promoter_DHS_genome)[2])]

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
      colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "promoter_DHS"

      ### enhancer_CAGE
      results_MetaSTAAR <- results_enhancer_CAGE_genome[,c(1,2,dim(results_enhancer_CAGE_genome)[2])]

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
      colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "enhancer_CAGE"

      ### enhancer_DHS
      results_MetaSTAAR <- results_enhancer_DHS_genome[,c(1,2,dim(results_enhancer_DHS_genome)[2])]

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
      colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "enhancer_DHS"

      ## ylim
      noncoding_minp <- min(genes_info_manhattan[,(dim(genes_info_manhattan)[2]-6):dim(genes_info_manhattan)[2]])
      min_y <- ceiling(-log10(noncoding_minp)) + 1
    }

    print("Q-Q plot")
    cex_point <- 1

    png(paste0(output_path,"gene_centric_noncoding_qqplot.png"), width = 9, height = 9, units = 'in', res = 600)

    ### UTR
    ## remove unconverged p-values
    observed <- sort(genes_info_manhattan$UTR[genes_info_manhattan$UTR < 1])

    lobs <- -(log10(observed))

    expected <- c(1:length(observed))
    lexp <- -(log10(expected / (length(expected)+1)))

    # par(mar=c(5,6,4,4))
    plot(lexp,lobs,pch=0, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
         xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
         font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)

    abline(0, 1, col="red",lwd=1)

    ### upstream
    ## remove unconverged p-values
    observed <- sort(genes_info_manhattan$upstream[genes_info_manhattan$upstream < 1])

    lobs <- -(log10(observed))

    expected <- c(1:length(observed))
    lexp <- -(log10(expected / (length(expected)+1)))

    par(new=T)
    # par(mar=c(5,6,4,4))
    plot(lexp,lobs,pch=1, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
         xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
         font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)

    abline(0, 1, col="red",lwd=1)

    ### downstream
    ## remove unconverged p-values
    observed <- sort(genes_info_manhattan$downstream[genes_info_manhattan$downstream < 1])

    lobs <- -(log10(observed))

    expected <- c(1:length(observed))
    lexp <- -(log10(expected / (length(expected)+1)))

    par(new=T)
    # par(mar=c(5,6,4,4))
    plot(lexp,lobs,pch=2, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
         xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
         font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)

    abline(0, 1, col="red",lwd=1)

    ### promoter_CAGE
    ## remove unconverged p-values
    observed <- sort(genes_info_manhattan$promoter_CAGE[genes_info_manhattan$promoter_CAGE < 1])

    lobs <- -(log10(observed))

    expected <- c(1:length(observed))
    lexp <- -(log10(expected / (length(expected)+1)))

    par(new=T)
    # par(mar=c(5,6,4,4))
    plot(lexp,lobs,pch=3, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
         xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
         font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)

    abline(0, 1, col="red",lwd=1)

    ### promoter_DHS
    ## remove unconverged p-values
    observed <- sort(genes_info_manhattan$promoter_DHS[genes_info_manhattan$promoter_DHS < 1])

    lobs <- -(log10(observed))

    expected <- c(1:length(observed))
    lexp <- -(log10(expected / (length(expected)+1)))

    par(new=T)
    # par(mar=c(5,6,4,4))
    plot(lexp,lobs,pch=4, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
         xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
         font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)

    abline(0, 1, col="red",lwd=1)

    ### enhancer_CAGE
    ## remove unconverged p-values
    observed <- sort(genes_info_manhattan$enhancer_CAGE[genes_info_manhattan$enhancer_CAGE < 1])

    lobs <- -(log10(observed))

    expected <- c(1:length(observed))
    lexp <- -(log10(expected / (length(expected)+1)))

    par(new=T)
    # par(mar=c(5,6,4,4))
    plot(lexp,lobs,pch=5, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
         xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
         font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)

    abline(0, 1, col="red",lwd=1)

    ### enhancer_DHS
    ## remove unconverged p-values
    observed <- sort(genes_info_manhattan$enhancer_DHS[genes_info_manhattan$enhancer_DHS < 1])

    lobs <- -(log10(observed))

    expected <- c(1:length(observed))
    lexp <- -(log10(expected / (length(expected)+1)))

    par(new=T)
    # par(mar=c(5,6,4,4))
    plot(lexp,lobs,pch=6, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
         xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
         font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)

    abline(0, 1, col="red",lwd=1)

    legend("topleft",legend=c("UTR","upstream","downstream","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"),ncol=1,bty="o",box.lwd=1,pch=0:6,cex=1,text.font=2)

    dev.off()
  }

  ## Q-Q plot for ncRNA genes
  if(QQ_plot)
  {
    if(!is.null(ncRNA_pos))
    {
      if(!manhattan_plot)
      {
        results_ncRNA_genome_temp <- data.frame(results_ncRNA_genome[,c(1,2,dim(results_ncRNA_genome)[2])],stringsAsFactors = FALSE)
        results_ncRNA_genome_temp[,2] <- as.numeric(results_ncRNA_genome_temp[,2])
        results_ncRNA_genome_temp[,1] <- as.character(results_ncRNA_genome_temp[,1])
        results_ncRNA_genome_temp[,3] <- as.numeric(results_ncRNA_genome_temp[,3])

        ncRNA_gene_pos_results <- dplyr::left_join(ncRNA_pos,results_ncRNA_genome_temp,by=c("chr"="Chr","ncRNA"="Gene.name"))
        ncRNA_gene_pos_results[is.na(ncRNA_gene_pos_results[,5]),5] <- 1
      }

      observed <- sort(ncRNA_gene_pos_results[ncRNA_gene_pos_results[,5] < 1,5])
      lobs <- -(log10(observed))

      expected <- c(1:length(observed))
      lexp <- -(log10(expected / (length(expected)+1)))

      ncRNA_minp <- min(ncRNA_gene_pos_results[,5])
      min_ncRNA_y <- ceiling(-log10(ncRNA_minp)) + 1

      print("ncRNA Q-Q plot")
      png(paste0(ncRNA_output_path,"gene_centric_ncRNA_qqplot.png"), width = 9, height = 9, units = 'in', res = 600)

      plot(lexp,lobs,pch=20, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_ncRNA_y),
           xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
           font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)

      abline(0, 1, col="red",lwd=1)

      dev.off()
    }
  }

}

