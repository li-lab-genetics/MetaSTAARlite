custom_MetaSTAARlite_cond <- function(chr,mask_name,
                                      sample.sizes,custom_sumstat_mask_list,custom_cov_mask_list,custom_cov_cond_mask_list,
                                      cov_maf_cutoff,rare_maf_cutoff=0.01,rv_num_cutoff=2,effect.cond = c("homogeneous","heterogeneous"),
                                      check_qc_label=FALSE,variant_type=c("SNV","Indel","variant"),
                                      Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,silent=FALSE){
  
  ## evaluate choices
  effect.cond <- match.arg(effect.cond)
  variant_type <- match.arg(variant_type)
  
  results <- c()
  
  mask_test_merge_cond <- MetaSTAARlite_merge_cond(chr,sample.sizes,custom_sumstat_mask_list,custom_cov_mask_list,custom_cov_cond_mask_list,
                                                   cov_maf_cutoff=cov_maf_cutoff,rare_maf_cutoff=rare_maf_cutoff,effect.cond=effect.cond,
                                                   check_qc_label=check_qc_label,variant_type=variant_type,
                                                   Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  if(length(mask_test_merge_cond$U_cond)>=2)
  {
    ## Annotation
    annotation_phred <- mask_test_merge_cond$annotation_phred
    pvalues <- 0
    try(pvalues <- MetaSTAAR_cond(mask_test_merge_cond,annotation_phred,rv_num_cutoff),silent=silent)
    
    if(inherits(pvalues,"list"))
    {
      results_temp <- rep(NA,3)
      results_temp[2] <- chr
      results_temp[1] <- paste0(as.character(mask_name),"_cond")
      results_temp[3] <- pvalues$num_variant
      
      
      results_temp <- c(results_temp,pvalues$results_MetaSTAAR_S_1_25,pvalues$results_MetaSTAAR_S_1_1,
                        pvalues$results_MetaSTAAR_B_1_25,pvalues$results_MetaSTAAR_B_1_1,pvalues$results_MetaSTAAR_A_1_25,
                        pvalues$results_MetaSTAAR_A_1_1,pvalues$results_ACAT_O_MS,pvalues$results_MetaSTAAR_O)
      
      results <- rbind(results,results_temp)
    }
  }
  
  if(!is.null(results))
  {
    colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
    colnames(results)[1:3] <- c("Mask name","Chr","#SNV")
    colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O-MS","MetaSTAAR-O")
  }
  
  return(results)
}
