.get.sampleid <- function(IP_BAM, Input_BAM, contrast_IP_BAM, contrast_Input_BAM) {
  
  if ((length(seq_along(contrast_IP_BAM)) == 0) & (length(seq_along(contrast_Input_BAM)) == 
                                                   0)) {
    IP_name <- paste0("IP", seq_along(IP_BAM))
    Input_name <- paste0("Input", seq_along(Input_BAM))
    referIP_name <- NULL
    referInput_name <- NULL
  }
  
  if ((length(seq_along(contrast_IP_BAM)) > 0) & (length(seq_along(contrast_Input_BAM)) > 
                                                  0)) {
    
    IP_name <- paste0("IP", seq_along(IP_BAM))
    Input_name <- paste0("Input", seq_along(Input_BAM))
    referIP_name <- paste0("refer_IP", seq_along(contrast_IP_BAM))
    referInput_name <- paste0("refer_Input", seq_along(contrast_Input_BAM))
  }
  sample_name <- list()
  sample_name$IP_name <- IP_name
  sample_name$Input_name <- Input_name
  sample_name$referIP_name <- referIP_name
  sample_name$referInput_name <- referInput_name
  return(sample_name)
}


