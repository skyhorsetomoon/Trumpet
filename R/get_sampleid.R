
.get.sampleid <- function(IP_BAM,Input_BAM,contrast_IP_BAM,contrast_Input_BAM) {
  # number of samples
  no_ip=length(IP_BAM)
  no_input=length(Input_BAM)
  no_treated_ip=length(contrast_IP_BAM)
  no_treated_input=length(contrast_Input_BAM)

  if ((no_treated_input == 0)&(no_treated_ip ==0)) {
    for(i in 1:no_ip)
    {
      IP_name<-paste0("IP",c(1:no_ip))
    }
    for(i in 1: no_input)
    {
      Input_name<-paste0("Input",c(1:no_input))
    }
    referIP_name<-NULL
    referInput_name<-NULL
    }
  
  if ((no_treated_input > 0)&(no_treated_ip >0) ){
    for(i in 1:no_ip)
    {
      IP_name<-paste0("IP",c(1:no_ip))
    }
    for(i in 1:no_input)
    {
      Input_name<-paste0("Input",c(1:no_input))
    }
    for(i in 1:no_treated_ip)
    {
      referIP_name<-paste0("refer_IP",c(1:no_treated_ip))
    }
    for(i in 1: no_treated_input)
    {
      referInput_name<-paste0("refer_Input",c(1:no_treated_input))
    }
    } 
  sample_name<-list()
  sample_name$IP_name<-IP_name
  sample_name$Input_name<-Input_name
  sample_name$referIP_name<-referIP_name
  sample_name$referInput_name<-referInput_name
  return(sample_name)
}