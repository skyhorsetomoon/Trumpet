.ctest_evluate<-function(result, IP_BAM, Input_BAM, contrast_IP_BAM, contrast_Input_BAM, condition1, condition2)
{
  s<-result[[1]]
  ind<-unique(s$pos)
  len <- length(ind)
  n <- nrow(s)
  se <- seq(1, n, len)
  sa<-s[,-(1:2)]
  sample_name<-.get.sampleid(IP_BAM, Input_BAM,contrast_IP_BAM,contrast_Input_BAM)
  if(length(sample_name)==2)
  {
    IP_groupname<-sample_name[[1]]
    Input_groupname<-sample_name[[2]]
    reference_IP_groupname<-NULL
    reference_Input_groupname<-NULL
  }
  if(length(sample_name)==4)
  {
    IP_groupname<-sample_name[[1]]
    Input_groupname<-sample_name[[2]]
    reference_IP_groupname<-sample_name[[3]]
    reference_Input_groupname<-sample_name[[4]]
  }
  ct_unit_bam<-function(group)
  {
    
    m<-vector()
    group<-as.matrix(group)
    v<-rep(0,length(tr_ve<-.trans_readsvector(group[,1], se, len, n)))
    for(i in seq_len(ncol(group)))
    {
      
      m<-.trans_readsvector(group[,i], se, len, n)
      v<-m+v
      
    }
    v<-as.matrix(v)
    return(v)
  }
  ##
  Ct_ev<-function(bam1,bam_name,bam2,maxfold)
  {
    
    if(ncol(bam1)>1)#bam1=group_IP
    {
      
      per<-vector(mode="numeric",length=0)
      percent<-vector(mode="numeric",length=0)
      for(i in 1:length(bam_name))
      {
        bam<-.trans_readsvector(bam1[,i], se, len, n)
        total_IP<-sum(bam)
        total_Input<-sum(bam2)
        for(i in seq_len(maxfold))
        {
          Ctest<-ctest(bam,bam2,total_IP,total_Input,FOLD = i,minimal_counts_in_fdr = 10)
          cest<-as.data.frame(Ctest)
          ce<-cest[cest$log.fdr<log(0.05),]
          per[i]<-length(ce$log.fdr)/length(cest$log.fdr)
        }
        percent<-cbind(percent,per)
        per<-vector(mode="numeric",length=0)
      }
      colnames(percent)<-bam_name
      percent<-as.data.frame(percent)
      fold<-seq_len(maxfold)
      com<-cbind(fold,percent)
      com<-as.data.frame(com)
      out<-list(com,percent)
      return(out)
    }
    if(ncol(bam1)==1)
    {
      per<-vector(mode="numeric",length=0)
      percent<-vector(mode="numeric",length=0)
      for(i in seq_len(length(bam_name)))
      {
        bam<-.trans_readsvector(bam2[,i], se, len, n)
        total_IP<-sum(bam1)
        total_Input<-sum(bam)
        per<-vector(mode="numeric",length=0)
        for(i in seq_len(maxfold))
        {
          Ctest<-ctest(bam1,bam,total_IP,total_Input,FOLD = i,minimal_counts_in_fdr = 10)
          cest<-as.data.frame(Ctest)
          ce<-cest[cest$log.fdr<log(0.05),]
          per[i]<-length(ce$log.fdr)/length(cest$log.fdr)
        }
        percent<-cbind(percent,per)
        per<-vector(mode="numeric",length=0)
      }
      colnames(percent)<-bam_name
      percent<-as.data.frame(percent)
      fold<-seq_len(maxfold)
      com<-cbind(fold,percent)
      com<-as.data.frame(com)
      out<-list(com,percent)
      return(out)
    }
  }
  ###
  if((length(reference_IP_groupname)==0)&(length(reference_Input_groupname)==0))
  {
    Group_IP<-sa[,(1:length(IP_groupname))]
    Group_IP<-as.matrix(Group_IP)
    Group_Input<-sa[,-(1:length(IP_groupname))]
    Group_Input<-as.matrix(Group_Input)
    Unit_Input<-ct_unit_bam(Group_Input)
    Unit_IP<-ct_unit_bam(Group_IP)
    ##IP
    out1<-suppressWarnings(Ct_ev(Group_IP,IP_groupname,Unit_Input,maxfold = 10))
    com1<-out1[[1]]
    percent1<-out1[2]
    com1<-melt(data=com1,id="fold",measure.vars=c(2:(ncol(Group_IP)+1)))
    ##Input
    out2<-suppressWarnings(Ct_ev(Unit_IP,Input_groupname,Group_Input,maxfold = 10))
    com2<-out2[[1]]
    percent2<-out2[[2]]
    com2<-melt(data=com2,id="fold",measure.vars=c(2:(ncol(Group_Input)+1)))
    com1<-as.data.frame(com1)
    fold<-com1$fold
    value<-com1$value
    variable<-com1$variable
    c_p1<-
      ggplot(data=com1,aes(x=fold,y=value,colour=variable))+
      geom_point(aes(shape=variable))+
      geom_line()+
      ylab("Percentgae of bins")+
      xlab("foldchange")+
      labs(title=paste("IP's foldchange under unified Input within",condition1,"condition"))
    com2<-as.data.frame(com2)
    fold<-com2$fold
    value<-com2$value
    variable<-com2$variable
    c_p2<-
      ggplot(data=com2,aes(x=fold,y=value,colour=variable))+
      geom_point(aes(shape=variable))+
      geom_line()+
      ylab("Percentgae of bins")+
      xlab("foldchange")+
      labs(title=paste("unified IP's foldchange under different Input within ",condition1,"condition"))
    .multiplot(c_p1,c_p2,cols=1)
  }
  else if((length(reference_IP_groupname)!=0)&(length(reference_Input_groupname)!=0))
  {
    group_IP<-sa[,(1:length(IP_groupname))]
    group_IP<-as.matrix(group_IP)
    Group_Input<-sa[,-(1:length(IP_groupname))]
    group_Input<-Group_Input[,-((length(Input_groupname)+1):ncol(Group_Input))]
    group_Input<-as.matrix(group_Input)
    ref_group<-Group_Input[,-(1:length(Input_groupname))]
    ref_group_IP<-ref_group[,1:length(reference_Input_groupname)]
    ref_group_IP<-as.matrix(ref_group_IP)
    ref_group_Input<-ref_group[,-(1:length(reference_Input_groupname))]
    ref_group_Input<-as.matrix(ref_group_Input)
    Unit_Input<-ct_unit_bam(group_Input)
    Unit_IP<-ct_unit_bam(group_IP)
    ref_unit_Input<-ct_unit_bam(ref_group_Input)
    ref_unit_IP<-ct_unit_bam(ref_group_IP)
    ##IP
    out1<-suppressWarnings(Ct_ev(group_IP,IP_groupname,Unit_Input,maxfold = 10))
    com1<-out1[[1]]
    percent1<-out1[[2]]
    com<-melt(data=com1,id="fold",measure.vars=c(2:(ncol(group_IP)+1)))
    ##Input
    out2<-suppressWarnings(Ct_ev(Unit_IP,Input_groupname,group_Input,maxfold = 10))
    com2<-out2[[1]]
    percent2<-out2[[2]]
    com2<-melt(data=com2,id="fold",measure.vars=c(2:(ncol(group_Input)+1)))
    ##refer_IP
    refer_out1<-suppressWarnings(Ct_ev(ref_group_IP,reference_IP_groupname,ref_unit_Input,maxfold = 10))
    refer_com1<-refer_out1[[1]]
    refer_percent1<-refer_out1[[2]]
    refer_com1<-melt(data=refer_com1,id="fold",measure.vars=c(2:(ncol(ref_group_IP)+1)))
    ##refer_Input
    refer_out2<-suppressWarnings(Ct_ev(ref_unit_IP,reference_Input_groupname,ref_group_Input,maxfold = 10))
    refer_com2<-refer_out2[[1]]
    refer_percent2<-refer_out2[[2]]
    refer_com2<-melt(data=refer_com2,id="fold",measure.vars=c(2:(ncol(ref_group_Input)+1)))
    com<-as.data.frame(com)
    fold<-com$fold
    value<-com$value
    variable<-com$variable
    c_p1<-
      ggplot(data=com,aes(x=fold,y=value,colour=variable))+
      geom_point(aes(shape=variable))+
      geom_line()+
      ylab("Percentgae of bins")+
      xlab("foldchange")+
      labs(title=paste("IP's foldchange under unified Input within",condition1,"condition"))
    com2<-as.data.frame(com2)
    fold<-com2$fold
    value<-com2$value
    variable<-com2$variable
    c_p2<-
      ggplot(data=com2,aes(x=fold,y=value,colour=variable))+
      geom_point(aes(shape=variable))+
      geom_line()+
      ylab("Percentgae of bins")+
      xlab("foldchange")+
      labs(title=paste("unified IP's foldchange under different Input within ",condition1,"condition"))
    refer_com1<-as.data.frame(refer_com1)
    fold<-refer_com1$fold
    value<-refer_com1$value
    variable<-refer_com1$variable
    refer_p1<-
      ggplot(data=refer_com1,aes(x=fold,y=value,colour=variable))+
      geom_point(aes(shape=variable))+
      geom_line()+
      ylab("Percentgae of bins")+
      xlab("foldchange")+
      labs(title=paste("refer_IP's foldchange under unified Input within" ,condition2,"condition"))
    refer_com2<-as.data.frame(refer_com2)
    fold<-refer_com2$fold
    value<-refer_com2$value
    variable<-refer_com2$variable
    refer_p2<-ggplot(refer_com2,aes(x=fold,y=value,colour=variable))+
      geom_point(aes(shape=variable))+
      geom_line()+
      ylab("Percentgae of bins")+xlab("foldchange")+
      labs(title=paste("unified refer_IP's foldchange under different Input within",condition2,"condition"))
    .multiplot(c_p1,c_p2,refer_p1,refer_p2,cols=2)
  }
}
