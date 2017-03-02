two_cond_ms<-function(result,condition1,condition2)
{
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
  if(is.null(reference_IP_groupname)&is.null(reference_Input_groupname))
  {
    print("Must provide two condition bam files using this function")
  }
  s<-result[[1]]
  ind<-unique(s$pos)
  len <- length(ind)
  n <- nrow(s)
  se <- seq(1, n, len)
  sa<-s[,-(1:2)]
  con_ms_f<-function(group1,group2,cond_name1,cond_name2)
  {
    f<-function(bam)
    {
      w<-bam[se]
      for (i in 2:length(ind))
      {
        se <- seq(i, n, len)
        w<- cbind(w, bam[se])
      }
      w<-as.matrix(w)
      w<-sapply(t(w),unlist)
      return(w)
    }
    com_bam<-function(group)
    {
      v<-vector()
      for(i in 1:ncol(group))
      {
        m<-vector()
        m<-f(group[,i])
        v<-rbind(v,m)
      }
      return(v)
    }
    p1<-com_bam(group_IP)
    p2<-com_bam(ref_group_IP)
    size1<-rowSums(p1)
    size_factor1 <- size1/exp(mean(log(size1)))
    p1<-apply(p1,2,function(x,a)x/a,a=size_factor1)
    size2 <- rowSums(p2)
    size_factor2<- size2/exp(mean(log(size2)))
    p2<-apply(p2,2,function(x,a)x/a,a=size_factor2)
    Mean<-apply(p1,2,mean)
    Sd<-apply(p1,2,sd)
    com1<-cbind(Mean,Sd)
    com1<-as.data.frame(com1)
    z<-which(com1$Mean==0)
    com1<-com1[-z,]
    Mean<-log10(com1$Mean)
    Sd<-log10(com1$Sd)
    com1<-cbind(Mean,Sd)
    com1<-as.data.frame(com1)
    Mean<-apply(p2,2,mean)
    Sd<-apply(p2,2,sd)
    com2<-cbind(Mean,Sd)
    com2<-as.data.frame(com2)
    z<-which(com2$Mean==0)
    com2<-com2[-z,]
    Mean<-log10(com2$Mean)
    Sd<-log10(com2$Sd)
    com2<-cbind(Mean,Sd)
    com2<-as.data.frame(com2)
    com<-rbind(com1,com2)
    com<-as.data.frame(com)
    ID<-rep(c(cond_name1,cond_name2),c(length(com1$Mean),length(com2$Mean)))
    com<-cbind(com,ID)
    com<-as.data.frame(com)
    return(com)
  }
  if((length(reference_IP_groupname)!=0)&(length(reference_Input_groupname)!=0)&((length(reference_IP_groupname)+length(reference_Input_groupname))<=2|(length(IP_groupname)+length(Input_groupname))<=2))
  {
    print("The number of samples in each condition should be more than three when using this function")
  }
  if((length(reference_IP_groupname)!=0)&(length(reference_Input_groupname)!=0)&((length(reference_IP_groupname)+length(reference_Input_groupname))>2))
  {
    group_IP<-sa[,(1:length(IP_groupname))]
    group_IP<-as.matrix(group_IP)
    Group_Input<-sa[,-(1:length(IP_groupname))]
    group_Input<-Group_Input[,-((length(Input_groupname)+1):ncol(Group_Input))]
    group_Input<-as.matrix(group_Input)
    ref_group<-Group_Input[,-(1:length(Input_groupname))]
    ref_group_IP<-ref_group[,1:length(reference_IP_groupname)]
    ref_group_IP<-as.matrix(ref_group_IP)
    ref_group_Input<-ref_group[,-(1:length(reference_IP_groupname))]
    ref_group_Input<-as.matrix(ref_group_Input)
    m_com1<-con_ms_f(group_IP,ref_group_IP,paste("IP group in",condition1,"condition"),paste("IP group under",condition2,"condition"))
    m_com2<-con_ms_f(ref_group_IP,ref_group_Input,paste("Input group in",condition1,"condition"),paste("Input group under",condition2,"condition"))
    lp1<-
      ggplot(m_com1,aes(Mean,Sd,colour=ID))+
      geom_smooth(aes(group =ID),span = 0.5)+
      geom_point(alpha = I(1 /200),size=0.002)+
      theme(title= element_text(size=14,color="black"))+
      labs(title=paste(" IP group's Mean-variance relationship within two condition"))
    lp2<-
      ggplot(m_com2,aes(Mean,Sd,colour=ID))+
      geom_smooth(aes(group =ID),span = 0.5)+
      geom_point(alpha = I(1 /200),size=0.002)+
      theme(title= element_text(size=14,color="black"))+
      labs(title=paste("Input group's Mean-variance relationship within two condition"))
    .multiplot(lp1,lp2,cols = 2)
  }
}

