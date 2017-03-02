ms_relation<-function(result,condition1,condition2)
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
  ms<-function(group1,group2,group1_name,group2_name,condition)
  {
    ms_f<-function(bam)
    {
      s1 <-bam[se]
      for (i in 2:length(ind))
      {
        se <- seq(i, n, len)
        s1 <- cbind(s1,bam[se])
        
      }
      row.mean<-apply(s1,1,mean)
      Mean<-row.mean
      Sd<-apply(s1,1,sd)
      IP<-cbind(Mean,Sd)
      IP<-as.data.frame(IP)
      z<-which(IP$Mean==0)
      IP<-IP[-z,]
      Mean<-log10(IP$Mean)
      Sd<-log10(IP$Sd)
      out<-list(Mean,Sd)
      return(out)
    }
    f<-function(bam)
    {
      w<-bam[se]
      for (i in 2:length(ind))
      {
        se <- seq(i, n, len)
        w<- cbind(w, bam[se])
      }
      return(w)
    }
    unit_ms_f<-function(group)
    {
      m<-matrix()
      v<-matrix(data=0,nrow=length(se),ncol = length(ind))
      group<-as.matrix(group)
      for(i in 1:ncol(group))
      {
        
        m<-f(group[,i])
        v<-m+v
      }
      v<-v/ncol(group)
      Mean<-apply(v,1,mean)
      Sd<-apply( v,1,sd)
      Input<-cbind(Mean,Sd)
      Input<-as.data.frame(Input)
      z<-which(Input$Mean==0)
      Input<-Input[-z,]
      Mean<-log10(Input$Mean)
      Sd<-log10(Input$Sd)
      out<-list(Mean,Sd)
      return(out)
    }
    IPa<-data.frame()
    ID<-vector()
    group1<-as.matrix(group1)
    for(i in 1:ncol(group1))
    {
      outa<-ms_f(group1[,i])
      Mean<-outa[[1]]
      Sd<-outa[[2]]
      id<-rep(group1_name[i],length(Sd))
      ID<-c(ID,id)
      id<-vector()
      IP<-cbind(Mean,Sd)
      IPb<-IP
      IPa<-rbind(IPa,IPb)
      IP<-matrix()
    }
    outb<-unit_ms_f(group2)
    Mean<-outb[[1]]
    Sd<-outb[[2]]
    Input<-cbind(Mean,Sd)
    Input<-as.data.frame(Input)
    Id<-rep(group2_name,length(Input$Sd))
    ID<-c(ID,Id)
    com<-rbind(IPa,Input)
    com<-cbind(com,ID)
    com<-as.data.frame(com)
    cod<-rep(condition,length(com$ID))
    com<-cbind(com,cod)
    com<-as.data.frame(com)
    names(com)<-c("Mean","Sd","ID","Condition")
    return(com)
  }
  if((length(reference_IP_groupname)==0)&(length(reference_Input_groupname)==0))
  {
    Group_IP<-sa[,(1:length(IP_groupname))]
    Group_Input<-sa[,-(1:length(IP_groupname))]
    coma<-ms(Group_IP,Group_Input,IP_groupname,"unified_Input",paste(condition1,"condition"))
    comb<-ms(Group_Input,Group_IP,Input_groupname,"unified_IP",paste(condition1,"condition"))
    m_p1<-
      ggplot(coma,aes(Mean,Sd,colour=ID))+facet_grid(~Condition)+
      geom_smooth(aes(group = ID),span = 0.5)+
      geom_point(alpha = I(1 /150),size=0.2)+
      theme(title= element_text(size=12,color="black"))+
      labs(title="Mean-variance relationship within IP group under one condition")
    m_p2<-
      ggplot(comb,aes(Mean,Sd,colour=ID))+facet_grid(~Condition)+
      geom_smooth(aes(group = ID),span = 0.5)+
      geom_point(alpha = I(1 /150),size=0.2)+
      theme(title= element_text(size=12,color="black"))+
      labs(title="Mean-variance relationship within Input group under one condition")
    .multiplot(m_p1,m_p2,cols = 1)
  }
  else if((length(reference_IP_groupname)!=0)&(length(reference_Input_groupname)!=0))
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
    com1<-ms(group_IP,group_Input,IP_groupname,"unified_Input",paste(condition1,"condition"))
    com2<-ms(ref_group_IP,ref_group_Input,reference_IP_groupname,"unified_Input",paste(condition2,"condition"))
    coma<-rbind(com1,com2)
    coma<-as.data.frame(coma)
    com3<-ms(group_Input,group_IP,Input_groupname,"unified_IP",paste(condition1,"condition"))
    com4<-ms(ref_group_Input,ref_group_IP,reference_Input_groupname,"unified_IP",paste(condition2,"condition"))
    comb<-rbind(com3,com4)
    comb<-as.data.frame(comb)
    m_p1<-
      ggplot(coma,aes(Mean,Sd,colour=ID))+facet_grid(~Condition)+
      geom_smooth(aes(group = ID),span = 0.5)+     
      geom_point(alpha = I(1 /150),size=0.2)+
      theme(title= element_text(size=12,color="black"))+
      labs(title="Mean-variance relationship within IP group under two condition")
    m_p2<-
      ggplot(comb,aes(Mean,Sd,colour=ID))+facet_grid(~Condition)+
      geom_smooth(aes(group = ID),span = 0.5)+
      geom_point(alpha = I(1 /150),size=0.2)+theme(title= element_text(size=12,color="black"))+
      labs(title="Mean-variance relationship within Input group under two condition")
    .multiplot(m_p1,m_p2,cols = 2)
  }
}