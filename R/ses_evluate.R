.ses_evluate<-function(result, IP_BAM, Input_BAM, contrast_IP_BAM, contrast_Input_BAM, condition1, condition2)
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
  Norm_bam<-function(bam)
  {
    s1<-.singleBAMreads(bam, se, len, n)
    row.sum<-apply(s1,1,sum)
    z<-which(row.sum>10)
    s2<-s1[z,]
    row.mean<-apply(s2,1,mean)
    for(i in seq_len(length(row.mean)))
    {
      if(row.mean[i]<2)
        row.mean[i]<-2
      
    }
    s3<-apply(s2,2,function(x,a)x/a,a=row.mean)
    s3<-as.matrix(s3)
    bam<-sapply(t(s3),unlist)
    return(bam)
  }##get IP
  Unit_bam<-function(group)
  {
    m<-matrix(nrow=length(se),ncol = length(ind))
    v<-matrix(data=0,nrow=length(se),ncol = length(ind))
    group<-as.matrix(group)
    for(i in 1:ncol(group))
    {
      
      m<-.singleBAMreads(group[,i],se, len, n)
      v<-m+v
    }
    v<-v/length(Input_groupname)
    row.sum<-apply(v,1,sum)
    z<-which(row.sum>10)
    Input1<-v[z,]
    row.mean<-apply(Input1,1,mean)
    for(i in seq_len(length(row.mean)))
    {
      if(row.mean[i]<2)
        row.mean[i]<-2
      
    }
    Input2<-apply(Input1,2,function(x,a)x/a,a=row.mean)
    Input2<-as.matrix(Input2)
    Input<-sapply(t(Input2),unlist)
    return(Input)
  }
  SES_IP<-function(group_bam,unit_bam,IP_group_name)
  {
    new<-data.frame()
    a<-vector(mode="numeric",length=0)
    b<-vector(mode="numeric",length=0)
    z<-vector(mode="numeric",length=0)
    ID<-vector()
    for(i in seq_len(length(IP_group_name)))
    {
      
      bam<-vector(mode="numeric",length=0)
      bam<-Norm_bam(group_bam[,i])
      Input1<-unit_bam
      M<-max(length(bam),length(Input1))
      zeroa1<-rep(0,(M-length(bam)))
      zeroa2<-rep(0,(M-length(Input1)))
      ip1<-c(zeroa1,bam)
      InPut1<-c(zeroa2,Input1)
      v1<-sort(ip1)
      v2<-sort(InPut1)
      x<-v1-min(v1)
      ip<-x/sum(x)
      cum_bam<-vector(mode="numeric",length=0)
      cum_bam<-cumsum(ip)
      newpos<-1:length(ip)
      pos<-newpos/length(ip)
      x1<-v2-min(v2)
      Input<-x1/sum(x1)
      unified_Input<-vector(mode="numeric",length=0)
      unified_Input<-cumsum(Input)
      c<-(unified_Input-cum_bam)
      a[i]<-pos[which.max(c)]
      z[i]<-max(c)
      com<-cbind(pos,cum_bam,unified_Input)
      com<-as.data.frame(com)
      new1<-melt(data=com,id="pos",value.name = "pro")
      new1<-as.data.frame(new1)
      var<-rep(c(IP_group_name[i],"unified_Input"),c(nrow(new1[(new1$variable)=="cum_bam",]),nrow(new1[(new1$variable)=="unified_Input",])))
      new1$variable<-var
      b[i]<-length(new1$pos)
      new2<-new1
      new1<-data.frame()
      new<-rbind(new,new2,new1)
      ID1<-rep(IP_group_name[i],b[i])
      ID2<-ID1
      ID1<-vector()
      ID<-c(ID,ID2,ID1)
    }
    new<-cbind(new,ID)
    new<-as.data.frame(new)
    Scale_factor<-z
    Scale_factor<-round(Scale_factor,2)
    p<-vector()
    for(i in seq_len(length(IP_group_name)))
    {
      
      p[i]<-paste(round((1-a[i])*100,2),"%")
    }
    Sample<-IP_group_name
    Enriched_percent<-p
    Enrich_table<-cbind(Sample,Enriched_percent,Scale_factor)
    Enrich_table<-as.data.frame(Enrich_table)
    unit<-list(new,a,Enrich_table)
    return(unit)
  }
  SES_Input<-function(group_Input,unit_IP,Input_group_name)
  {
    new<-data.frame()
    for(i in seq_len(length(Input_group_name)))
    {
      bam<-vector(mode="numeric",length=0)
      bam<-Norm_bam(group_Input[,i])
      M<-max(length(bam),length(unit_IP))
      zeroa1<-rep(0,(M-length(bam)))
      zeroa2<-rep(0,(M-length(unit_IP)))
      input1<-c(zeroa1,bam)
      IP1<-c(zeroa2,unit_IP)
      v1<-sort(input1)
      v2<-sort(IP1)
      x<-v1-min(v1)
      input<-x/sum(x)
      cum_bam<-vector(mode="numeric",length=0)
      cum_bam<-cumsum(input)
      newpos<-1:length(input)
      pos<-newpos/length(input)
      x1<-v2-min(v2)
      IP<-x1/sum(x1)
      unified_IP<-vector(mode="numeric",length=0)
      unified_IP<-cumsum(IP)
      com<-cbind(pos,cum_bam,unified_IP)
      com<-as.data.frame(com)
      new1<-melt(data=com,id="pos",value.name = "pro")
      new1<-as.data.frame(new1)
      var<-rep(c(Input_group_name[i],"unified_IP"),c(nrow(new1[(new1$variable)=="cum_bam",]),nrow(new1[(new1$variable)=="unified_IP",])))
      new1$variable<-var
      new2<-new1
      new1<-data.frame()
      new<-rbind(new,new2,new1)
    }
    new<-as.data.frame(new)
    return(new)
  }
  if((length(reference_IP_groupname)==0)&(length(reference_Input_groupname)==0))
  {
    Group_IP<-sa[,(1:length(IP_groupname))]
    Group_IP<-as.matrix(Group_IP)
    Group_Input<-sa[,-(1:length(IP_groupname))]
    Group_Input<-as.matrix(Group_Input)
    Unit_Input<-Unit_bam(Group_Input)
    Unit_IP<-Unit_bam(Group_IP)
    out<-SES_IP(Group_IP,Unit_Input,IP_groupname)
    newa<-out[[1]]
    a<-out[[2]]
    Enrich_table<-out[[3]]
    new1<-SES_Input(Group_Input,Unit_IP,Input_groupname)
    vline<-data.frame(ID=IP_groupname,pos=a)
    
    newa<-as.data.frame(newa)
    pos<-newa$pos
    pro<-newa$pro
    variable<-newa$variable
    p1<-ggplot(data=newa,aes(x=pos,y=pro,colour=variable))+
      geom_line()+facet_wrap(~ID)+
      geom_vline(aes(xintercept = pos),vline)+
      theme(title= element_text(size=10,color="black"))+
      labs(x="Percentage of bins",y="Percentage of tags",title = paste("IP cumulative percentage enrichment within",condition1,"condition"))
    
    new1<-as.data.frame(new1)
    pos<-new1$pos
    pro<-new1$pro
    variable<-new1$variable
    p2<-
      ggplot(data=new1,aes(x=pos,y=pro,colour=variable))+
      geom_line()+
      theme(title= element_text(size=10,color="black"))+
      labs(x="Percentage of bins",y="Percentage of tags",title =paste("Input cumulative percentage enrichment within",condition1,"condition"))
    .multiplot(p1, p2, Enrich_table,cols=1)
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
    Unit_Input<-Unit_bam(group_Input)
    Unit_IP<-Unit_bam(group_IP)
    ref_unit_Input<-Unit_bam(ref_group_Input)
    ref_unit_IP<-Unit_bam(ref_group_IP)
    
    out<-SES_IP(group_IP,Unit_Input,IP_groupname)
    refer_out<-SES_IP(ref_group_IP,ref_unit_Input,reference_IP_groupname)
    newa<-out[[1]]
    a<-out[[2]]
    Enrich_table<-out[[3]]
    newb<-refer_out[[1]]
    b<-refer_out[[2]]
    refer_Enrich_table<-refer_out[[3]]
    new1<-SES_Input(group_Input,Unit_IP,Input_groupname)
    new2<-SES_Input(ref_group_Input,ref_unit_IP,reference_Input_groupname)
    vline1<-data.frame(ID=IP_groupname,pos=a)
    
    newa<-as.data.frame(newa)
    pos<-newa$pos
    pro<-newa$pro
    variable<-newa$variable
    p1<-
      ggplot(data=newa,aes(x=pos,y=pro,colour=variable))+
      geom_line()+facet_wrap(~ID)+
      geom_vline(aes(xintercept = pos),vline1)+
      theme(title= element_text(size=10,color="black"))+
      labs(x="Percentage of bins",y="Percentage of tags",title = paste("IP cumulative percentage enrichment within",condition1,"condition"))
    
    new1<-as.data.frame(new1)
    pos<-new1$pos
    pro<-new1$pro
    variable<-new1$variable
    p2<-
      ggplot(data=new1,aes(x=pos,y=pro,colour=variable))+geom_line()+
      theme(title= element_text(size=10,color="black"))+
      labs(x="Percentage of bins",y="Percentage of tags",title = paste("Input cumulative percentage enrichment within",condition1,"condition"))
    vline2<-data.frame(ID=reference_IP_groupname,pos=b)
    
    newb<-as.data.frame(newb)
    pos<-newb$pos
    pro<-newb$pro
    variable<-newb$variable
    p3<-
      ggplot(data=newb,aes(x=pos,y=pro,colour=variable))+
      geom_line()+facet_wrap(~ID)+
      geom_vline(aes(xintercept = pos),vline2)+
      theme(title= element_text(size=10,color="black"))+
      labs(x="Percentage of bins",y="Percentage of tags",title = paste("refer_IP cumulative percentage enrichment within",condition2,"condition"))
    
    new2<-as.data.frame(new2)
    pos<-new2$pos
    pro<-new2$pro
    variable<-new2$variable
    p4<-
      ggplot(data=new2,aes(x=pos,y=pro,colour=variable))+geom_line()+
      theme(title= element_text(size=10,color="black"))+
      labs(x="Percentage of bins",y="Percentage of tags",title = paste("refer_Input cumulative percentage enrichment within",condition2,"condition"))
    .multiplot(p1,p3,cols=1)
    .multiplot( p2,p4,cols=1)
    tab<-list(Enrich_table,refer_Enrich_table)
    return(tab)
  }
}
