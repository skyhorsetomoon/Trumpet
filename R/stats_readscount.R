.stats_readscount<-function(result, IP_BAM, Input_BAM, contrast_IP_BAM, contrast_Input_BAM)
{
  s<-result[[1]]
  ind<-unique(s$pos)
  len <- length(ind)
  n <- nrow(s)
  se <- seq(1, n, len)
  sa<-s[,-(1:2)]
  k<-c(0,1,10,100,1000,10000,100000)
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
  
  f<-function(bam)
  {
    s1<-.singleBAMreads(bam, se, len, n)
    vbam<-sapply(s1,unlist)
    r<-apply(s1,1,sum)
    genzero<-vector(mode="numeric",length=0)
    perzero<-vector(mode="numeric",length=0)
    for(i in 2:length(k))
    {
      if(max( vbam)>=k[i-1])
      {
        z<-which(vbam<k[i]&vbam>=k[i-1])
        perzero[i]<-length(z)/length(vbam)
      }
      if(max( vbam)<k[i-1])
      {
        perzero[i]<-0
      }
      if(max(r)>=k[i-1])
      {
        z1<-which(r<k[i]&r>=k[i-1])
        if(length(z1)!=0)
        {
          if(length(z1)==1)
          {
            gbam<-s1[(z1),]
            gbam<-t(as.matrix(gbam))
            genzero[i]<-nrow(gbam)/nrow(s1)
          }
          else
          {
            gbam<-s1[(z1),]
            gbam<-as.matrix(gbam)
            genzero[i]<-nrow(gbam)/nrow(s1)
          }
        }
        if(length(z1)==0)
        {
          genzero[i]<-0
        }
      }
      if(max( vbam)<k[i-1])
      {
        genzero[i]<-0
      }
    }
    genzero<-genzero[-1]
    perzero<-perzero[-1]
    unit<-list(genzero,perzero)
    return(unit)
  }
  
  p<-list()
  q<-list()
  for(i in seq_len(ncol(sa)))
  {
    p<-f(sa[,i])
    q<-c(q,p)
  }
  
  bin_p<-data.frame()
  trans_z<-data.frame()
  se1<-seq(1,length(q),2)
  for(i in seq_len(length(se1)))
  {
    geom_z<-t(q[[se1[i]]])
    trans_z<-rbind(trans_z,geom_z)
  }
  se2<-seq(2,length(q),2)
  for(i in seq_len(length(se2)))
  {
    bin_z<-t(q[[se2[i]]])
    bin_p<-rbind(bin_p,bin_z)
  }
  name<-colnames(s)
  rownames(trans_z)<-name[-(1:2)]
  rownames(bin_p)<-name[-(1:2)]
  IP_name<-c(IP_groupname,reference_IP_groupname)
  Input_name<-c(Input_groupname,reference_Input_groupname)
  trans_z<-as.data.frame(trans_z)
  bin_p<-as.data.frame(bin_p)
  IP_geom<-trans_z[IP_name,]
  IP_bin<-bin_p[IP_name,]
  Input_geom<-trans_z[Input_name,]
  Input_bin<-bin_p[Input_name,]
  p_geom<-vector()
  for(i in seq_len(nrow(IP_geom)))
  {
    p_geom<-rbind(p_geom,paste(round(IP_geom[i,]*100,3),"%"))
  }
  p_bin<-vector()
  for(i in seq_len(nrow(IP_bin)))
  {
    p_bin<-rbind(p_bin,paste(round(IP_bin[i,]*100,3),"%"))
  }
  rownames(p_geom)<-IP_name
  colnames(p_geom)<-c("0","1~10","10~100","100~1000","1000~10000","10000~100000")
  rownames(p_bin)<-IP_name
  colnames(p_bin)<-c("0","1~10","10~100","100~1000","1000~10000","10000~100000")
  
  pt_geom<-vector()
  for(i in seq_len(nrow(Input_geom)))
  {
    pt_geom<-rbind(pt_geom,paste(round(Input_geom[i,]*100,3),"%"))
  }
  rownames(pt_geom)<-Input_name
  colnames(pt_geom)<-c("0","1~10","10~100","100~1000","1000~10000","10000~100000")
  
  pt_bin<-vector()
  for(i in seq_len(nrow(Input_bin)))
  {
    pt_bin<-rbind(pt_bin,paste(round(Input_bin[i,]*100,3),"%"))
  }
  rownames(pt_bin)<-Input_name
  colnames(pt_bin)<-c("0","1~10","10~100","100~1000","1000~10000","10000~100000")
  readscount_summary<-result[[2]]
  readscount_summary<-as.data.frame(readscount_summary)
  transform_table<-result[[3]]
  p_statics<-list(transform_table,readscount_summary,p_geom,pt_geom,p_bin,pt_bin)
  return(p_statics)
}
