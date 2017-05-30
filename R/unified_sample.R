.unified_sample<-function(group,se,ind,len,n)
{
  m<-matrix(nrow = length(se),ncol = length(ind))
  v<-matrix(data=0,nrow=length(se),ncol = length(ind))
  group<-as.matrix(group)
  for(i in seq_len(ncol(group))){
    m<-.singleBAMreads(group[,i],se,len,n)
    v<-m+v
  }
  v<-v/length(ncol(group))  
  return(v)
}
