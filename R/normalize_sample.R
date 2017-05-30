.normalize_sample<-function(s1){
  row.sum<-rowSums(s1)
  z<-which(row.sum>10)
  s2<-s1[z,]
  row.mean<-rowMeans(s2)
  for(i in seq_len(length(row.mean))){
    if(row.mean[i]<2)
      row.mean[i]<-2
  }
  s3<-apply(s2,2,function(x,a)x/a, a=row.mean)
  s3<-as.matrix(s3)
  return(s3)
}
