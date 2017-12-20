.normalize_sample<-function(s1,se, len, n){
  row.sum<-rowSums(s1)
  z<-which(row.sum>0)
  s2<-s1[z,]
  row.mean<-rowMeans(s2)
  ls2 <- which(row.mean<10)
  select_s2 <- s2[-ls2,]
  select_rowmean <- row.mean[-ls2]
  s3<-apply(select_s2,2,function(x,a)x/a, a=select_rowmean)
  ls3 <- sapply(t(s3), unlist)
  return(ls3)
}
