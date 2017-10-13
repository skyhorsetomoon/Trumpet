.Unit_bam <- function(group, se,ind,len,n) {
  v <- .unified_sample(group,se,ind,len,n)
  Input1 <- .normalize_sample(v)
  Input <- sapply(t(Input1), unlist)
  return(Input)
}

.Norm_bam <- function(bam, se, len, n) {
  s1 <- .singleBAMreads(bam, se, len, n)
  s3 <- .normalize_sample(s1)
  bam <- sapply(t(s3), unlist)
  return(bam)
}  ##get IP

.SES_IPenrich <- function(group_bam, unit_bam, IP_group_name,se, len, n){
  bknuread <- vector(mode = "numeric", length = 0)
  signuread <- vector(mode = "numeric", length = 0)
  new <- data.frame()
  a <- vector(mode = "numeric", length = 0)
  b <- vector(mode = "numeric", length = 0)
  z <- vector(mode = "numeric", length = 0)
  ID <- vector()
  for (i in seq_len(length(IP_group_name))) {
    
    bam <- vector(mode = "numeric", length = 0)
    bam <- .Norm_bam(group_bam[, i], se, len, n)
    Input1 <- unit_bam
    bam_read <- .singleBAMreads(group_bam[, i], se, len, n)
    vbam_read <- .trans_readsvector(group_bam[, i], se, len, n)
    nubin <- length(vbam_read)
    vbam_read <- sort(vbam_read)
    M <- max(length(bam), length(Input1))
    zeroa1 <- rep(0, (M - length(bam)))
    zeroa2 <- rep(0, (M - length(Input1)))
    ip1 <- c(zeroa1, bam)
    InPut1 <- c(zeroa2, Input1)
    v1 <- sort(ip1)
    nuv1 <- length(v1)
    v2 <- sort(InPut1)
    nuv2 <- length(v2)
    x <- v1 - min(v1)
    ip <- x/sum(x)
    cum_bam <- vector(mode = "numeric", length = 0)
    cum_bam <- cumsum(ip)
    newpos <- 1:length(ip)
    pos <- newpos/length(ip)
    x1 <- v2 - min(v2)
    Input <- x1/sum(x1)
    unified_Input <- vector(mode = "numeric", length = 0)
    unified_Input <- cumsum(Input)
    c <- (unified_Input - cum_bam)
    a[i] <- pos[which.max(c)]
    bknuread[i] <- sum(vbam_read[1:(nubin*round(a[i],2))])
    signuread[i] <- sum(vbam_read[((nubin*round(a[i],2))):nubin])
    
    z[i] <- max(c)
  }
  backgroundread <- round(bknuread/(bknuread+signuread)* 100, 2) 
  Signal_readcount <- round(signuread/(bknuread+signuread)* 100, 2)
  Scale_factor <- z
  Scale_factor <- round(Scale_factor, 3)
  p <- vector()
  for (i in seq_len(length(IP_group_name))) {
    p[i] <- round((1 - a[i])* 100 , 2)
    
  }
  Enrichment_region <- p
  Enrich_table <- cbind(Enrichment_region, Signal_readcount,  Scale_factor)
  Enrichtable <- as.data.frame(Enrich_table)
  return(Enrichtable)
}

.SES_enrich_DF <- function(result, IP_BAM, Input_BAM, contrast_IP_BAM, contrast_Input_BAM){
  
  s <- result[[1]]
  ind <- unique(s$pos)
  len <- length(ind)
  n <- nrow(s)
  se <- seq(1, n, len)
  sa <- s[, -(1:2)]
  sample_name <- .get.sampleid(IP_BAM, Input_BAM, contrast_IP_BAM, contrast_Input_BAM)
  if (length(sample_name) == 2) {
    IP_groupname <- sample_name[[1]]
    Input_groupname <- sample_name[[2]]
    reference_IP_groupname <- NULL
    reference_Input_groupname <- NULL
  }
  if (length(sample_name) == 4) {
    IP_groupname <- sample_name[[1]]
    Input_groupname <- sample_name[[2]]
    reference_IP_groupname <- sample_name[[3]]
    reference_Input_groupname <- sample_name[[4]]
  }
  
  if ((length(reference_IP_groupname) == 0) & (length(reference_Input_groupname) == 
                                               0)) {
    Group_IP <- sa[, (seq_len(length(IP_groupname)))]
    Group_IP <- as.matrix(Group_IP)
    Group_Input <- sa[, -(seq_len(length(IP_groupname)))]
    Group_Input <- as.matrix(Group_Input)
    Unit_Input <- .Unit_bam(Group_Input, se,ind,len,n)
    ##get the user enrich data
    Enrichtable <- .SES_IPenrich(Group_IP, Unit_Input, IP_groupname, se, len, n)
    new_sample <- c(IP_groupname)
    
  } else if ((length(reference_IP_groupname) != 0) & (length(reference_Input_groupname) != 
                                                      0)) {
    group_IP <- sa[, (seq_len(length(IP_groupname)))]
    group_IP <- as.matrix(group_IP)
    Group_Input <- sa[, -(seq_len(length(IP_groupname)))]
    group_Input <- Group_Input[, -((length(Input_groupname) + 1):ncol(Group_Input))]
    group_Input <- as.matrix(group_Input)
    ref_group <- Group_Input[, -(seq_len(length(Input_groupname)))]
    ref_group_IP <- ref_group[, seq_len(length(reference_IP_groupname))]
    ref_group_IP <- as.matrix(ref_group_IP)
    ref_group_Input <- ref_group[, -(seq_len(length(reference_IP_groupname)))]
    ref_group_Input <- as.matrix(ref_group_Input)
    Unit_Input <- .Unit_bam(group_Input, se,ind,len,n)
    ref_unit_Input <- .Unit_bam(ref_group_Input, se,ind,len,n)
    
    EnrichTable <- .SES_IPenrich (group_IP, Unit_Input, IP_groupname, se, len, n)
    refer_Enrichtable <- .SES_IPenrich (ref_group_IP, ref_unit_Input, reference_IP_groupname, se, len, n)
    ##get the user enrich data
    Enrichtable <- rbind(EnrichTable, refer_Enrichtable)
    new_sample <- c(IP_groupname, reference_IP_groupname)
  }
  dircs<-system.file("extdata", "experiment_data.Rdata", package="Trumpet")
  load(dircs)
  enrichtable <- experiment_data
  old_sample <- rep("known data set",61)
  IP_sample <- c(old_sample, new_sample)
  enrichtable  <- cbind(old_sample,enrichtable)
  enrich_region <- enrichtable[,1:2]
  enrich_region <- as.data.frame(enrich_region)
  
  p1 <- ggplot(enrich_region, aes(x=Enrichment_region,fill=old_sample))+
    geom_density(alpha = 0.6,adjust=1)+
    geom_vline(data=Enrichtable, aes(xintercept=Enrichment_region,colour=new_sample), show.legend = TRUE, size=0.75)+
    theme(axis.title.x =element_text(size=9), axis.title.y=element_text(size=9),
          title = element_text(size = 9),
          legend.key.height=unit(0.5,'cm'),
          legend.key.width=unit(0.5,'cm'),
          legend.text=element_text(size=9),
          legend.title=element_text(size=9))+
    labs(x="Percentage of enrichment region(%)",  title="The distribution of IP samples' enrichment region ")
  
  signalread <- enrichtable[,c(1,3)]
  signalread <- as.data.frame(signalread)
  p2 <- ggplot(signalread, aes(x=Signal_readcount,fill=old_sample))+
    geom_density(alpha = 0.6,adjust=1)+
    geom_vline(data=Enrichtable, aes(xintercept=Signal_readcount,colour=new_sample), show.legend = TRUE, size=0.5)+
    theme(axis.title.x =element_text(size=9), axis.title.y=element_text(size=9),
          title = element_text(size = 9),
          legend.key.height=unit(0.5,'cm'),
          legend.key.width=unit(0.5,'cm'),
          legend.text=element_text(size=9),
          legend.title=element_text(size=9))+
    labs(x="Percentage of signal reads count(%)",  title="The distribution of IP samples' signal reads count ")
  scalefactor <- enrichtable[, c(1,4)]
  p3 <- ggplot(scalefactor, aes(x=Scale_factor,fill=old_sample))+
    geom_density(alpha = 0.6,adjust=1)+
    geom_vline(data=Enrichtable, aes(xintercept=Scale_factor,colour=new_sample), show.legend = TRUE, size=0.5)+
    theme(axis.title.x =element_text(size=9), axis.title.y=element_text(size=9),
          title = element_text(size = 9),
          legend.key.height=unit(0.5,'cm'),
          legend.key.width=unit(0.5,'cm'),
          legend.text=element_text(size=9),
          legend.title=element_text(size=9))+
    labs(x="Scale factor",  title="The distribution of IP samples' scale factor")
  .multiplot(p1, p2, p3, cols = 3)
  
}
