.Norm_bam <- function(bam, se, len, n) {
  s1 <- .singleBAMreads(bam, se, len, n)
  s3 <- .normalize_sample(s1)
  bam <- sapply(t(s3), unlist)
  return(bam)
}  ##get IP


.Unit_bam <- function(group,se,ind,len,n) {
  v <- .unified_sample(group,se,ind,len,n)
  Input1 <- .normalize_sample(v)
  Input <- sapply(t(Input1), unlist)
  return(Input)
}

.SES_IP <- function(group_bam, unit_bam, IP_group_name, se, len, n) {
  new <- data.frame()
  a <- vector(mode = "numeric", length = 0)
  b <- vector(mode = "numeric", length = 0)
  z <- vector(mode = "numeric", length = 0)
  ID <- vector()
  for (i in seq_len(length(IP_group_name))) {
    
    bam <- vector(mode = "numeric", length = 0)
    bam <- .Norm_bam(group_bam[, i], se, len, n)
    Input1 <- unit_bam
    if(length(bam)==length(Input1)){
      ip1 <- bam
      InPut1 <- Input1
    }
    M <- max(length(bam), length(Input1))
    if(length(bam) < M && length(Input1)==M){
      bam_sample <- sample(bam, size = (M-length(bam)),replace = TRUE)
      ip1 <- c(bam_sample, bam)
      InPut1 <-  Input1
    }
    if(length(Input1) < M && length(bam)==M){
      Input1_sample <- sample(Input1, size = (M-length(Input1)),replace = TRUE)
      ip1 <- bam  
      InPut1 <- c(Input1_sample, Input1)
    }
    v1 <- sort(ip1)
    v2 <- sort(InPut1)
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
    z[i] <- max(c)
    com <- cbind(pos, cum_bam, unified_Input)
    com <- as.data.frame(com)
    new1 <- melt(data = com, id = "pos", value.name = "pro")
    new1 <- as.data.frame(new1)
    var <- rep(c(IP_group_name[i], "unified_Input"), c(nrow(new1[(new1$variable) == 
                                                                   "cum_bam", ]), nrow(new1[(new1$variable) == "unified_Input", 
                                                                                            ])))
    new1$variable <- var
    b[i] <- length(new1$pos)
    new2 <- new1
    new1 <- data.frame()
    new <- rbind(new, new2, new1)
    ID1 <- rep(IP_group_name[i], b[i])
    ID2 <- ID1
    ID1 <- vector()
    ID <- c(ID, ID2, ID1)
  }
  new <- cbind(new, ID)
  new <- as.data.frame(new)
  Scale_factor <- z
  Scale_factor <- round(Scale_factor, 2)
  p <- vector()
  for (i in seq_len(length(IP_group_name))) {
    
    p[i] <- paste(round((1 - a[i]) * 100, 2), "%")
  }
  Sample <- IP_group_name
  Enrichment_region <- p
  Enrich_table <- cbind(Sample, Enrichment_region, Scale_factor)
  Enrich_table <- as.data.frame(Enrich_table)
  unit <- list(new, a, Enrich_table)
  return(unit)
}

.SES_Input <- function(group_Input, unit_IP, Input_group_name,se,len,n) {
  new <- data.frame()
  for (i in seq_len(length(Input_group_name))) {
    bam <- vector(mode = "numeric", length = 0)
    bam <- .Norm_bam(group_Input[, i],se,len,n)
    if(length(bam)==length(unit_IP)){
      input1 <- bam
      IP1 <- unit_IP
    }
    M <- max(length(bam), length(unit_IP))
    if(length(bam) < M && length(unit_IP)==M){
      bam_sample <- sample(bam, size = (M-length(bam)),replace = TRUE)
      input1 <- c(bam_sample, bam)
      IP1 <-  unit_IP
    }
    if(length(unit_IP) < M && length(bam)==M){
      IP_sample <- sample(unit_IP, size = (M-length(unit_IP)),replace = TRUE)
      input1 <- bam  
      IP1 <- c(IP_sample, unit_IP)
    }
    v1 <- sort(input1)
    v2 <- sort(IP1)
    x <- v1 - min(v1)
    input <- x/sum(x)
    cum_bam <- vector(mode = "numeric", length = 0)
    cum_bam <- cumsum(input)
    newpos <- 1:length(input)
    pos <- newpos/length(input)
    x1 <- v2 - min(v2)
    IP <- x1/sum(x1)
    unified_IP <- vector(mode = "numeric", length = 0)
    unified_IP <- cumsum(IP)
    com <- cbind(pos, cum_bam, unified_IP)
    com <- as.data.frame(com)
    new1 <- melt(data = com, id = "pos", value.name = "pro")
    new1 <- as.data.frame(new1)
    var <- rep(c(Input_group_name[i], "unified_IP"), c(nrow(new1[(new1$variable) == 
                                                                   "cum_bam", ]), nrow(new1[(new1$variable) == "unified_IP", 
                                                                                            ])))
    new1$variable <- var
    new2 <- new1
    new1 <- data.frame()
    new <- rbind(new, new2, new1)
  }
  new <- as.data.frame(new)
  return(new)
}

.ESES_evaluate <- function(result, IP_BAM, Input_BAM, contrast_IP_BAM, contrast_Input_BAM, 
                         condition1, condition2) {
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
    Unit_Input <- .Unit_bam(Group_Input,se,ind,len,n)
    Unit_IP <- .Unit_bam(Group_IP,se,ind,len,n)
    out <- .SES_IP(Group_IP, Unit_Input, IP_groupname,se,len,n)
    newa <- out[[1]]
    a <- out[[2]]
    Enrich_table <- out[[3]]
    new1 <- .SES_Input(Group_Input, Unit_IP, Input_groupname,se,len,n)
    vline <- data.frame(ID = IP_groupname, pos = a)
    
    colnames(newa)<-c("pos","Sample","pro","ID")
    pos <- newa$pos
    pro <- newa$pro
    Sample <- newa$Sample
    p1 <- ggplot(data = newa, aes(x = pos, y = pro, colour = Sample )) + 
      geom_line() +
      facet_grid(ID~.) +
      geom_vline(aes(xintercept = pos), vline) +
      theme(axis.title.x =element_text(size=12), axis.title.y=element_text(size=12),
            title = element_text(size = 12),
            plot.title = element_text(hjust = 0.5),
            legend.key.height=unit(0.5,'cm'),
            legend.key.width=unit(0.25,'cm'),
            legend.text=element_text(size=10),
            legend.title=element_text(size=10))+ 
      labs(x = "Percentage of Bins", y = "Percentage of Reads",title = paste("Cumulative percentage enrichment of IP within" , condition1))
    
    colnames(new1)<-c("pos","Sample","pro")
    new1 <- as.data.frame(new1)
    pos <- new1$pos
    pro <- new1$pro
    Sample <- new1$Sample
    p2 <- ggplot(data = new1, aes(x = pos, y = pro, colour = Sample)) + 
      geom_line() +
      theme(axis.title.x =element_text(size=12), axis.title.y=element_text(size=12),
            title = element_text(size = 12),
            plot.title = element_text(hjust = 0.5),
            legend.key.height=unit(0.5,'cm'),
            legend.key.width=unit(0.25,'cm'),
            legend.text=element_text(size=10),
            legend.title=element_text(size=10))+ 
      labs(x = "Percentage of Bins", y = "Percentage of Reads",title = paste("Cumulative percentage enrichment of Input within" , condition1))
    
    .multiplot(p1, p2, cols = 1)
    colnames(Enrich_table) <- c("Sample ID", "Percent of Region Enriched with Signal", "Scale Factor")
    return(Enrich_table)
    
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
    Unit_Input <- .Unit_bam(group_Input,se,ind,len,n)
    Unit_IP <- .Unit_bam(group_IP,se,ind,len,n)
    ref_unit_Input <- .Unit_bam(ref_group_Input,se,ind,len,n)
    ref_unit_IP <- .Unit_bam(ref_group_IP,se,ind,len,n)
    
    out <- .SES_IP(group_IP, Unit_Input, IP_groupname,se,len,n)
    refer_out <- .SES_IP(ref_group_IP, ref_unit_Input, reference_IP_groupname,se,len,n)
    newa <- out[[1]]
    a <- out[[2]]
    Enrich_table <- out[[3]]
    newb <- refer_out[[1]]
    b <- refer_out[[2]]
    refer_Enrich_table <- refer_out[[3]]
    new1 <- .SES_Input(group_Input, Unit_IP, Input_groupname,se,len,n)
    new2 <- .SES_Input(ref_group_Input, ref_unit_IP, reference_Input_groupname,se,len,n)
    vline1 <- data.frame(ID = IP_groupname, pos = a)
    vline2 <- data.frame(ID = reference_IP_groupname, pos = b)
    
    colnames(newa)<-c("pos","Sample","pro","ID")
    pos <- newa$pos
    pro <- newa$pro
    Sample <- newa$Sample
    p1 <- ggplot(data = newa, aes(x = pos, y = pro, colour = Sample )) + 
      geom_line() +
      facet_grid(ID~.) +
      geom_vline(aes(xintercept = pos), vline1) +
      theme(axis.title.x =element_text(size=9), axis.title.y=element_text(size=9),
            title = element_text(size = 9),
            plot.title = element_text(hjust = 0.5),
            legend.key.height=unit(0.5,'cm'),
            legend.key.width=unit(0.25,'cm'),
            legend.text=element_text(size=9),
            legend.title=element_text(size=9))+ 
      labs(x = "Percentage of Bins", y = "Percentage of Reads",title = paste("Cumulative percentage enrichment of IP within" , condition1))
    
    colnames(new1)<-c("pos","Sample","pro")
    new1 <- as.data.frame(new1)
    pos <- new1$pos
    pro <- new1$pro
    Sample <- new1$Sample
    p2 <- ggplot(data = new1, aes(x = pos, y = pro, colour = Sample)) + 
      geom_line() +
      theme(axis.title.x =element_text(size=9), axis.title.y=element_text(size=9),
            title = element_text(size = 9),
            plot.title = element_text(hjust = 0.5),
            legend.key.height=unit(0.5,'cm'),
            legend.key.width=unit(0.25,'cm'),
            legend.text=element_text(size=9),
            legend.title=element_text(size=9))+ 
      labs(x = "Percentage of Bins", y = "Percentage of Reads",title = paste("Cumulative percentage enrichment of Input within" , condition1))
    
    
    colnames(newb)<-c("pos","Sample","pro","ID")
    newb <- as.data.frame(newb)
    pos <- newb$pos
    pro <- newb$pro
    Sample <- newb$Sample
    p3 <- ggplot(data = newb, aes(x = pos, y = pro, colour = Sample)) + 
      geom_line() + 
      facet_grid(ID~.) + 
      geom_vline(aes(xintercept = pos), vline2) + 
      theme(axis.title.x =element_text(size=9), axis.title.y=element_text(size=9),
            title = element_text(size = 9),
            plot.title = element_text(hjust = 0.5),
            legend.key.height=unit(0.5,'cm'),
            legend.key.width=unit(0.25,'cm'),
            legend.text=element_text(size=9),
            legend.title=element_text(size=9))+ 
      labs(x = "Percentage of Bins", y = "Percentage of Reads", title = paste("Cumulative percentage enrichment of Refer_IP within", 
                                                                              condition2))
    colnames(new2)<-c("pos","Sample","pro")
    new2 <- as.data.frame(new2)
    pos <- new2$pos
    pro <- new2$pro
    Sample <- new2$Sample
    p4 <- ggplot(data = new2, aes(x = pos, y = pro, colour = Sample)) + 
      geom_line() + 
      theme(axis.title.x =element_text(size=9), axis.title.y=element_text(size=9),
            title = element_text(size = 9),
            plot.title = element_text(hjust = 0.5),
            legend.key.height=unit(0.5,'cm'),
            legend.key.width=unit(0.25,'cm'),
            legend.text=element_text(size=9),
            legend.title=element_text(size=9))+ 
      labs(x = "Percentage of bins", y = "Percentage of Reads", title = paste("Cumulative percentage enrichment of Refer_Input within", condition2))
    .multiplot(p1, p3,p2, p4, cols = 2)
    tab <- rbind(Enrich_table, refer_Enrich_table)
    colnames(tab) <- c("Sample ID", "Percent of Region Enriched with Signal", "Scale Factor")
    return(tab)
  }
}
