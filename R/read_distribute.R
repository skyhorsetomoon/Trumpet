.read_cover <- function(read_count, ind) {
  r <- rowSums(read_count)
  read_count <- cbind(read_count, r)
  read_count <- as.data.frame(read_count)
  read_count <- read_count[order(read_count[, (length(ind) + 1)]),                            ]
  read_count <- read_count[, -(length(ind) + 1)]
  s3 <- .normalize_sample(read_count)
  b <- vector(mode = "numeric", length = 0)
  c <- vector(mode = "numeric", length = 0)
  d <- vector(mode = "numeric", length = 0)
  for (i in seq_len(length(ind))) {
    b[i] <- quantile(s3[, i], 0.25)
    c[i] <- quantile(s3[, i], 0.5)
    d[i] <- quantile(s3[, i], 0.75)
    
  }
  qv <- cbind(b, c, d)
  pos <- seq(0.025, 2.975, 0.05)
  dt <- cbind(pos, qv)
  dt <- as.data.frame(dt)
  colnames(dt) <- c("pos", "25%", "50%", "75%")
  df <- melt(dt, id = "pos")
  df <- data.frame(df)
  return(df)
}

.read_distribute <- function(GENOME, UCSC_TABLE_NAME, GENE_ANNO_GTF, TXDB, result, IP_BAM, Input_BAM, contrast_IP_BAM, 
                             contrast_Input_BAM, condition1, condition2) {
  
  # download the annotation
  if (suppressWarnings((!is.na(GENOME)) & (!is.na(UCSC_TABLE_NAME)) & 
                       is.na(TXDB) & is.na(GENE_ANNO_GTF))) {
    op <- options(warn = (-1))
    txdb = makeTxDbFromUCSC(genome = GENOME, tablename = UCSC_TABLE_NAME)
    options(op)
  }
  if (suppressWarnings(!is.na(GENE_ANNO_GTF) & is.na(TXDB))) {
    op <- options(warn = (-1))
    txdb <- makeTxDbFromGFF(GENE_ANNO_GTF, format = "gtf")
    options(op)
  }
  
  # use provided annotation data file
  if (suppressWarnings(!is.na(TXDB))) {
    txdb <- loadDb(TXDB)
  }
  
  ## get component
  utr5 <- fiveUTRsByTranscript(txdb, use.names=TRUE)
  cds <- cdsBy(txdb, by = "tx",  use.names=TRUE)
  utr3 <- threeUTRsByTranscript(txdb, use.names=TRUE)
  
  ## get componet name
  utr5_name <- names(utr5)
  utr5_name <- as.character(utr5_name)
  cds_name <- names(cds)
  cds_name <- as.character(cds_name)
  utr3_name <- names(utr3)
  utr3_name <- as.character(utr3_name)
  ## Select tx 
  s <- result[[1]]
  max_transcript_name <- as.character(unique(s$txid))
  ## Select component
  select_utr5_name <- as.numeric(match(utr5_name, max_transcript_name)) 
  select_utr5_name <- select_utr5_name[-which(is.na(select_utr5_name))] 
  select_utr5 <- utr5[select_utr5_name]
  select_cds_name <- as.numeric(match( cds_name, max_transcript_name))
  select_cds_name <- select_cds_name[-which(is.na(select_cds_name))] 
  select_cds <- cds[select_cds_name]
  select_utr3_name <- as.numeric(match(utr3_name, max_transcript_name))
  select_utr3_name <- select_utr3_name[-which(is.na(select_utr3_name))]
  select_utr3 <- utr3[select_utr3_name]
  
  select_name <- c(names(select_utr5), names(select_cds), names(select_utr3))
  last_select_name <- rownames(as.data.frame( which(table(select_name)==3))) 
  last_cds_num <- match(last_select_name, names(select_cds))
  last_cds <- select_cds[last_cds_num ]                     
  last_utr5_num <- match(last_select_name, names(select_utr5))
  last_utr5 <- select_utr5[last_utr5_num]
  last_utr3_num <- match(last_select_name, names(select_utr3))
  last_utr3 <- select_utr3[last_utr3_num]
  ## component size
  .sum_comb <- function(i, component){
    
    sum_width <- sum(width(component[[i]]))
    
  }
  
  utr5_size <- unlist(lapply(1:length(last_utr5), .sum_comb, component=last_utr5)) 
  cds_size <- unlist(lapply(1:length(last_cds), .sum_comb, component=last_cds)) 
  utr3_size <- unlist(lapply(1:length(last_utr3), .sum_comb, component=last_utr3)) 
  ## component size factor
  utr5.SF <- round((median(utr5_size)/median(cds_size)), 2)
  utr3.SF <- round((median(utr3_size)/median(cds_size)), 2) 
  
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
    group <- sa[, seq_len(length(IP_groupname) + length(Input_groupname))]
    p <- data.frame()
    q <- data.frame()
    group <- as.matrix(group)
    for (i in seq_len(ncol(group))) {
      d <- .singleBAMreads(group[, i], se, len, n)
      p <- .read_cover(d, ind)
      q <- rbind(q, p)
    }
    Group_IP <- q[seq_len(length(IP_groupname) * length(unique(q$pos)) * 
                            (nrow(q[q$pos == (unique(q$pos)[1]), ]))/ncol(group)), ]
    Group_Input <- q[-(seq_len(nrow(Group_IP))), ]
    fr_num <- length(unique(q$pos)) * (nrow(q[q$pos == (unique(q$pos)[1]), 
                                              ]))/ncol(group)
    group_IP <- group[, seq_len(length(IP_groupname))]
    group_IP <- as.matrix(group_IP)
    group_pt <- group[, -(seq_len(ncol(group_IP)))]
    group_pt <- as.matrix(group_pt)
    unit_IP <- .unified_sample(group_IP,se,ind,len,n)
    unit_Input <- .unified_sample(group_pt, se,ind,len,n)
    unified_IP <- .read_cover(unit_IP, ind)
    unified_Input <- .read_cover(unit_Input, ind)
    id_name1 <- c(IP_groupname, "Unified Input")
    ID1 <- rep(id_name1, rep(fr_num, length(id_name1)))
    id_name2 <- c(Input_groupname, "Unified IP")
    ID2 <- rep(id_name2, rep(fr_num, length(id_name2)))
    df1 <- rbind(Group_IP, unified_Input)
    df1 <- cbind(df1, ID1)
    df1 <- as.data.frame(df1)
    df2 <- rbind(Group_Input, unified_IP)
    df2 <- cbind(df2, ID2)
    df2 <- as.data.frame(df2)
    
    colnames(df1)<-c("pos","Quantile","value","ID1")
    pos <- df1$pos
    value <- df1$value
    Quantile <- df1$Quantile
    pos_utr5 <- unique(pos)[1:20]
    pos_cds <- unique(pos)[21:40]
    pos_utr3 <- unique(pos)[41:60]
    rescale_utr5 <- rescale(pos_utr5, c(1-utr5.SF,1), from = c(0,1))
    rescale_utr3 <- rescale(pos_utr3, to=c(2,2+utr3.SF), from = c(2,3))
    pos <- rep( c( rescale_utr5,  pos_cds, rescale_utr3), 3*length(unique(df1$ID1)))
    df1 <- cbind(pos, df1[,-1])
    
    p1 <- ggplot(df1, aes(pos, value, colour = Quantile)) + 
      geom_line() + 
      facet_grid(ID1~.) +  
      geom_vline(xintercept = 1:2, linetype = "dotted", colour="red") + 
      annotate("rect", xmin = pos[1], xmax = pos[length(unique(pos))/3], ymin = -0.2, ymax = -0.1, alpha = 0.99, colour = "black") + 
      annotate("rect", xmin = 1.006, xmax = 1.996, ymin = -0.25, ymax = -0.05, alpha = 0.99, colour = "black") + 
      annotate("rect", xmin = pos[(length(unique(pos))*2/3+1)]-0.007, xmax = pos[length(unique(pos))], ymin = -0.2, ymax = -0.1, alpha = 0.99, colour = "black") + 
      annotate("text", x = (pos[length(unique(pos))/6]+0.005), y= -0.4, label = "5'UTR", size = 3) + 
      annotate("text", x = 1.5, y = -0.4, label = "CDS", size = 3) + 
      annotate("text", x = 2.25, y = -0.4, label = "3'UTR", size = 3) + 
      theme(axis.title.x =element_text(size=9), axis.title.y=element_text(size=9),
            title = element_text(size = 9),
            plot.title = element_text(hjust = 0.5),
            legend.key.height=unit(0.5,'cm'),
            legend.key.width=unit(0.25,'cm'),
            legend.text=element_text(size=9),
            legend.title=element_text(size=9))+
      labs(title = paste("IP and Unified Input Reads Coverage within", condition1), x = "mRNA", y = "Normalized Reads Density") 
    
    colnames(df2)<-c("pos","Quantile","value","ID2")
    pos <- df2$pos
    value <- df2$value
    Quantile <- df2$Quantile
    pos_utr5 <- unique(pos)[1:20]
    pos_cds <- unique(pos)[21:40]
    pos_utr3 <- unique(pos)[41:60]
    rescale_utr5 <- rescale(pos_utr5, c(1-utr5.SF,1), from = c(0,1))
    rescale_utr3 <- rescale(pos_utr3, to=c(2,2+utr3.SF), from = c(2,3))
    pos <- rep( c( rescale_utr5,  pos_cds, rescale_utr3), 3*length(unique(df2$ID2)))
    df2 <- cbind(pos, df2[,-1])
    
    p2 <- ggplot(df2, aes(pos, value, colour = Quantile)) +
      geom_line() + 
      facet_grid(ID2~.) +  
      geom_vline(xintercept = 1:2, linetype = "dotted", colour="red") + 
      annotate("rect", xmin = pos[1], xmax = pos[length(unique(pos))/3], ymin = -0.2, ymax = -0.1, alpha = 0.99, colour = "black") + 
      annotate("rect", xmin = 1.006, xmax = 1.996, ymin = -0.25, ymax = -0.05, alpha = 0.99, colour = "black") + 
      annotate("rect", xmin = pos[(length(unique(pos))*2/3+1)]-0.007, xmax = pos[length(unique(pos))], ymin = -0.2, ymax = -0.1, alpha = 0.99, colour = "black") + 
      annotate("text", x = (pos[length(unique(pos))/6]+0.005), y= -0.4, label = "5'UTR", size = 3) + 
      annotate("text", x = 1.5, y = -0.4, label = "CDS", size = 3) + 
      annotate("text", x = 2.25, y = -0.4, label = "3'UTR", size = 3) + 
      theme(axis.title.x =element_text(size=9), axis.title.y=element_text(size=9),
            title = element_text(size = 9),
            plot.title = element_text(hjust = 0.5),
            legend.key.height=unit(0.5,'cm'),
            legend.key.width=unit(0.25,'cm'),
            legend.text=element_text(size=9),
            legend.title=element_text(size=9))+
      labs(title = paste("Input and Unified IP Reads Coverage within",condition1), x = "mRNA", y = "Normalized Reads Density") 
    
    .multiplot(p1, p2, cols = 1)
  } else if ((length(reference_IP_groupname) != 0) & (length(reference_Input_groupname) != 
                                                      0)) {
    
    group_one <- sa[, seq_len(length(IP_groupname) + length(Input_groupname))]
    p1 <- data.frame()
    q1 <- data.frame()
    group_one <- as.matrix(group_one)
    for (i in seq_len(ncol(group_one))) {
      q <- .singleBAMreads(group_one[, i], se, len, n)
      p1 <- .read_cover(q, ind)
      q1 <- rbind(q1, p1)
    }
    Group_IP <- q1[seq_len((length(IP_groupname) * length(unique(q1$pos)) * 
                              (nrow(q1[q1$pos == (unique(q1$pos)[1]), ]))/ncol(group_one))), 
                   ]
    Group_Input <- q1[-(seq_len(nrow(Group_IP))), ]
    group_two <- sa[, -(seq_len(ncol(group_one)))]
    p2 <- data.frame()
    q2 <- data.frame()
    group_two <- as.matrix(group_two)
    for (i in seq_len(ncol(group_two))) {
      p <- .singleBAMreads(group_two[, i], se, len, n)
      p2 <- .read_cover(p, ind)
      q2 <- rbind(q2, p2)
    }
    refer_Group_IP <- q2[seq_len((length(reference_IP_groupname) * 
                                    length(unique(q2$pos)) * (nrow(q2[q2$pos == (unique(q2$pos)[1]), 
                                                                      ]))/ncol(group_two))), ]
    refer_Group_Input <- q2[-(seq_len(nrow(refer_Group_IP))), ]
    fr_num1 <- length(unique(q1$pos)) * (nrow(q1[q1$pos == (unique(q1$pos)[1]), 
                                                 ]))/ncol(group_one)
    fr_num2 <- length(unique(q2$pos)) * (nrow(q2[q2$pos == (unique(q2$pos)[1]), 
                                                 ]))/ncol(group_two)
    group_one_IP <- group_one[, seq_len(length(IP_groupname))]
    group_one_IP <- as.matrix(group_one_IP)
    group_one_pt <- group_one[, -(seq_len(ncol(group_one_IP)))]
    group_two_IP <- group_two[, seq_len(length(reference_IP_groupname))]
    group_two_IP <- as.matrix(group_two_IP)
    group_two_pt <- group_two[, -(seq_len(ncol(group_two_IP)))]
    unit_IP <- .unified_sample(group_one_IP, se,ind,len,n)
    unit_Input <- .unified_sample(group_one_pt, se,ind,len,n)
    refer_unit_IP <- .unified_sample(group_two_IP, se,ind,len,n)
    refer_unit_pt <- .unified_sample(group_two_pt, se,ind,len,n)
    unified_IP <- .read_cover(unit_IP, ind)
    unified_Input <- .read_cover(unit_Input, ind)
    refer_unified_IP <- .read_cover(refer_unit_IP, ind)
    refer_unified_Input <- .read_cover(refer_unit_pt, ind)
    id_name1 <- c(IP_groupname, "Unified Input")
    ID1 <- rep(id_name1, rep(fr_num1, length(id_name1)))
    id_name2 <- c(Input_groupname, "Unified IP")
    ID2 <- rep(id_name2, rep(fr_num1, length(id_name2)))
    refer_id_name1 <- c(reference_IP_groupname, "Unified Input")
    refer_ID1 <- rep(refer_id_name1, rep(fr_num2, length(refer_id_name1)))
    refer_id_name2 <- c(reference_Input_groupname, "Unified IP")
    refer_ID2 <- rep(refer_id_name2, rep(fr_num2, length(refer_id_name2)))
    df1 <- rbind(Group_IP, unified_Input)
    df1 <- cbind(df1, ID1)
    df1 <- as.data.frame(df1)
    df2 <- rbind(Group_Input, unified_IP)
    df2 <- cbind(df2, ID2)
    df2 <- as.data.frame(df2)
    Df1 <- rbind(refer_Group_IP, refer_unified_Input)
    Df1 <- cbind(Df1, refer_ID1)
    Df1 <- as.data.frame(Df1)
    Df2 <- rbind(refer_Group_Input, refer_unified_IP)
    Df2 <- cbind(Df2, refer_ID2)
    Df2 <- as.data.frame(Df2)
    
    colnames(df1)<-c("pos","Quantile","value","ID1")
    pos <- df1$pos
    value <- df1$value
    Quantile <- df1$Quantile
    pos_utr5 <- unique(pos)[1:20]
    pos_cds <- unique(pos)[21:40]
    pos_utr3 <- unique(pos)[41:60]
    rescale_utr5 <- rescale(pos_utr5, c(1-utr5.SF,1), from = c(0,1))
    rescale_utr3 <- rescale(pos_utr3, to=c(2,2+utr3.SF), from = c(2,3))
    pos <- rep( c( rescale_utr5,  pos_cds, rescale_utr3), 3*length(unique(df1$ID1)))
    df1 <- cbind(pos, df1[,-1])
    
    p1 <- ggplot(df1, aes(pos, value, colour = Quantile)) + 
      geom_line() + 
      facet_grid(ID1~.) +  
      geom_vline(xintercept = 1:2, linetype = "dotted", colour="red") + 
      annotate("rect", xmin = pos[1], xmax = pos[length(unique(pos))/3], ymin = -0.2, ymax = -0.1, alpha = 0.99, colour = "black") + 
      annotate("rect", xmin = 1.006, xmax = 1.996, ymin = -0.25, ymax = -0.05, alpha = 0.99, colour = "black") + 
      annotate("rect", xmin = pos[(length(unique(pos))*2/3+1)]-0.007, xmax = pos[length(unique(pos))], ymin = -0.2, ymax = -0.1, alpha = 0.99, colour = "black") + 
      annotate("text", x = (pos[length(unique(pos))/6]+0.005), y= -0.4, label = "5'UTR", size = 3) + 
      annotate("text", x = 1.5, y = -0.4, label = "CDS", size = 3) + 
      annotate("text", x = 2.25, y = -0.4, label = "3'UTR", size = 3) + 
      theme(axis.title.x =element_text(size=9), axis.title.y=element_text(size=9),
            title = element_text(size = 9),
            plot.title = element_text(hjust = 0.5),
            legend.key.height=unit(0.5,'cm'),
            legend.key.width=unit(0.25,'cm'),
            legend.text=element_text(size=9),
            legend.title=element_text(size=9))+
      labs(title = paste("IP and Unified Input Reads Coverage within", condition1), x = "mRNA", y = "Normalized Reads Density") 
    
    colnames(df2)<-c("pos","Quantile","value","ID2")
    pos <- df2$pos
    value <- df2$value
    Quantile <- df2$Quantile
    pos_utr5 <- unique(pos)[1:20]
    pos_cds <- unique(pos)[21:40]
    pos_utr3 <- unique(pos)[41:60]
    rescale_utr5 <- rescale(pos_utr5, c(1-utr5.SF,1), from = c(0,1))
    rescale_utr3 <- rescale(pos_utr3, to=c(2,2+utr3.SF), from = c(2,3))
    pos <- rep( c( rescale_utr5,  pos_cds, rescale_utr3), 3*length(unique(df2$ID2)))
    df2 <- cbind(pos, df2[,-1])
    
    p2 <- ggplot(df2, aes(pos, value, colour = Quantile)) +
      geom_line() + 
      facet_grid(ID2~.) +  
      geom_vline(xintercept = 1:2, linetype = "dotted", colour="red") + 
      annotate("rect", xmin = pos[1], xmax = pos[length(unique(pos))/3], ymin = -0.2, ymax = -0.1, alpha = 0.99, colour = "black") + 
      annotate("rect", xmin = 1.006, xmax = 1.996, ymin = -0.25, ymax = -0.05, alpha = 0.99, colour = "black") + 
      annotate("rect", xmin = pos[(length(unique(pos))*2/3+1)]-0.007, xmax = pos[length(unique(pos))], ymin = -0.2, ymax = -0.1, alpha = 0.99, colour = "black") + 
      annotate("text", x = (pos[length(unique(pos))/6]+0.005), y= -0.4, label = "5'UTR", size = 3) + 
      annotate("text", x = 1.5, y = -0.4, label = "CDS", size = 3) + 
      annotate("text", x = 2.25, y = -0.4, label = "3'UTR", size = 3) + 
      theme(axis.title.x =element_text(size=9), axis.title.y=element_text(size=9),
            title = element_text(size = 9),
            plot.title = element_text(hjust = 0.5),
            legend.key.height=unit(0.5,'cm'),
            legend.key.width=unit(0.25,'cm'),
            legend.text=element_text(size=9),
            legend.title=element_text(size=9))+
      labs(title = paste("Input and Unified IP Reads Coverage within",condition1), x = "mRNA", y = "Normalized Reads Density") 
    
    colnames(Df1)<-c("pos","Quantile","value","Refer_ID1")
    pos <- Df1$pos
    value <- Df1$value
    Quantile <- Df1$Quantile
    pos_utr5 <- unique(pos)[1:20]
    pos_cds <- unique(pos)[21:40]
    pos_utr3 <- unique(pos)[41:60]
    rescale_utr5 <- rescale(pos_utr5, c(1-utr5.SF,1), from = c(0,1))
    rescale_utr3 <- rescale(pos_utr3, to=c(2,2+utr3.SF), from = c(2,3))
    pos <- rep( c( rescale_utr5,  pos_cds, rescale_utr3), 3*length(unique(Df1$Refer_ID1)))
    Df1 <- cbind(pos, Df1[,-1])
    
    p3 <- ggplot(Df1, aes(pos, value, colour = Quantile)) + 
      geom_line() + 
      facet_grid(Refer_ID1~.) + 
      geom_vline(xintercept = 1:2, linetype = "dotted", colour="red") + 
      annotate("rect", xmin = pos[1], xmax = pos[length(unique(pos))/3], ymin = -0.2, ymax = -0.1, alpha = 0.99, colour = "black") + 
      annotate("rect", xmin = 1.006, xmax = 1.996, ymin = -0.25, ymax = -0.05, alpha = 0.99, colour = "black") + 
      annotate("rect", xmin = pos[(length(unique(pos))*2/3+1)]-0.007, xmax = pos[length(unique(pos))], ymin = -0.2, ymax = -0.1, alpha = 0.99, colour = "black") + 
      annotate("text", x = (pos[length(unique(pos))/6]+0.005), y= -0.4, label = "5'UTR", size = 3) + 
      annotate("text", x = 1.5, y = -0.4, label = "CDS", size = 3) + 
      annotate("text", x = 2.25, y = -0.4, label = "3'UTR", size = 3) + 
      theme(axis.title.x =element_text(size=9), axis.title.y=element_text(size=9),
            title = element_text(size = 9),
            plot.title = element_text(hjust = 0.5),
            legend.key.height=unit(0.5,'cm'),
            legend.key.width=unit(0.25,'cm'),
            legend.text=element_text(size=9),
            legend.title=element_text(size=9))+
      labs(title = paste("Refer_IP and Unified Input Reads Coverage within", condition2),x = "mRNA", y = "Normalized Reads Density") 
    
    colnames(Df2)<-c("pos","Quantile","value","Refer_ID2")
    pos <- Df2$pos
    value <- Df2$value
    Quantile <- Df2$Quantile
    pos_utr5 <- unique(pos)[1:20]
    pos_cds <- unique(pos)[21:40]
    pos_utr3 <- unique(pos)[41:60]
    rescale_utr5 <- rescale(pos_utr5, c(1-utr5.SF,1), from = c(0,1))
    rescale_utr3 <- rescale(pos_utr3, to=c(2,2+utr3.SF), from = c(2,3))
    pos <- rep( c( rescale_utr5,  pos_cds, rescale_utr3), 3*length(unique(Df2$Refer_ID2)))
    Df2 <- cbind(pos, Df2[,-1])
    
    
    p4 <- ggplot(Df2, aes(pos, value, colour = Quantile)) +
      geom_line() + 
      facet_grid(Refer_ID2~.) +
      geom_vline(xintercept = 1:2, linetype = "dotted", colour="red") + 
      annotate("rect", xmin = pos[1], xmax = pos[length(unique(pos))/3], ymin = -0.2, ymax = -0.1, alpha = 0.99, colour = "black") + 
      annotate("rect", xmin = 1.006, xmax = 1.996, ymin = -0.25, ymax = -0.05, alpha = 0.99, colour = "black") + 
      annotate("rect", xmin = pos[(length(unique(pos))*2/3+1)]-0.007, xmax = pos[length(unique(pos))], ymin = -0.2, ymax = -0.1, alpha = 0.99, colour = "black") + 
      annotate("text", x = (pos[length(unique(pos))/6]+0.005), y= -0.4, label = "5'UTR", size = 3) + 
      annotate("text", x = 1.5, y = -0.4, label = "CDS", size = 3) + 
      annotate("text", x = 2.25, y = -0.4, label = "3'UTR", size = 3) + 
      theme(axis.title.x =element_text(size=9), axis.title.y=element_text(size=9),
            title = element_text(size = 9),
            plot.title = element_text(hjust = 0.5),
            legend.key.height=unit(0.5,'cm'),
            legend.key.width=unit(0.25,'cm'),
            legend.text=element_text(size=9),
            legend.title=element_text(size=9))+
      labs(title = paste("Refer_Input and Unified IP Reads Coverage within", condition2),x = "mRNA", y = "Normalized Reads Density")
    
    .multiplot(p1, p2, p3, p4, cols = 2)
  }
}
