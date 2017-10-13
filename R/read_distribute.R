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
  colnames(dt) <- c("pos", "quantile=25%", "quantile=50%", "quantile=75%")
  df <- melt(dt, id = "pos")
  df <- data.frame(df)
  return(df)
}

.read_distribute <- function(result, IP_BAM, Input_BAM, contrast_IP_BAM, 
                             contrast_Input_BAM, condition1, condition2) {
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
    id_name1 <- c(IP_groupname, "unified Input")
    ID1 <- rep(id_name1, rep(fr_num, length(id_name1)))
    id_name2 <- c(Input_groupname, "unified IP")
    ID2 <- rep(id_name2, rep(fr_num, length(id_name2)))
    df1 <- rbind(Group_IP, unified_Input)
    df1 <- cbind(df1, ID1)
    df1 <- as.data.frame(df1)
    df2 <- rbind(Group_Input, unified_IP)
    df2 <- cbind(df2, ID2)
    df2 <- as.data.frame(df2)
    
    colnames(df1)<-c("pos","Quantile","value","ID1")
    pos <- df1$posample
    value <- df1$value
    Quantile <- df1$Quantile
    p1 <- ggplot(df1, aes(pos, value, colour = Quantile)) + 
      geom_line() + 
      facet_wrap(~ID1) +  
      annotate("text", x = 0.5, y = -0.4, label = "5'UTR", size = 4) + 
      annotate("text", x = 1.5, y = -0.4, label = "CDS", size = 4) + 
      annotate("text", x = 2.5, y = -0.4, label = "3'UTR", size = 4) + 
      geom_vline(xintercept = 1:2, linetype = "dotted") + 
      annotate("rect", xmin = 0.025, xmax = 0.975, ymin = -0.2, ymax = -0.1, alpha = 0.99, colour = "black") + 
      annotate("rect", xmin = 1.025, xmax = 1.975, ymin = -0.2, ymax = -0.1, alpha = 0.99, colour = "black") + 
      annotate("rect", xmin = 2.025, xmax = 2.975, ymin = -0.2, ymax = -0.1, alpha = 0.99, colour = "black") + 
      theme(axis.title.x =element_text(size=12), axis.title.y=element_text(size=12),
            title = element_text(size = 12),
            legend.key.height=unit(0.5,'cm'),
            legend.key.width=unit(0.25,'cm'),
            legend.text=element_text(size=10),
            legend.title=element_text(size=10))+
      labs(title = paste("IP and Unified Input Reads Coverage within", condition1), x = "mRNA", y = "Normalized Reads Density")
    
    colnames(df2)<-c("pos","Quantile","value","ID2")
    pos <- df2$pos
    value <- df2$value
    Quantile <- df2$Quantile
    p2 <- ggplot(df2, aes(pos, value, colour = Quantile)) +
      geom_line() + 
      facet_wrap(~ID2) +  
      annotate("text", x = 0.5, y = -0.4, label = "5'UTR", size = 4) + 
      annotate("text", x = 1.5, y = -0.4, label = "CDS", size = 4) + 
      annotate("text", x = 2.5, y = -0.4, label = "3'UTR", size = 4) + 
      geom_vline(xintercept = 1:2, linetype = "dotted") +
      annotate("rect", xmin = 0.025, xmax = 0.975, ymin = -0.2, ymax = -0.1, alpha = 0.99, colour = "black") +
      annotate("rect", xmin = 1.025, xmax = 1.975, ymin = -0.2, ymax = -0.1, alpha = 0.99, colour = "black") + 
      annotate("rect", xmin = 2.025, xmax = 2.975, ymin = -0.2, ymax = -0.1, alpha = 0.99, colour = "black") + 
      theme(axis.title.x =element_text(size=12), axis.title.y=element_text(size=12), 
            title = element_text(size = 12),
            legend.key.height=unit(0.5,'cm'),
            legend.key.width=unit(0.25,'cm'),
            legend.text=element_text(size=10),
            legend.title=element_text(size=10))+
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
    id_name1 <- c(IP_groupname, "unified Input")
    ID1 <- rep(id_name1, rep(fr_num1, length(id_name1)))
    id_name2 <- c(Input_groupname, "unified IP")
    ID2 <- rep(id_name2, rep(fr_num1, length(id_name2)))
    refer_id_name1 <- c(reference_IP_groupname, "unified Input")
    refer_ID1 <- rep(refer_id_name1, rep(fr_num2, length(refer_id_name1)))
    refer_id_name2 <- c(reference_Input_groupname, "unified IP")
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
    pos <- df1$posample
    value <- df1$value
    Quantile <- df1$Quantile
    p1 <- ggplot(df1, aes(pos, value, colour = Quantile)) + 
      geom_line() + 
      facet_grid(ID1~.) +  
      annotate("text", x = 0.5, y = -0.6, label = "5'UTR", size = 3) + 
      annotate("text", x = 1.5, y = -0.6, label = "CDS", size = 3) + 
      annotate("text", x = 2.5, y = -0.6, label = "3'UTR", size = 3) + 
      geom_vline(xintercept = 1:2, linetype = "dotted") + 
      annotate("rect", xmin = 0.025, xmax = 0.975, ymin = -0.3, ymax = -0.1, alpha = 0.99, colour = "black") + 
      annotate("rect", xmin = 1.025, xmax = 1.975, ymin = -0.3, ymax = -0.1, alpha = 0.99, colour = "black") + 
      annotate("rect", xmin = 2.025, xmax = 2.975, ymin = -0.3, ymax = -0.1, alpha = 0.99, colour = "black") + 
      theme(axis.title.x =element_text(size=9), axis.title.y=element_text(size=9),
            title = element_text(size = 9),
            legend.key.height=unit(0.5,'cm'),
            legend.key.width=unit(0.25,'cm'),
            legend.text=element_text(size=9),
            legend.title=element_text(size=9))+
      labs(title = paste("IP and Unified Input Reads Coverage within", condition1), x = "mRNA", y = "Normalized Reads Density") 
    
    colnames(df2)<-c("pos","Quantile","value","ID2")
    pos <- df2$pos
    value <- df2$value
    Quantile <- df2$Quantile
    p2 <- ggplot(df2, aes(pos, value, colour = Quantile)) +
      geom_line() + 
      facet_grid(ID2~.) +  
      annotate("text", x = 0.5, y = -0.6, label = "5'UTR", size = 3) + 
      annotate("text", x = 1.5, y = -0.6, label = "CDS", size = 3) + 
      annotate("text", x = 2.5, y = -0.6, label = "3'UTR", size = 3) + 
      geom_vline(xintercept = 1:2, linetype = "dotted") +
      annotate("rect", xmin = 0.025, xmax = 0.975, ymin = -0.3, ymax = -0.1, alpha = 0.99, colour = "black") +
      annotate("rect", xmin = 1.025, xmax = 1.975, ymin = -0.3, ymax = -0.1, alpha = 0.99, colour = "black") + 
      annotate("rect", xmin = 2.025, xmax = 2.975, ymin = -0.3, ymax = -0.1, alpha = 0.99, colour = "black") + 
      theme(axis.title.x =element_text(size=9), axis.title.y=element_text(size=9),
            title = element_text(size = 9),
            legend.key.height=unit(0.5,'cm'),
            legend.key.width=unit(0.25,'cm'),
            legend.text=element_text(size=9),
            legend.title=element_text(size=9))+
      labs(title = paste("Input and Unified IP Reads Coverage within",condition1), x = "mRNA", y = "Normalized Reads Density") 
    
    colnames(Df1)<-c("pos","Quantile","value","refer_ID1")
    pos <- Df1$pos
    value <- Df1$value
    Quantile <- Df1$Quantile
    p3 <- ggplot(Df1, aes(pos, value, colour = Quantile)) + 
      geom_line() + 
      facet_grid(refer_ID1~.) + 
      annotate("text", x = 0.5, y = -0.6, label = "5'UTR", size = 3) + 
      annotate("text", x = 1.5, y = -0.6, label = "CDS", size = 3) + 
      annotate("text", x = 2.5, y = -0.6, label = "3'UTR", size = 3) + 
      geom_vline(xintercept = 1:2, linetype = "dotted") +
      annotate("rect", xmin = 0.025, xmax = 0.975, ymin = -0.3, ymax = -0.1, alpha = 0.99, colour = "black") +
      annotate("rect", xmin = 1.025, xmax = 1.975, ymin = -0.3, ymax = -0.1, alpha = 0.99, colour = "black") + 
      annotate("rect", xmin = 2.025, xmax = 2.975, ymin = -0.3, ymax = -0.1, alpha = 0.99, colour = "black") +
      theme(axis.title.x =element_text(size=9), axis.title.y=element_text(size=9), 
            title = element_text(size = 9),
            legend.key.height=unit(0.5,'cm'),
            legend.key.width=unit(0.25,'cm'),
            legend.text=element_text(size=9),
            legend.title=element_text(size=9))+
      labs(title = paste("refer_IP and Unified Input Reads Coverage within", condition2),x = "mRNA", y = "Normalized Reads Density") 
    
    colnames(Df2)<-c("pos","Quantile","value","refer_ID2")
    pos <- Df2$pos
    value <- Df2$value
    Quantile <- Df2$Quantile
    p4 <- ggplot(Df2, aes(pos, value, colour = Quantile)) +
      geom_line() + 
      facet_grid(refer_ID2~.) +
      annotate("text", x = 0.5, y = -0.6, label = "5'UTR", size = 3) + 
      annotate("text", x = 1.5, y = -0.6, label = "CDS", size = 3) + 
      annotate("text", x = 2.5, y = -0.6, label = "3'UTR", size = 3) + 
      geom_vline(xintercept = 1:2, linetype = "dotted") + 
      annotate("rect", xmin = 0.025, xmax = 0.975, ymin = -0.3, ymax = -0.1, alpha = 0.99, colour = "black") + 
      annotate("rect", xmin = 1.025, xmax = 1.975, ymin = -0.3, ymax = -0.1, alpha = 0.99, colour = "black") + 
      annotate("rect", xmin = 2.025, xmax = 2.975, ymin = -0.3, ymax = -0.1, alpha = 0.99, colour = "black") + 
      theme(axis.title.x =element_text(size=9), axis.title.y=element_text(size=9),
            title = element_text(size = 9),
            legend.key.height=unit(0.5,'cm'),
            legend.key.width=unit(0.25,'cm'),
            legend.text=element_text(size=9),
            legend.title=element_text(size=9))+
      labs(title = paste("refer_Input and Unified IP Reads Coverage within", condition2),x = "mRNA", y = "Normalized Reads Density")
    
    .multiplot(p1, p2, p3, p4, cols = 2)
  }
}
