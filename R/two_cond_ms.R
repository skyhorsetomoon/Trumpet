.two_cond_ms <- function(result, IP_BAM, Input_BAM, contrast_IP_BAM, contrast_Input_BAM, 
                         condition1, condition2) {
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
  if (is.null(reference_IP_groupname) & is.null(reference_Input_groupname)) {
    print("Must provide two condition bam files using this function")
  }
  s <- result[[1]]
  ind <- unique(s$pos)
  len <- length(ind)
  n <- nrow(s)
  se <- seq(1, n, len)
  sa <- s[, -(1:2)]
  con_ms_f <- function(group1, group2, cond_name1, cond_name2) {
    com_bam <- function(group) {
      v <- vector()
      for (i in seq_len(ncol(group))) {
        m <- vector()
        m <- .trans_readsvector(group[, i], se, len, n)
        v <- rbind(v, m)
      }
      return(v)
    }
    meanSDrel <- function(group) {
      p <- com_bam(group)
      size <- rowSums(p)
      size_factor <- size/exp(mean(log(size)))
      q <- apply(p, 2, function(x, a) x/a, a = size_factor)
      Mean <- apply(q, 2, mean)
      SD <- apply(q, 2, sd)
      com <- cbind(Mean, SD)
      com <- as.data.frame(com)
      z <- which(com$Mean == 0)
      com <- com[-z, ]
      Mean <- log10(com$Mean)
      SD <- log10(com$SD)
      com <- cbind(Mean, SD)
      com <- as.data.frame(com)
      return(com)
    }
    com1 <- meanSDrel(group1)
    com2 <- meanSDrel(group2)
    com <- rbind(com1, com2)
    com <- as.data.frame(com)
    ID <- rep(c(cond_name1, cond_name2), c(length(com1$Mean), length(com2$Mean)))
    com <- cbind(com, ID)
    com <- as.data.frame(com)
    return(com)
  }
  if ((length(reference_IP_groupname) != 0) & (length(reference_Input_groupname) != 
                                               0) & ((length(reference_IP_groupname) + length(reference_Input_groupname)) <= 
                                                     2 | (length(IP_groupname) + length(Input_groupname)) <= 2)) {
    print("The number of samples in each condition should be more than three when using this function")
  }
  if ((length(reference_IP_groupname) != 0) & (length(reference_Input_groupname) != 
                                               0) & ((length(reference_IP_groupname) + length(reference_Input_groupname)) > 
                                                     2)) {
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
    m_com1 <- con_ms_f(group_IP, ref_group_IP, paste("IP group in", 
                                                     condition1, "condition"), paste("IP group under", condition2, 
                                                                                     "condition"))
    m_com2 <- con_ms_f(group_Input, ref_group_Input, paste("Input group in", 
                                                           condition1, "condition"), paste("Input group under", condition2, 
                                                                                           "condition"))
    m_com1 <- as.data.frame(m_com1)
    Mean <- m_com1$Mean
    SD <- m_com1$SD
    ID <- m_com1$ID
    lp1 <- ggplot(m_com1, aes(Mean, SD, colour = ID)) + geom_smooth(aes(group = ID), 
                                                                    span = 0.5) + geom_point(alpha = I(1/200), size = 0.002) + 
      theme(title = element_text(size = 14, color = "black")) + labs(title = paste(" IP group's Mean-SD relationship within two condition"))
    m_com2 <- as.data.frame(m_com2)
    Mean <- m_com2$Mean
    SD <- m_com2$SD
    ID <- m_com2$ID
    lp2 <- ggplot(m_com2, aes(Mean, SD, colour = ID)) + geom_smooth(aes(group = ID), 
                                                                    span = 0.5) + geom_point(alpha = I(1/200), size = 0.002) + 
      theme(title = element_text(size = 14, color = "black")) + labs(title = paste("Input group's Mean-SD relationship within two condition"))
    .multiplot(lp1, lp2, cols = 2)
  }
}
