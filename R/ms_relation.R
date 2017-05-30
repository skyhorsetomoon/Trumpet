.ms_relation <- function(result, IP_BAM, Input_BAM, contrast_IP_BAM, contrast_Input_BAM, 
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
  ms <- function(group1, group2, group1_name, group2_name, condition) {
    IPa <- data.frame()
    ID <- vector()
    group1 <- as.matrix(group1)
    for (i in seq_len(ncol(group1))) {
      s1 <- .singleBAMreads(group1[, i], se, len, n)
      outa <- .getmeanSD(s1)
      Mean <- outa[[1]]
      SD <- outa[[2]]
      id <- rep(group1_name[i], length(SD))
      ID <- c(ID, id)
      id <- vector()
      IP <- cbind(Mean, SD)
      IPb <- IP
      IPa <- rbind(IPa, IPb)
      IP <- matrix()
    }
    v <- .unified_sample(group2,se,ind,len,n)
    outb <- .getmeanSD(v)
    Mean <- outb[[1]]
    SD <- outb[[2]]
    Input <- cbind(Mean, SD)
    Input <- as.data.frame(Input)
    Id <- rep(group2_name, length(Input$SD))
    ID <- c(ID, Id)
    com <- rbind(IPa, Input)
    com <- cbind(com, ID)
    com <- as.data.frame(com)
    cod <- rep(condition, length(com$ID))
    com <- cbind(com, cod)
    com <- as.data.frame(com)
    names(com) <- c("Mean", "SD", "ID", "Condition")
    return(com)
  }
  if ((length(reference_IP_groupname) == 0) & (length(reference_Input_groupname) == 
                                               0)) {
    Group_IP <- sa[, (seq_len(length(IP_groupname)))]
    Group_Input <- sa[, -(seq_len(length(IP_groupname)))]
    coma <- ms(Group_IP, Group_Input, IP_groupname, "unified_Input", 
               paste(condition1, "condition"))
    comb <- ms(Group_Input, Group_IP, Input_groupname, "unified_IP", 
               paste(condition1, "condition"))
    coma <- as.data.frame(coma)
    Mean <- coma$Mean
    SD <- coma$SD
    ID <- coma$ID
    m_p1 <- ggplot(coma, aes(Mean, SD, colour = ID)) + facet_grid(~Condition) + 
      geom_smooth(aes(group = ID), span = 0.5) + geom_point(alpha = I(1/150), 
                                                            size = 0.2) + theme(title = element_text(size = 12, color = "black")) + 
      labs(title = "Mean-SD relationship within IP group under one condition")
    comb <- as.data.frame(comb)
    Mean <- comb$Mean
    SD <- comb$SD
    ID <- comb$ID
    m_p2 <- ggplot(comb, aes(Mean, SD, colour = ID)) + facet_grid(~Condition) + 
      geom_smooth(aes(group = ID), span = 0.5) + geom_point(alpha = I(1/150), 
                                                            size = 0.2) + theme(title = element_text(size = 12, color = "black")) + 
      labs(title = "Mean-SD relationship within Input group under one condition")
    .multiplot(m_p1, m_p2, cols = 1)
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
    com1 <- ms(group_IP, group_Input, IP_groupname, "unified_Input", 
               paste(condition1, "condition"))
    com2 <- ms(ref_group_IP, ref_group_Input, reference_IP_groupname, 
               "unified_Input", paste(condition2, "condition"))
    coma <- rbind(com1, com2)
    coma <- as.data.frame(coma)
    com3 <- ms(group_Input, group_IP, Input_groupname, "unified_IP", 
               paste(condition1, "condition"))
    com4 <- ms(ref_group_Input, ref_group_IP, reference_Input_groupname, 
               "unified_IP", paste(condition2, "condition"))
    comb <- rbind(com3, com4)
    comb <- as.data.frame(comb)
    Mean <- coma$Mean
    SD <- coma$SD
    ID <- coma$ID
    m_p1 <- ggplot(coma, aes(Mean, SD, colour = ID)) + facet_grid(~Condition) + 
      geom_smooth(aes(group = ID), span = 0.5) + geom_point(alpha = I(1/150), 
                                                            size = 0.2) + theme(title = element_text(size = 12, color = "black")) + 
      labs(title = "Mean-SD relationship within IP group under two condition")
    Mean <- comb$Mean
    SD <- comb$SD
    ID <- comb$ID
    m_p2 <- ggplot(comb, aes(Mean, SD, colour = ID)) + facet_grid(~Condition) + 
      geom_smooth(aes(group = ID), span = 0.5) + geom_point(alpha = I(1/150), 
                                                            size = 0.2) + theme(title = element_text(size = 12, color = "black")) + 
      labs(title = "Mean-SD relationship within Input group under two condition")
    .multiplot(m_p1, m_p2, cols = 2)
  }
}

