.hc <- function(group_bam, group_name, se, len, n) {
  m <- vector()
  v <- vector()
  group_bam <- as.matrix(group_bam)
  for (i in seq_len(ncol(group_bam))) {
    
    m <- .trans_readsvector(group_bam[, i], se, len, n)
    v <- rbind(m, v)
  }
  rownames(v) <- group_name
  v <- as.data.frame(v)
  size <- rowSums(v)
  size_factor <- size/exp(mean(log(size)))
  v <- apply(v, 2, function(x, a) x/a, a = size_factor)
  return(v)
}

.hcluster <- function(result, IP_BAM, Input_BAM, contrast_IP_BAM, contrast_Input_BAM) {
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
                                               0) & ((length(IP_groupname) + length(Input_groupname)) < 3)) {
    print("The number of IP_BAM or Input_BAM must be more than three.")
  }
  if ((length(reference_IP_groupname) == 0) & (length(reference_Input_groupname) == 
                                               0) & ((length(IP_groupname) + length(Input_groupname)) > 2)) {
    Group_IP <- sa[, (seq_len(length(IP_groupname)))]
    Group_IP <- as.matrix(Group_IP)
    Group_Input <- sa[, -(seq_len(length(IP_groupname)))]
    Group_Input <- as.matrix(Group_Input)
    if (length(IP_groupname) > 2 & length(Input_groupname) > 2) {
      hct_IP <- .hc(Group_IP, IP_groupname, se, len, n)
      ##PCA IP plot
      phc <- prcomp(hct_IP)
      pcd <- phc$x
      main_pca <- pcd[,1:2]
      main_pca <- as.data.frame(main_pca)
      pc1 <- main_pca$PC1
      pc2 <- main_pca$PC2
      
      y_fold <- max(pc2)*1
      x_fold <- max(pc1)*1
      if((max(pc1)>max(pc2))){
        y_fold <- max(pc2)*2.5
      }
      
      if((max(pc1)<max(pc2))){
        
        x_fold <- max(pc1)*2.5
      }
      
      p_IP <- ggplot(main_pca, aes(PC1, PC2)) +geom_point(size=4)+ geom_text(label=rownames(main_pca),colour="red",size=4,vjust = 0, hjust = 2, angle = 270)+
        xlim(min(pc1),x_fold )+
        ylim(min(pc2), y_fold)
      
      
      hcA <- hclust(dist(hct_IP, method = "euclidean"), method = "complete")
      ## PCA Input plot
      hct_Input <- .hc(Group_Input, Input_groupname, se, len, n)
      phc <- prcomp(hct_Input)
      pcd <- phc$x
      main_pca <- pcd[,1:2]
      main_pca <- as.data.frame(main_pca)
      pc1 <- main_pca$PC1
      pc2 <- main_pca$PC2
      
      y_fold <- max(pc2)*1
      x_fold <- max(pc1)*1
      if((max(pc1)>max(pc2))){
        y_fold <- max(pc2)*2.5
      }
      
      if((max(pc1)<max(pc2))){
        
        x_fold <- max(pc1)*2.5
      }
      
      p_Input <- ggplot(main_pca, aes(PC1, PC2)) +geom_point(size=4)+ geom_text(label=rownames(main_pca),colour="red",size=4,vjust = 0, hjust = 2, angle = 270)+
        xlim(min(pc1),x_fold )+
        ylim(min(pc2), y_fold)
      
      hcB <- hclust(dist(hct_Input, method = "euclidean"), method = "complete")
      nf <- layout(matrix(c(1,2),1,2,byrow = TRUE), c(1,1), c(1,1), TRUE)
      par(mar = c(5,5,1,1))
      plot(hcA, hang = -1, main = "IP group cluster")
      par(mar = c(5,5,1,1))
      plot(hcB, hang = -1, main = "Input group cluster")
      
      suppressWarnings(.multiplot(p_IP, p_Input,cols = 2))
    }
    if (length(IP_groupname) > 2 & length(Input_groupname) <= 2) {
      hct_IP <- .hc(Group_IP, IP_groupname, se, len, n)
      hcA <- hclust(dist(hct_IP), method = "complete")
      print("The number of Input_BAM are less three, so you can't get the Input group cluster")
      plot(hcA, hang = -1, main = "IP group cluster")
      ## PCA IP
      phc <- prcomp(hct_IP)
      pcd <- phc$x
      main_pca <- pcd[,1:2]
      main_pca <- as.data.frame(main_pca)
      pc1 <- main_pca$PC1
      pc2 <- main_pca$PC2
      
      y_fold <- max(pc2)*1
      x_fold <- max(pc1)*1
      if((max(pc1)>max(pc2))){
        y_fold <- max(pc2)*2.5
      }
      
      if((max(pc1)<max(pc2))){
        
        x_fold <- max(pc1)*2.5
      }
      
      p_IP <- ggplot(main_pca, aes(PC1, PC2)) +geom_point(size=3)+ geom_text(label=rownames(main_pca),colour="red",size=2.5,vjust = 0.3, hjust = 1.3, angle = 270)+
        xlim(min(pc1),x_fold )+
        ylim(min(pc2), y_fold)
      
      
      
      suppressWarnings(.multiplot(p_IP, cols = 1))
      
      
      
    }
    if (length(IP_groupname) <= 2 & length(Input_groupname) > 2) {
      hct_Input <- .hc(Group_Input, Input_groupname, se, len, n)
      hcB <- hclust(dist(hct_Input), method = "complete")
      print("The number of IP_BAM are less three, so you can't get the IP group cluster")
      plot(hcB, hang = -1, main = "Input group cluster")
      ## PCA Input
      phc <- prcomp(hct_Input)
      pcd <- phc$x
      main_pca <- pcd[,1:2]
      main_pca <- as.data.frame(main_pca)
      pc1 <- main_pca$PC1
      pc2 <- main_pca$PC2
      
      y_fold <- max(pc2)*1
      x_fold <- max(pc1)*1
      if((max(pc1)>max(pc2))){
        y_fold <- max(pc2)*2.5
      }
      
      if((max(pc1)<max(pc2))){
        
        x_fold <- max(pc1)*2.5
      }
      
      p_Input <- ggplot(main_pca, aes(PC1, PC2)) +geom_point(size=3)+ geom_text(label=rownames(main_pca),colour="red",size=2.5,vjust = 0.3, hjust = 1.3, angle = 270)+
        xlim(min(pc1),x_fold )+
        ylim(min(pc2), y_fold)
      
      
      suppressWarnings(.multiplot( p_Input, cols = 1))
      
      
      
    }
  } else if ((length(reference_IP_groupname) != 0) & (length(reference_Input_groupname) != 
                                                      0)) {
    group_IP <- sa[, (seq_len(length(IP_groupname)))]
    group_IP <- as.matrix(group_IP)
    Group_Input <- sa[, -(seq_len(length(IP_groupname)))]
    Group_Input <- as.matrix(Group_Input)
    group_Input <- Group_Input[, -((length(Input_groupname) + 1):ncol(Group_Input))]
    group_Input <- as.matrix(group_Input)
    ref_group <- Group_Input[, -(seq_len(length(Input_groupname)))]
    ref_group_IP <- ref_group[, seq_len(length(reference_IP_groupname))]
    ref_group_IP <- as.matrix(ref_group_IP)
    ref_group_Input <- ref_group[, -(seq_len(length(reference_IP_groupname)))]
    ref_group_Input <- as.matrix(ref_group_Input)
    IP_group <- cbind(group_IP, ref_group_IP)
    group_IPname <- c(IP_groupname, reference_IP_groupname)
    Input_group <- cbind(group_Input, ref_group_Input)
    group_Inputname <- c(Input_groupname, reference_Input_groupname)
    if (((length(IP_groupname) + length(reference_IP_groupname)) < 
         3) & ((length(Input_groupname) + length(reference_Input_groupname)) < 
               3)) {
      print("The number of IP group samples or Input group samples must be more than three.")
      
    }
    
    if (((length(IP_groupname) + length(reference_IP_groupname)) > 
         2) & ((length(Input_groupname) + length(reference_Input_groupname)) < 
               3)) {
      hctA <- .hc(IP_group, group_IPname, se, len, n)
      hcA <- hclust(dist(hctA), method = "complete")
      p1 <- plot(hcA, hang = -1, main = "IP group cluster")
      print("The total number of Input_BAM and contrast_Input_BAM are less three, so you can't get the Input group cluster")
      ## PCA
      phc <- prcomp(hctA)
      pcd <- phc$x
      main_pca <- pcd[,1:2]
      main_pca <- as.data.frame(main_pca)
      pc1 <- main_pca$PC1
      pc2 <- main_pca$PC2
      
      y_fold <- max(pc2)*1
      x_fold <- max(pc1)*1
      if((max(pc1)>max(pc2))){
        y_fold <- max(pc2)*2.5
      }
      
      if((max(pc1)<max(pc2))){
        
        x_fold <- max(pc1)*2.5
      }
      
      p_IP <- ggplot(main_pca, aes(PC1, PC2)) +geom_point(size=3)+ geom_text(label=rownames(main_pca),colour="red",size=2.5,vjust = 0.3, hjust = 1.3, angle = 270)+
        xlim(min(pc1),x_fold )+
        ylim(min(pc2), y_fold)
      
      
      
      suppressWarnings(.multiplot(p_IP, cols = 1))
    }
    if (((length(IP_groupname) + length(reference_IP_groupname)) < 
         3) & ((length(Input_groupname) + length(reference_Input_groupname)) > 
               2)) {
      hctB <- .hc(Input_group, group_Inputname, se, len, n)
      hcB <- hclust(dist(hctB), method = "complete")
      plot(hcB, hang = -1, main = "Input group cluster")
      print("The total number of IP_BAM and contrast_IP_BAM are less three, so you can't get the IP group cluster")
      ## PCA Input
      phc <- prcomp(hctB)
      pcd <- phc$x
      main_pca <- pcd[,1:2]
      main_pca <- as.data.frame(main_pca)
      pc1 <- main_pca$PC1
      pc2 <- main_pca$PC2
      
      y_fold <- max(pc2)*1
      x_fold <- max(pc1)*1
      if((max(pc1)>max(pc2))){
        y_fold <- max(pc2)*2.5
      }
      
      if((max(pc1)<max(pc2))){
        
        x_fold <- max(pc1)*2.5
      }
      
      p_Input <- ggplot(main_pca, aes(PC1, PC2)) +geom_point(size=3)+ geom_text(label=rownames(main_pca),colour="red",size=2.5,vjust = 0.3, hjust = 1.3, angle = 270)+
        xlim(min(pc1),x_fold )+
        ylim(min(pc2), y_fold)
      
      
      suppressWarnings(.multiplot( p_Input, cols = 1))
      
    }
    if (((length(IP_groupname) + length(reference_IP_groupname)) > 
         2) & ((length(Input_groupname) + length(reference_Input_groupname)) > 
               2)) {
      hctA <- .hc(IP_group, group_IPname, se, len, n)
      hcA <- hclust(dist(hctA, method = "euclidean"), method = "complete")
      ##PCA plot
      phc <- prcomp(hctA)
      pcd <- phc$x
      main_pca <- pcd[,1:2]
      main_pca <- as.data.frame(main_pca)
      pc1 <- main_pca$PC1
      pc2 <- main_pca$PC2
      y_fold <- max(pc2)*1
      x_fold <- max(pc1)*1
      if((max(pc1)>max(pc2))){
        y_fold <- max(pc2)*2.5
      }
      
      if((max(pc1)<max(pc2))){
        
        x_fold <- max(pc1)*2.5
      }
      
      
      p_IP <- ggplot(main_pca, aes(PC1, PC2)) +geom_point(size=3)+ geom_text(label=rownames(main_pca),colour="red",size=2.5,vjust = 0.3, hjust = 1.3, angle = 270)+
        xlim(min(pc1),x_fold )+
        ylim(min(pc2), y_fold)
      
      
      
      hctB <- .hc(Input_group, group_Inputname, se, len, n)
      hcB <- hclust(dist(hctB, method = "euclidean"), method = "complete")
      ##PCA plot
      phc <- prcomp(hctB)
      pcd <- phc$x
      main_pca <- pcd[,1:2]
      main_pca <- as.data.frame(main_pca)
      pc1 <- main_pca$PC1
      pc2 <- main_pca$PC2
      
      y_fold <- max(pc2)*1
      x_fold <- max(pc1)*1
      if((max(pc1)>max(pc2))){
        y_fold <- max(pc2)*2.5
      }
      
      if((max(pc1)<max(pc2))){
        
        x_fold <- max(pc1)*2.5
      }
      
      p_Input <- ggplot(main_pca, aes(PC1, PC2)) +geom_point(size=3)+ geom_text(label=rownames(main_pca),colour="red",size=2.5,vjust = 0.3, hjust = 1.3, angle = 270)+
        xlim(min(pc1),x_fold )+
        ylim(min(pc2), y_fold)
      
      
      nf <- layout(matrix(c(1,2),1,2,byrow = TRUE), c(1,1), c(1,1), TRUE)
      par(mar = c(5,5,1,1))
      plot(hcA, hang = -1, main = "IP group cluster")
      par(mar = c(5,5,1,1))
      plot(hcB, hang = -1, main = "Input group cluster")
      suppressWarnings(.multiplot(p_IP, p_Input,cols = 2))
      
      
    }
    
  }
}
