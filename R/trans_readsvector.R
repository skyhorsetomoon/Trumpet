.trans_readsvector <- function(bam, se, len, n) {
  w <- vector(mode = "numeric", length = 0)
  for (i in seq_len(len)) {
    se <- seq(i, n, len)
    w <- cbind(w, bam[se])
  }
  w <- as.matrix(w)
  w <- sapply(t(w), unlist)
  return(w)
}
