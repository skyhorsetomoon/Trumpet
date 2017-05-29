.getmeanSD <- function(matr) {
  Mean <- apply(matr, 1, mean)
  SD <- apply(matr, 1, sd)
  Input <- cbind(Mean, SD)
  Input <- as.data.frame(Input)
  z <- which(Input$Mean == 0)
  Input <- Input[-z, ]
  Mean <- log10(Input$Mean)
  SD <- log10(Input$SD)
  out <- list(Mean, SD)
  return(out)
}
