#R script for correlation plot
#Francesca Bottacini
#f.bottacini@umail.ucc.ie


library(corrplot)
library(RColorBrewer)

## Read input matrix with 1row and 1col as headers
matrix <- read.delim("matrix.txt", row.names=1)

tM <- t(matrix)

#Calculate significance

cor.mtest <- function(tM, ...) {
  tM <- as.matrix(tM)
  n <- ncol(tM)
  p.tM<- matrix(NA, n, n)
  diag(p.tM) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(tM[, i], tM[, j], ...)
      p.tM[i, j] <- p.tM[j, i] <- tmp$p.value
    }
  }
  colnames(p.tM) <- rownames(p.tM) <- colnames(tM)
  p.tM
}
# matrix of the p-value of the correlation
p.tM <- cor.mtest(tM)
head(p.tM[, 1:5])

#Corrplot

tM <- cor(tM)
tM[is.na(x = tM)] <- 0
#tM2 <- na.omit(tM)

#corrplot(tM, type="upper", method = "circle")

corrplot(tM, type="upper", method = "circle", tl.col = "black", cl.ratio = 0.1, cl.align = "r", col = brewer.pal(n = 8, name = "PuOr"))

#corrplot(tM, type="upper", method = "circle", p.tM = p.tM, sig.level = 0.01, insig = "blank", na.label = "square", na.label.col = "grey")
#corrplot(tM, type="upper", method = "circle", p.tM = p.tM, sig.level = 0.01, insig = "blank", na.label = "NA", na.label.col = "grey")
#corrplot.mixed(cor(tM), lower="circle", upper="color", tl.pos="lt", diag="n", order="hclust", hclust.method="complete", lab)


