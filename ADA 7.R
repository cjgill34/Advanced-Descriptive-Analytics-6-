library(Hmisc) 
library(leaps)
library(MASS)
library(NbClust)
library(dendextend)
library(ggplot2)
library(dplyr)

#Question 15.2 

pharmaceutical_data <- read.csv("/Users/Cassie Gill/OneDrive/SCMA 854/CGill HW7/SCMA 854 HW7/Pharmaceuticals.csv", header=TRUE)
row.names(pharmaceutical_data) <- pharmaceutical_data[,2]
pharmaceutical_data <- pharmaceutical_data[, -c(1,2,12,13,14)]
pharmaceutical_data.norm <- sapply(pharmaceutical_data, scale)

pharmaceutical_data.norm

nc <- NbClust(pharmaceutical_data.norm, distance="euclidean", min.nc=2, max.nc=10, method="average")

table(nc$Best.n[1,])
barplot(table(nc$Best.n[1,]), xlab="Number of Clusters", ylab="Number of criteria", main="Number of clusters chosen by criteria")
d <- dist(pharmaceutical_data.norm)

fit.average <- hclust(d, method="average")

plot(fit.average, hang = -1, cex=0.8, main="average linkage clustering")

clusters <- cutree(fit.average, k=6)

aggregate(pharmaceutical_data.norm, by=list(cluster=clusters), median)
rect.hclust(fit.average, k=6)

# The best number of clusters is 3 based on majority the tree is cut to 6
# because the dissimilarity increases between clusters when we go from 5 to 6

dend <- as.dendrogram(fit.average)
dend <- rotate(dend, 1:21)


dend <- color_branches(dend, k=6)

dend <- hang.dendrogram(dend,hang_height=0.1)

plot(dend, 
     main = "Clustered Pharma data set", 
     horiz =  TRUE,  nodePar = list(cex = .007))


clust.means <- function(x, res.clust, groups)
{
  if(!is.matrix(x))
    x <- as.matrix(x)
  means <- tapply(x, list(rep(cutree(res.clust, groups), ncol(x)),
                          col(x)),
                  mean)
  dimnames(means) <- list(NULL, dimnames(x)[[2]])
  return(as.data.frame(means))
}


pharma_centroids<-clust.means(pharmaceutical_data.norm, fit.average, 6)

y<-apply(as.matrix(pharma_centroids),2,as.double )


plot(c(0), xaxt = 'n', ylab = "", type = "l", ylim = c(min(y), max(y)), xlim = c(0,9))
axis(1, at = c(1:9), labels = names(pharmaceutical_data))
for (i in c(1:6))
  lines(y[i,],  lty = i, lwd = 2,
        col = ifelse(i %in% c(1),"blue",
                     (ifelse(i %in% c(2),"green",
                             (ifelse(i %in% c(3),"red",
                                     (ifelse(i %in% c(4,5),"black","dark grey"))))))))

text(x = 0.5, y = pharma_centroids[, 1], labels = paste("Cluster", c(1:6)))


pharmaceutical_data
pc<-pharmaceutical_data[,1]
row.names(pharmaceutical_data) <- pharmaceutical_data[,1]
pharmaceutical_data <- pharmaceutical_data[, -c(1,2,12,13,14)]
pharmaceutical_data.norm <- sapply(pharmaceutical_data, scale)
set.seed(1)

devAskNewPage(ask=TRUE)

nc <- NbClust(pharmaceutical_data.norm, min.nc=2, 
              max.nc=10, method="kmeans")
table(nc$Best.n[1,])

barplot(table(nc$Best.n[1,]), xlab="Number of Clusters", ylab="Number of criteria", main="Number of clusters chosen by criteria")

# The clusters are reasonable but cluster 5 could be high risk low recovery

wssplot <- function(pharmaceutical_data.norm, nc=10, seed=42) {
  wss <- (nrow(pharmaceutical_data.norm)-1)*sum(apply(pharmaceutical_data.norm, 2, var)) 
  for (i in 2:nc) {
    set.seed(42) 
    wss[i] <- sum(kmeans(pharmaceutical_data.norm, centers=i)$withinss)
  } 
  plot(1:nc, wss, type="b", xlab="Number of clusters", ylab="Within groups sum of squares")
}
wssplot(pharmaceutical_data.norm,nc=10)

fit.km <- kmeans(pharmaceutical_data.norm, 4, nstart=10)
fit.km$size
fit.km$centers
fit.km$withinss

fit.km$centers

fit.km$withinss

fit.km$centers

dist(fit.km$centers)

plot(c(0), xaxt = 'n', ylab = "", type = "l", ylim = c(min(fit.km$centers), max(fit.km$centers)), xlim = c(0, 8))
axis(1, at = c(1:9), labels = names(pharmaceutical_data))
for (i in c(1:4))
  lines(fit.km$centers[i,], lty = i, lwd = 2, 
        col = ifelse(i %in% c(1),"black",
                     (ifelse(i %in% c(2),"blue",
                             (ifelse(i %in% c(3),"green",
                                     (ifelse(i %in% c(4),"red","dark grey"))))))))
text(x = 0.5, y = fit.km$centers[, 1], labels = paste("Cluster", c(1:4))) 


