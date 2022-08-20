library(DT)
library(factoextra)
library(EFAtools)
library(psych)
library(PerformanceAnalytics)
library(readxl)
library(ggplot2)
library(cluster)
library(tidyverse)
library(mclust)
library(dplyr)

################################ PCA ##########################################
data=read.csv("D:/wine-clustering.csv", sep=",", header=TRUE)
datatable(data)

#Uji Multikolinearitas
chart.Correlation(data)
#ada yg lebih dr 0.5, maka dikatakan terdapat multiko

data = subset(data, select = -c(Ash) )

#Uji Asumsi
#KMO Test
KMO(data)

##Asumsi Homogenitas
bartlett.test(data)

##Lakukan pengencekan komponen utama
data.pr <- prcomp(data, center = TRUE, scale = TRUE)
data.pr
summary(data.pr)

screeplot(data.pr, type = "l", npcs = 7, main = "Screeplot")
abline(h = 1, col="red", lty=7)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=7, cex=0.6)
cumpro <- cumsum(data.pr$sdev^2 / sum(data.pr$sdev^2))

#Visualisasi eigen value
plot(cumpro[0:7], xlab = "PC #", ylab = "Amount of explained variance",main = "Cumulative variance plot")
abline(v = 2, col="blue", lty=7)
abline(h = 0.6450, col="blue", lty=7)

#cek PC
round(data.pr$rotation[,1:4],3)
round(data.pr$sdev^2,3) 

library("factoextra")
fviz_pca_var(data.pr, col.var="coord",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE, # Avoid text overlapping
             axes = c(1, 2)) # choose PCs to plot
fviz_pca_ind(data.pr, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             col.ind = "black", 
             palette = "jco", 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Covid") +
  ggtitle("Clustering according to PC") +
  theme(plot.title = element_text(hjust = 0.5))

##data hasil PCA
pca_fix=data.pr$x[,1:3]
pca_fix

new_Data=as.data.frame(pca_fix)
new_Data
bartlett.test(new_Data)

#standarisasi data
datafix <- scale(new_Data) 
datafix
write.csv(datafix, "D:/DataBuatClusterAD.csv", row.names=FALSE)


############################# CLUSTER ######################################
data <- read.csv("D:/DataBuatClusterAD.csv",sep = ",")
summary(data)

#Menentukan k Optimum
library(factoextra)
# Elbow method = optimum pada antara 2 titik yang grafiknya melandai
fviz_nbclust(data, kmeans, method = "wss")
# Silhouette method = optimum pada titik tertinggi (hasil 2 metode ini dibandingkan dan dipilih yang paling sesuai kriteria optimum masing2 metode)
fviz_nbclust(data, kmeans, method = "silhouette")

#Analisis Cluster
#K=3
kmeans_1 <- kmeans(data,3)
databind_1 <- cbind(data, Cluster = kmeans_1$cluster)
head(databind_1)
options(scipen = 9)
centers_1 <- data.frame(kmeans_1$centers)
centers_1$CLUSTER <- 1:nrow(centers_1)
centers_1s <- data.frame(t(centers_1[-3]))
colnames(centers_1s) <- centers_1[,3]
centers_1s
g1 <- fviz_cluster(kmeans_1, geom=c("point","text"), data=data) + ggtitle("k = 3") + theme_light()
plot(g1)
library(dplyr)
data %>%
  group_by(kmeans_1$cluster) %>%
  summarise_all(funs(mean = mean, median = median, sd = sd)) %>%
  as.data.frame()
##########################ICD RATE###################
icdrate = function(Data, nc, c)
{
  n = dim(Data)[1]
  p = dim(Data)[2]
  X = Data[,1:(p-1)]
  Group = Data[,p]
  p = dim(X)[2]
  Mean.X = matrix(ncol = p, nrow = (nc+1))
  for (i in 1:nc)
  {
    for (j in 1:p)
    {
      Mean.X[i,j] = mean(X[which(Group==i),j])
      Mean.X[(nc+1),j] = mean(X[,j])
    }
  }
  SST = matrix(ncol=p, nrow=n)
  for (i in 1:n)
  {
    for (j in 1:p)
    {
      SST[i,j] = (X[i,j] - Mean.X[(nc+1),j])^2
    }
  }
  SST = sum(sum(SST))
  SSE = matrix(ncol=p, nrow=n)
  for (i in 1:n)
  {
    for (j in 1:p)
    {
      for (k in 1:nc)
      {
        if (Group[i]==k)
        {
          SSE[i,j] = (X[i,j] - Mean.X[k,j])^2
        }
      }
    }
  }
  SSE = sum(sum(SSE))
  Rsq = (SST-SSE)/SST
  icdrate = 1-Rsq
  Pseudof = (Rsq/(c-1))/((icdrate)/(nc-c))
  list(Rsq=Rsq, icdrate=icdrate, pseudof=Pseudof)
}

icdrate(databind_1,length(databind_1),3)

#K = 4
kmeans_2 <- kmeans(data,4)
databind_2 <- cbind(data, Cluster = kmeans_2$cluster)
head(databind_2)
centers_2 <- data.frame(kmeans_2$centers)
centers_2$CLUSTER <- 1:nrow(centers_2)
centers_2s <- data.frame(t(centers_2[-18]))
colnames(centers_2s) <- centers_2[,18]
centers_2s
g2 <- fviz_cluster(kmeans_2, geom=c("point","text"), data=data) + ggtitle("k = 4") + theme_light()
plot(g2)
data %>%
  group_by(kmeans_2$cluster) %>%
  summarise_all(funs(mean = mean, median = median, sd = sd)) %>%
  as.data.frame()
##########################ICD RATE###################
icdrate = function(Data, nc, c)
{
  n = dim(Data)[1]
  p = dim(Data)[2]
  X = Data[,1:(p-1)]
  Group = Data[,p]
  p = dim(X)[2]
  Mean.X = matrix(ncol = p, nrow = (nc+1))
  for (i in 1:nc)
  {
    for (j in 1:p)
    {
      Mean.X[i,j] = mean(X[which(Group==i),j])
      Mean.X[(nc+1),j] = mean(X[,j])
    }
  }
  SST = matrix(ncol=p, nrow=n)
  for (i in 1:n)
  {
    for (j in 1:p)
    {
      SST[i,j] = (X[i,j] - Mean.X[(nc+1),j])^2
    }
  }
  SST = sum(sum(SST))
  SSE = matrix(ncol=p, nrow=n)
  for (i in 1:n)
  {
    for (j in 1:p)
    {
      for (k in 1:nc)
      {
        if (Group[i]==k)
        {
          SSE[i,j] = (X[i,j] - Mean.X[k,j])^2
        }
      }
    }
  }
  SSE = sum(sum(SSE))
  Rsq = (SST-SSE)/SST
  icdrate = 1-Rsq
  Pseudof = (Rsq/(c-1))/((icdrate)/(nc-c))
  list(Rsq=Rsq, icdrate=icdrate, pseudof=Pseudof)
}

icdrate(databind_2,length(databind_2),2)

#K = 5
kmeans_3 <- kmeans(data,5)
databind_3 <- cbind(data, Cluster = kmeans_3$cluster)
head(databind_3)
centers_3 <- data.frame(kmeans_3$centers)
centers_3$CLUSTER <- 1:nrow(centers_3)
centers_3s <- data.frame(t(centers_3[-18]))
colnames(centers_3s) <- centers_3[,18]
centers_3s
g3 <- fviz_cluster(kmeans_3, geom=c("point","text"), data=data) + ggtitle("k = 5") + theme_light()
plot(g3)
data %>%
  group_by(kmeans_3$cluster) %>%
  summarise_all(funs(mean = mean, median = median, sd = sd)) %>%
  as.data.frame()
##########################ICD RATE###################
icdrate = function(Data, nc, c)
{
  n = dim(Data)[1]
  p = dim(Data)[2]
  X = Data[,1:(p-1)]
  Group = Data[,p]
  p = dim(X)[2]
  Mean.X = matrix(ncol = p, nrow = (nc+1))
  for (i in 1:nc)
  {
    for (j in 1:p)
    {
      Mean.X[i,j] = mean(X[which(Group==i),j])
      Mean.X[(nc+1),j] = mean(X[,j])
    }
  }
  SST = matrix(ncol=p, nrow=n)
  for (i in 1:n)
  {
    for (j in 1:p)
    {
      SST[i,j] = (X[i,j] - Mean.X[(nc+1),j])^2
    }
  }
  SST = sum(sum(SST))
  SSE = matrix(ncol=p, nrow=n)
  for (i in 1:n)
  {
    for (j in 1:p)
    {
      for (k in 1:nc)
      {
        if (Group[i]==k)
        {
          SSE[i,j] = (X[i,j] - Mean.X[k,j])^2
        }
      }
    }
  }
  SSE = sum(sum(SSE))
  Rsq = (SST-SSE)/SST
  icdrate = 1-Rsq
  Pseudof = (Rsq/(c-1))/((icdrate)/(nc-c))
  list(Rsq=Rsq, icdrate=icdrate, pseudof=Pseudof)
}

icdrate(databind_3,length(databind_3),2)
