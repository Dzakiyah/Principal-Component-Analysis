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

data=read.csv("D:/cvdfix.csv", sep=",", header=TRUE)
data=data[,3:8]
datatable(data)


#Uji Multikolinearitas
chart.Correlation(data)
#ada yg lebih dr 0.5, maka dikatakan terdapat multiko

#Uji Asumsi
#KMO Test
KMO(data)

##Asumsi Homogenitas
bartlett.test(data)

##Lakukan pengencekan komponen utama
data.pr <- prcomp(data, center = TRUE, scale = TRUE)
data.pr
summary(data.pr)

screeplot(data.pr, type = "l", npcs = 5, main = "Screeplot")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)
cumpro <- cumsum(data.pr$sdev^2 / sum(data.pr$sdev^2))

#Visualisasi eigen value
plot(cumpro[0:5], xlab = "PC #", ylab = "Amount of explained variance",main = "Cumulative variance plot")
abline(v = 2, col="blue", lty=5)
abline(h = 0.96719, col="blue", lty=5)

#cek PC
round(data.pr$rotation[,1:2],2)
round(data.pr$sdev^2,2) 

library("factoextra")
win.graph()
fviz_pca_var(data.pr, col.var="coord",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE, # Avoid text overlapping
             axes = c(1, 2)) # choose PCs to plot
win.graph()
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
pca_fix=data.pr$x[,1:2]
pca_fix

new_Data=as.data.frame(pca_fix)
new_Data
bartlett.test(new_Data)

#standarisasi data
datafix <- scale(new_Data) 
datafix
write.csv(datafix, "D:/DataBuatClusterAD.csv", row.names=FALSE)

km_fit1 = kmeans(datafix,centers = 2,iter.max = 300 )
win.graph()
fviz_cluster(km_fit1, data = datafix)
new_data <- cbind(data,data.frame(km_fit1$cluster))
new_data$km_fit1.cluster <- as.factor(new_data$km_fit1.cluster)

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

icdrate(new_data,length(new_data),2)

new_data
library(writexl)
write_xlsx(new_data,"D:/2ClusterPCA.xlsx")

##PCA tahap 2
data1=data[4:152,]
##Lakukan pengencekan komponen utama
data1.pr <- prcomp(data1, center = TRUE, scale = TRUE)
data1.pr
summary(data1.pr)

screeplot(data1.pr, type = "l", npcs = 5, main = "Screeplot")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)
cumpro1 <- cumsum(data1.pr$sdev^2 / sum(data1.pr$sdev^2))

#Visualisasi eigen value
plot(cumpro1[0:5], xlab = "PC #", ylab = "Amount of explained variance",main = "Cumulative variance plot")
abline(v = 2, col="blue", lty=5)
abline(h = 0.8962, col="blue", lty=5)

#cek PC
round(data1.pr$rotation[,1:2],2)
round(data1.pr$sdev^2,2) 

library("factoextra")
win.graph()
fviz_pca_var(data1.pr, col.var="coord",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE, # Avoid text overlapping
             axes = c(1, 2)) # choose PCs to plot
win.graph()
fviz_pca_ind(data1.pr, geom.ind = "point", pointshape = 21, 
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
pca_fix1=data1.pr$x[,1:2]
pca_fix1

new_Data1=as.data.frame(pca_fix1)
new_Data1

#standarisasi data
datafix1 <- scale(new_Data1) 
datafix1

km_fit2 = kmeans(datafix1,centers = 2,iter.max = 300 )
win.graph()
fviz_cluster(km_fit2, data = datafix1)
new_data1 <- cbind(data1,data.frame(km_fit2$cluster))
new_data1$km_fit2.cluster <- as.factor(new_data1$km_fit2.cluster)

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

icdrate(new_data1,length(new_data1),2)
new_data1
write_xlsx(new_data1,"D:/2ClusterPCA1.xlsx")

##PCA tahap 3
data2=new_data1[new_data1$km_fit2.cluster=="2",]
data2=data2[,-6]
data2

##Lakukan pengencekan komponen utama
data2.pr <- prcomp(data2, center = TRUE, scale = TRUE)
data2.pr
summary(data2.pr)

screeplot(data2.pr, type = "l", npcs = 5, main = "Screeplot")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)
cumpro2 <- cumsum(data2.pr$sdev^2 / sum(data2.pr$sdev^2))

#Visualisasi eigen value
plot(cumpro2[0:5], xlab = "PC #", ylab = "Amount of explained variance",main = "Cumulative variance plot")
abline(v = 2, col="blue", lty=5)
abline(h = 0.8758, col="blue", lty=5)

#cek PC
round(data2.pr$rotation[,1:2],2)
round(data2.pr$sdev^2,2) 

library("factoextra")
win.graph()
fviz_pca_var(data2.pr, col.var="coord",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE, # Avoid text overlapping
             axes = c(1, 2)) # choose PCs to plot
win.graph()
fviz_pca_ind(data2.pr, geom.ind = "point", pointshape = 21, 
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
pca_fix2=data2.pr$x[,1:2]
pca_fix2

new_Data2=as.data.frame(pca_fix2)
new_Data2

#standarisasi data
datafix2 <- scale(new_Data2) 
datafix2

km_fit3 = kmeans(datafix2,centers = 2,iter.max = 300 )
win.graph()
fviz_cluster(km_fit3, data = datafix2)
new_data2 <- cbind(data2,data.frame(km_fit3$cluster))
new_data2$km_fit3.cluster <- as.factor(new_data2$km_fit3.cluster)

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

icdrate(new_data2,length(new_data2),2)
new_data2
write_xlsx(new_data2,"D:/2ClusterPCA2.xlsx")

##PCA tahap 4
data3=new_data2[new_data2$km_fit3.cluster=="1",]
data3=data3[,-6]
data3

##Lakukan pengencekan komponen utama
data3.pr <- prcomp(data3, center = TRUE, scale = TRUE)
data3.pr
summary(data3.pr)

screeplot(data3.pr, type = "l", npcs = 5, main = "Screeplot")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)
cumpro3 <- cumsum(data3.pr$sdev^2 / sum(data3.pr$sdev^2))

#Visualisasi eigen value
plot(cumpro3[0:5], xlab = "PC #", ylab = "Amount of explained variance",main = "Cumulative variance plot")
abline(v = 2, col="blue", lty=5)
abline(h = 0.8086, col="blue", lty=5)

#cek PC
round(data3.pr$rotation[,1:2],2)
round(data3.pr$sdev^2,2) 

library("factoextra")
win.graph()
fviz_pca_var(data3.pr, col.var="coord",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE, # Avoid text overlapping
             axes = c(1, 2)) # choose PCs to plot
win.graph()
fviz_pca_ind(data3.pr, geom.ind = "point", pointshape = 21, 
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
pca_fix3=data3.pr$x[,1:2]
pca_fix3

new_Data3=as.data.frame(pca_fix3)
new_Data3

#standarisasi data
datafix3 <- scale(new_Data3) 
datafix3

km_fit4 = kmeans(datafix3,centers = 2,iter.max = 300 )
win.graph()
fviz_cluster(km_fit4, data = datafix3)
new_data3 <- cbind(data3,data.frame(km_fit4$cluster))
new_data3$km_fit4.cluster <- as.factor(new_data3$km_fit4.cluster)

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

icdrate(new_data3,length(new_data3),4)
new_data3
write_xlsx(new_data3,"D:/2ClusterPCA3.xlsx")
