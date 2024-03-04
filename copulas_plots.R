## Filename: copulas_plots.R
library(copula)
library(dplyr)
library(plotly)
source("Arc_copulas.R")
tau <- 0.7
bi_Clayton <- Arc_copula(p = 2, tau = tau)$clayton.cop
bi_Frank <- Arc_copula(p = 2, tau = tau)$frank.cop
bi_Gumbel <-  Arc_copula(p = 2, tau = tau)$gumbel.cop
bi_Normal <- Arc_copula(p = 2, tau = tau)$normal.cop
set.seed(1234);
nx <- 200
X.frank <- rCopula(nx, copula = bi_Frank)
X.gumbel <- rCopula(nx, copula = bi_Gumber)
X.clayton <- rCopula(nx, copula = bi_Clayton)
X.normal <- rCopula(nx, copula = bi_Normal)

par(mfrow = c(2, 2), cex.axis = 1.2,
    cex.lab = 1.2, mar = c(4, 4, 2,1), cex.main = 1.2
)
plot(X.frank, main = "Frank copula", xlab = "", ylab = "")
plot(X.gumbel, main = "Gumbel copula", xlab = "", ylab = "")
plot(X.clayton, main = "Clayton copula", xlab = "", ylab = "")
plot(X.normal, main = "Normal copula", xlab = "", ylab = "")

#Frank copulas
bi_Frank_cop <- copula::mvdc(bi_Frank, margins = c("unif", "unif"), 
                             paramMargins=list(
                               list(min = 0, max = 1),list(min = 0, max = 1)
                             )
)
x1 <- x2 <- seq(0, 1, by = 0.005)
v<-c()
for (i in x1){
  for (j in x2){
    v<-c(v,copula::dMvdc(c(i,j), bi_Frank_cop) )
  }
}
f<-t(matrix(v,nrow=100,byrow=TRUE))
mycolor <- matrix(c(0, "#FFFFFF", 
                    1, "#2171B5"), nrow = 2, byrow = T)
plotly::plot_ly(x=x1,y=x2,z=f,type = "contour",colorscale = mycolor, 
                reversescale = F, showscale = F)%>%
  plotly::layout(xaxis=list(title="x1"),
                 yaxis=list(title="x2"),
                 title=paste("Frank copulas - tau = 0.7")
  )
############
# Gumbel copulas
bi_Gumbel_cop <- copula::mvdc(bi_Gumbel, margins = c("unif", "unif"), 
                              paramMargins=list(
                                list(min = 0, max = 1) ,list(min = 0, max = 1))
)
x1 <- x2 <- seq(-3, 3, length= 100)
v<-c()
for (i in x1){
  for (j in x2){
    v<-c(v,copula::dMvdc(c(i,j), bi_Gumbel_cop))
  }
}
f<-t(matrix(v,nrow=100,byrow=TRUE))
plotly::plot_ly(x=x1,y=x2,z=f,type = "contour",colorscale = mycolor, 
                reversescale = F, showscale = F)%>%
  plotly::layout(xaxis=list(title="x1"),
                 yaxis=list(title="x2"),
                 title=paste("Gumbel copulas - tau = 0.7"))
############
# Clayton copulas
bi_Clayton_cop <- copula::mvdc(bi_Clayton, margins = c("unif", "unif"), 
                               paramMargins=list(
                                 list(min = 0, max = 1) ,list(min = 0, max = 1)
                               )
)
x1 <- x2 <- seq(0, 1, by = 0.01)
v<-c()
for (i in x1){
  for (j in x2){
    print(c(x1, x2))
    v<-c(v,copula::dMvdc(c(i,j), bi_Clayton_cop))
  }
}
f<-t(matrix(v,nrow=100,byrow=TRUE))

plotly::plot_ly(x=x1,y=x2,z=f,type = "contour",colorscale = mycolor, 
                reversescale = F, showscale = F)%>%
  plotly::layout(xaxis=list(title="x1"),
                 yaxis=list(title="x2"),
                 title=paste("Clayton copulas - tau = 0.7")
  )
############
# Normal copulas
bi_Normal_cop <- copula::mvdc(bi_Normal, margins = c("unif", "unif"), 
                              paramMargins=list(
                                list(min = 0, max = 1) ,list(min = 0, max = 1)
                              )
)
x1 <- x2 <- seq(0, 1, by = 0.01)
v<-c()
for (i in x1){
  for (j in x2){
    v<-c(v,copula::dMvdc(c(i,j), bi_Normal_cop))
  }
}
f<-t(matrix(v,nrow=100,byrow=TRUE))

plotly::plot_ly(x=x1,y=x2,z=f,type = "contour", colorscale = mycolor,
                reversescale = F, showscale = F)%>%
  plotly::layout(xaxis=list(title="x1"),
                 yaxis=list(title="x2"),title=paste("Normal copulas - tau = 0.7"))
###########################################################
