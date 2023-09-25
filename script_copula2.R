library("copula")
library("tidyverse")
library("GGally")
require(MSQC);data("bimetal1");data("bimetal2")

# DESCRIPTIVE
summary(bimetal1[,1]);sd(bimetal1[,1])
summary(bimetal1[,2]);sd(bimetal1[,2])
summary(bimetal1[,3]);sd(bimetal1[,3])

# Skewness -0.5 < x < 0.5 is SYMMETRIC
skewness(bimetal1[,1]);skewness(bimetal1[,2]);skewness(bimetal1[,3])
# Kurtosis close to 3 is MESOKURTIC (data concentrated around the mean)
kurtosis(bimetal1[,1]);kurtosis(bimetal1[,2]);kurtosis(bimetal1[,3])


##---VISUAL DATAMINING---##
cor_fun <- function(data, mapping, method="kendall", ndp=2, sz=8, stars=TRUE){ #ndp is to adjust the number of decimals
  
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  corr <- cor.test(x, y, method=method)
  est <- corr$estimate
  lb.size <- sz 
  
  if(stars){
    stars <- c("***", "**", "*", "")[findInterval(corr$p.value, c(0, 0.001, 0.01, 0.05, 1))]
    lbl <- paste0(round(est, ndp), stars)
  }else{
    lbl <- round(est, ndp)
  }
  
  ggplot(data=data, mapping=mapping) +
    annotate("text", x=mean(x, na.rm=TRUE), y=mean(y, na.rm=TRUE), label=lbl, size=lb.size)+
    theme(panel.grid = element_blank(), panel.background=element_rect(fill="snow1")) 
  
}

colfunc<-colorRampPalette(c("darkblue","cyan","yellow","red"))

my_fn <- function(data, mapping){
  p <- ggplot(data = data, mapping = mapping) + 
    stat_density2d(aes(fill=..density..), geom="tile", contour = FALSE) +
    scale_fill_gradientn(colours = colfunc(100)) + theme_classic()
}

ggpairs(as.data.frame(bimetal1[,1:3]), columns = c(1,2,3),
        lower=list(continuous=my_fn),
        diag=list(continuous=wrap("densityDiag", fill="gray92")), #densityDiag is a function
        upper=list(continuous=cor_fun)) + theme(panel.background=element_rect(fill="white")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, color = "black"),text = element_text(size = 16)) + 
  theme(axis.text.y = element_text(angle = 0, vjust = 1 , color = "black"),text = element_text(size = 16))

##-----------------------------------------------------##

shapiro.test(bimetal1[,1]);shapiro.test(bimetal1[,2]); shapiro.test(bimetal1[,3])
mvnormtest::mshapiro.test(t(cbind(bimetal1[,1:3])))

cor(bimetal1[,1:3],method = "pearson")
cor(bimetal1[,1:3],method = "kendal")

# DATA FIT
T2=qcc::mqcc(bimetal1[,1:3],type="T2",newdata = bimetal2[,1:3],
        add.stats=TRUE,confidence.level=0.95)
summary(T2)

# Building the SPC multivariate w/ Copula
# VERSION 1
s1=fitdistrplus::fitdist(bimetal1[,1],"norm");plot(s1)
s2=fitdistrplus::fitdist(bimetal1[,2],"norm");plot(s2)
s3=fitdistrplus::fitdist(bimetal1[,3],"norm");plot(s3)

trinorm.cop <- normalCopula(0.2, dim = 3, dispstr = "ar1")     # 0.2 (prior rho)
myD <- mvdc(trinorm.cop, margins = c("norm", "norm", "norm"), paramMargins =
              list(list(mean = s1$estimate[1], sd = s1$estimate[2]), 
                   list(mean = s2$estimate[1], sd = s2$estimate[2]),
                   list(mean = s3$estimate[1], sd = s3$estimate[2])))
fitted.copula <- fitMvdc(bimetal1[,1:3], myD, start = c(21,1,40,1,15,1,0.5), 
                   optim.control = list(trace = TRUE, maxit = 2000))
coef(fitted.copula)
coef(fitted.copula)["rho.1"]                                  # estimated rho

loglikMvdc(c(21.0161,0.2964,40.0161,0.1337,15.1921,0.3215,0.605), 
           bimetal1[,1:3], myD)

##
dataFULL=cbind(c(bimetal1[,1],bimetal2[,1]),
               c(bimetal1[,2],bimetal2[,2]),
               c(bimetal1[,3],bimetal2[,3]))

# PLOT COPULA DENSITY
my_dist=mvdc(normalCopula(coef(fitted.copula)["rho.1"], dim = 3, dispstr = "ar1"), 
         margins = c("norm", "norm", "norm"), 
         paramMargins = list(list(mean = s1$estimate[1], sd = s1$estimate[2]), 
                             list(mean = s2$estimate[1], sd = s2$estimate[2]),
                             list(mean = s3$estimate[1], sd = s3$estimate[2])))

TR=quantile(dMvdc(bimetal1[,1:3], my_dist),0.05)  # Copula Tolerance Region (Phase 1)
plot(dMvdc(dataFULL, my_dist),type="b",ylab="Copula Density"); abline(h=TR,col="red")
abline(v=28.5,lty=2);text(15,7.2, "PHASE I");text(43,7.2, "PHASE II")

# PLOT MARGINALS
plot(dataFULL[,1],type="b", ylab="Deflection level")
abline(h=qnorm(0.95,mean = s1$estimate[1], sd = s1$estimate[2]),col="red")
abline(v=28.5,lty=2);text(15,21.75, "PHASE I");text(43,21.75, "PHASE II")

plot(dataFULL[,2],type="b", ylab="Curvature level")
abline(h=qnorm(0.95,mean = s2$estimate[1], sd = s2$estimate[2]),col="red")
abline(v=28.5,lty=2);text(12,40.4, "PHASE I");text(40,40.4, "PHASE II")

plot(dataFULL[,3],type="b", ylab="Resistivity level")
abline(h=qnorm(0.95,mean = s3$estimate[1], sd = s3$estimate[2]),col="red")
abline(v=28.5,lty=2);text(13,16.3, "PHASE I");text(40,16.3, "PHASE II")
