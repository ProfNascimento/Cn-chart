library(copula)
require(MSQC);data("water1");data("water2")

#-----------------------------------------------------------#
## Build the bivariate distribution analysis
## (PART I)
#-----------------------------------------------------------#
# DESCRIPTIVE
summary(1/water1[,1]);sd(1/water1[,1])
summary(1/sqrt(water1[,2]));sd(1/sqrt(water1[,2]))

shapiro.test(1/water1[,1])
shapiro.test(1/sqrt(water1[,2]))
mvnormtest::mshapiro.test(t(cbind(1/water1[,1],1/sqrt(water1[,2]))))

cor(1/water1[,1],1/sqrt(water1[,2]),method = "kendal")

# DATA FIT
s1=fitdistrplus::fitdist(1/water1[,1],"norm");plot(s1)
s2=fitdistrplus::fitdist(1/sqrt(water1[,2]),"logis");plot(s2)

# Build the bivariate distribution (Copula)
## (PART I)
my_dist <- mvdc(claytonCopula(param = 2.1567, dim = 2), margins = c("norm","logis"), 
                paramMargins = list(list(mean = s1$estimate[1], sd = s1$estimate[2]), 
                                    list(location = s2$estimate[1], scale = s2$estimate[2])))

# Generate random sample observations from the multivariate distribution
v1 <- rMvdc(5000, my_dist)
# Compute the density
pdf_mvd <- dMvdc(v1, my_dist)
# Compute the CDF
cdf_mvd <- pMvdc(v1, my_dist)

# 3D plain scatterplot of the generated bivariate distribution
par(mfrow = c(1, 2))
# scatterplot3d::scatterplot3d(1/water1[,1],1/sqrt(water1[,2]), pdf_mvd, color="red", main="Density", xlab = "u1", ylab="u2", zlab="pMvdc",pch=".")
# scatterplot3d::scatterplot3d(1/water1[,1],1/sqrt(water1[,2]), cdf_mvd, color="red", main="CDF", xlab = "u1", ylab="u2", zlab="pMvdc",pch=".")
persp(my_dist, dMvdc, xlim = c(0.08, 0.2), ylim=c(0, 5), main = "Density")
contour(my_dist, dMvdc, xlim = c(0.08, 0.2), ylim=c(0, 5), main = "Contour plot")
dev.off()

##
# Plot copula density with histograms
#nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE) 
#par(mar=c(3,3,1,1)) 
contour(my_dist, dMvdc, xlim = c(0.08, 0.2), ylim=c(0, 5), main = "")
abline(v=qnorm(0.975,mean=s1$estimate[1],sd=s1$estimate[2]),col="red",lty=2)
abline(v=qnorm(0.025,mean=s1$estimate[1],sd=s1$estimate[2]),col="red",lty=2)
abline(h=qlogis(0.975,location=2.896, scale=0.445),col="red",lty=2)
abline(h=qlogis(0.025,location=2.896, scale=0.445),col="red",lty=2)
#par(mar=c(0,3,1,1)) 
hist(1/water1[-29,1], breaks=8, xlim=c(0.08,0.2), main="", xaxt = "n", yaxt = "n") 
#par(mar=c(3,0,1,1),crt=45) 
hist(1/sqrt(water1[-29,2]), breaks=6, xlim=c(0, 5), main="", yaxt = "n", xaxt = "n")
#dev.off()

##
contour(my_dist, dMvdc, xlim = c(0.08, 0.2), ylim=c(0, 5), main = "",levels=c(1.1) )
points(1/water1[,1],1/sqrt(water1[,2]),pch=16,cex=1.5,col="black")
points(1/water2[,1],1/sqrt(water2[,2]),pch=16,cex=1.5,col="blue")

# PLOT COPULA DENSITY
plot(dMvdc(as.matrix(cbind(1/c(water1[,1],water2[,1]),1/sqrt(c(water1[,2],water2[,2])))), 
           my_dist),type="b",ylab="Copula Density"); abline(h=1,col="red")
abline(v=30.5,lty=2);text(15,20, "PHASE I");text(45,20, "PHASE II")

# PLOT MARGINALS
plot(dnorm(1/c(water1[,1],water2[,1]),mean = s1$estimate[1], sd = s1$estimate[2]),type="b",ylim=c(0,25), 
     ylab=paste0("Normal(",round(s1$estimate[1],2),",",round(s1$estimate[2],2),") Density" ))
abline(h=pnorm(0.95,mean = s1$estimate[1], sd = s1$estimate[2]),col="red")
abline(v=30.5,lty=2);text(15,25, "PHASE I");text(45,25, "PHASE II")

plot(dlogis(1/sqrt(c(water1[,2],water2[,2])),location = s2$estimate[1], scale = s2$estimate[2]),type="b", ylim=c(0,0.6),
     ylab=paste0("Logistic(",round(s2$estimate[1],2),",",round(s2$estimate[2],2),") Density" ))
abline(h=plogis(0.95,location = s2$estimate[1], scale = s2$estimate[2]),col="red")
abline(v=30.5,lty=2);text(15,0.6, "PHASE I");text(45,0.6, "PHASE II")
#-----------------------------------------------------------#
## Build the bivariate distribution (with conditional mean)
## (PART II)
#-----------------------------------------------------------#
fit_lm=lm(water1[,1]~water1[,3:5])
summary(fit_lm)
# OPT MODEL
fit2_lm=lm(water1[,1]~water1[,4])
summary(fit2_lm)

s11=fitdistrplus::fitdist(residuals(fit2_lm),"norm");plot(s11)
s12=fitdistrplus::fitdist(water1[,2],"logis");plot(s12)

library(VineCopula)
u <- pobs(as.matrix(cbind(residuals(fit2_lm),water1[,2])))[,1]
v <- pobs(as.matrix(cbind(residuals(fit2_lm),water1[,2])))[,2]
selectedCopula <- BiCopSelect(u,v,familyset=NA)
selectedCopula

# Build the bivariate distribution (with conditional mean)
my_dist2 <- mvdc(VC2copula::tawnT1Copula(param = c(selectedCopula$par,selectedCopula$par2)), margins = c("norm","logis"), 
                paramMargins = list(list(mean = s11$estimate[1], sd = s11$estimate[2]), 
                                    list(location = s12$estimate[1], scale = s12$estimate[2])))

# Generate random sample observations from the multivariate distribution
v2 <- rMvdc(5000, my_dist2)
# Compute the density
pdf_mvd2 <- dMvdc(v2, my_dist2)
# Compute the CDF
cdf_mvd2 <- pMvdc(v2, my_dist2)

# 3D plain scatterplot of the generated bivariate distribution
par(mfrow = c(1, 2))
persp(my_dist2, dMvdc, ylim = c(0, 0.3), xlim=c(-1.5, 2), main = "Density")
contour(my_dist2, dMvdc, ylim = c(0, 0.3), xlim=c(-1.5, 2), main = "Contour plot")
dev.off()

res2=water2[,1]-(fit2_lm$coefficients[1]+fit2_lm$coefficients[2]*water2[,4])
plot(dMvdc(as.matrix(cbind(c(residuals(fit2_lm),res2),c(water1[,2],water2[,2]) )), my_dist2),type="b",
     ylab="Copula Density")
abline(h=1,col="red")
abline(v=30.5,lty=2);text(15,8, "PHASE I");text(45,8, "PHASE II")


##
contour(my_dist2, dMvdc, ylim = c(0, 0.27), xlim=c(-1.7, 1.5), main = "",levels=c(1) )
points(residuals(fit2_lm),water1[,2],pch=16,cex=1.5,col="black")
points(res2,water2[,2],pch=16,cex=1.5,col="blue")

# MARGINAL LM
plot(c(residuals(fit2_lm),res2),type="b")
abline(h=qnorm(0.95,mean = s11$estimate[1], sd = s11$estimate[2]),col="red",lty=2)
abline(h=qnorm(0.05,mean = s11$estimate[1], sd = s11$estimate[2]),col="red",lty=2)
abline(v=30.5,lty=2);text(5,1.1, "PHASE I");text(35,1.1, "PHASE II")
car::qqPlot(c(residuals(fit2_lm),res2))
# MARGINAL OXYGEN (ENDOGENOUS VARIABLE)
plot(c(water1[,4],water2[,4]),type="b")

# MARGINAL PO4
plot(dlogis(c(water1[,2],water2[,2]),location = s12$estimate[1], scale = s12$estimate[2]),type="b",ylim=c(0,10),
     ylab=paste0("Logistic(",round(s12$estimate[1],2),",",round(s12$estimate[2],2),") Density" ))
abline(h=plogis(0.95,location = s12$estimate[1], scale = s12$estimate[2]),col="red")
abline(v=30.5,lty=2);text(15,10, "PHASE I");text(45,10, "PHASE II")

