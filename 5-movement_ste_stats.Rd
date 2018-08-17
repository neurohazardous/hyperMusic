library(car)
setwd('/Users/hectorOrozco/Desktop/hM_analysis/hm_nuri/')
mydata = read.csv("movement_ste_horizontal.csv")  # read csv file 
movement.idata = read.csv("movement_levels_r.csv")
# Future Hector: I did most of the data wrangling with excel. 
movement.mat <- as.matrix(mydata) #get ste


# Dyad 1
movement.D1 <- as.matrix(movement.mat[1:8,])
movement.D1.mlm <- lm(movement.D1~1)
movement.D1.aov <- Anova(movement.D1.mlm, idata=movement.idata, idesign=~Direction*Duo,type="III")
summary(movement.D1.aov,multivariate=F)



# Dyad 2
movement.D2 <- as.matrix(movement.mat[9:16,])
movement.D2.mlm <- lm(movement.D2~1)
movement.D2.aov <- Anova(movement.D2.mlm, idata=movement.idata, idesign=~Direction*Duo,type="III")
summary(movement.D2.aov,multivariate=F)

# Both dyads together
movement.mlm <- lm(movement.mat~1)
movement.aov <- Anova(movement.mlm, idata=movement.idata, idesign=~Direction*Duo,type="III")
summary(movement.aov,multivariate=F)

# The way we wanted to originally
movement.1.mean <- colMeans(movement.mat[1:8,])
movement.2.mean <- colMeans(movement.mat[9:16,])
movement.means.mat <- rbind(movement.1.mean, movement.2.mean)
movement.mlm <- lm(movement.means.mat~1)
movement.aov <- Anova(movement.mlm, idata=movement.idata, idesign=~Direction*Duo,type="III")
summary(movement.aov,multivariate=F)


# Post hoc shit
movement.H <- as.matrix(movement.mat[,c(1,3)])
movement.P <- as.matrix(movement.mat[,c(2,4)])

movement.H.mean <- rowMeans(movement.H)
movement.P.mean <- rowMeans(movement.P)
t.test(movement.P.mean, movement.H.mean, paired=TRUE)