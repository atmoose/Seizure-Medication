setwd("C:/Users/icamo/OneDrive/Documents/Clemson/Fall 2019/MATH 8850/mid-term project")
set.seed(123)

seizure2 <- read.csv("seizure2.csv")

seizure2$id <- 1:dim(seizure2)[1]

library(tidyr)
library(glmmML)
library(lme4)


# Format the dataset into long
long <- seizure2 %>%
  gather('base', 'y1', 'y2', 'y3', 'y4', key = "period", value = "count")
map = setNames(0:4, c('base','y1','y2','y3','y4'))
long[,4] <- map[unlist(long[,4])]
# add in length and cage
idinds<- unlist( lapply( unique(long$id), FUN = l<- function(i) which(long$id == i)[1]))
long <- long[order(long$id),]
long[,c(6,7)] <- c(rep(c(8,2,2,2,2), 57), long$age - mean(long$age[idinds]))
names(long)[c(6,7)] <- c("length","cage")

# Order the data
col.order <- c("id", "age", "cage", "trt", "period", "length", "count")
long <- long[,col.order]
long$fperiod <- factor(long$period)

head(long)


# Make a matrix where there are indicators for the periods & trt
# Removing the baseline b/c it is absorbed into the intercept
new_x<- (model.matrix(~trt * fperiod-1, data=long)  )[,-c(2)]

head(new_x)

nlong <- data.frame(id=long$id, age=long$age, cage=long$cage, count=long$count, new_x,length=long$length)
head(nlong)
AQuad20 <- glmer(count ~ fperiod1 + fperiod2 + fperiod3 + fperiod4 + trt + 
                   trt.fperiod1 + trt.fperiod2 + trt.fperiod3 + trt.fperiod4 + (1|id), 
                 family="poisson", data= nlong, offset=(I(log(length))), nAGQ = 20)


# Bootstrap variance
iter = 500
boot.beta <- matrix(NA, ncol = 10, nrow = iter)
# Ignoring the individuals
for(i in 1:iter){
  samp <- sample(1:dim(nlong)[1], dim(nlong)[1], replace = TRUE)
  model <- glmer(count ~ fperiod1 + fperiod2 + fperiod3 + fperiod4 + trt + 
                        trt.fperiod1 + trt.fperiod2 + trt.fperiod3 + trt.fperiod4 + (1|id),
                        family="poisson", data= nlong[samp,], offset=(I(log(length))), nAGQ = 20)
  boot.beta[i,] <- summary(model)$coefficients[,1]
  if(i%%10 == 0){
    print(i)
  }
}



hat_beta<- summary(AQuad20)$coefficients[,1]

V_boot <- var(boot.beta)


#########################
# Defining our groups
#########################

# baseline of individuals who were NOT treated
x_base_p<- c(1,0,0,0,0,0,0,0,0,0)

# baseline of individuals who were treated
x_base_t<- c(1,0,0,0,0,1,0,0,0,0)

#untreated at time 1
x_p_1 <- c(1,1,0,0,0,0,0,0,0,0)
#untreated at time 2
x_p_2 <- c(1,0,1,0,0,0,0,0,0,0)
#untreated at time 3
x_p_3 <- c(1,0,0,1,0,0,0,0,0,0)
#untreated at time 4
x_p_4 <- c(1,0,0,0,1,0,0,0,0,0)

#treated at time 1
x_t_1 <- c(1,1,0,0,0,1,1,0,0,0)
#treated at time 2
x_t_2 <- c(1,0,1,0,0,1,0,1,0,0)
#treated at time 3
x_t_3 <- c(1,0,0,1,0,1,0,0,1,0)
#treated at time 4
x_t_4 <- c(1,0,0,0,1,1,0,0,0,1)

# Matrix of above groups
x <- t(matrix( t(c(x_base_p - x_base_t, 
                   x_base_p - x_p_1, x_p_1 - x_p_2, x_p_2 - x_p_3, x_p_3 - x_p_4, 
                   x_base_t - x_t_4)),
               nrow = 10))

# Function that calculates all 12 CI's and outputs a matrix of them
CI_Boot<- function(x){ 
  CI_Matrix <- matrix(0,length(x[ ,1]),3)
  for(i in 1:length(x[ ,1])){
    se<- sqrt(x[i, ] %*% V_boot %*% x[i, ])
    estimate<- exp(x[i, ] %*% hat_beta) 
    CI<- exp( c(x[i, ]%*% hat_beta - 1.96 * se , x[i, ] %*%hat_beta + 1.96 * se ) ) 
    CI_Matrix[i,1] <- estimate
    CI_Matrix[i,2:3] <- CI
  }
  colnames(CI_Matrix) <- c("Mean Estimate", "95% LCL", "95% UCL")
  rownames(CI_Matrix) <- c("Baseline Trt vs Plac", "Baseline Plac - Plac 1", 
                           "Plac 1 - Plac 2", "Plac 2 - Plac 3",
                           "Plac 3 - Plac 4", "Baseline Trt - Trt 4")
  return(CI_Matrix)
}
CI_Boot(x)





# AVG number of seizures
par(mfrow=c(1,2))
#first look at those in placebo group 
seizure.plac <- subset(long , trt == 0)
plot(	c(mean(seizure.plac$count[seizure.plac$fperiod == 0]), sum(mean(seizure.plac$count[seizure.plac$fperiod == 1]), mean(seizure.plac$count[seizure.plac$fperiod == 2]),
        mean(seizure.plac$count[seizure.plac$fperiod == 3]), mean(seizure.plac$count[seizure.plac$fperiod == 4]))) ~ c(0:1),
      xaxt="n", ylab="Avg Number of Seizures", xlab="Time Period", type="l", 
      ylim=c(0, 50), main="Placebo")  
axis(side=1, at= 0:1,labels =c("pre*", "post*") )  

#then the progabide group 
seizure.trt <- subset(long, trt == 1)
plot(c(mean(seizure.trt$count[seizure.trt$fperiod == 0]), sum(mean(seizure.trt$count[seizure.trt$fperiod == 1]), mean(seizure.trt$count[seizure.trt$fperiod == 2]),
       mean(seizure.trt$count[seizure.trt$fperiod == 3]), mean(seizure.trt$count[seizure.trt$fperiod == 4]))) ~ c(0:1),
      xaxt="n", ylab="Avg Number of Seizures", xlab="Time Period", type="l", 
      ylim=c(0, 50), main="Progabide")  
axis(side=1, at= 0:1,labels =c("pre*", "post*") )  

#########
# Only the difference between pre-treatment and afer week 4 of treatment
plot(	c(mean(seizure.plac$count[seizure.plac$fperiod == 0])/4, mean(seizure.plac$count[seizure.plac$fperiod == 4])) ~ c(0:1),
      xaxt="n", ylab="Avg Number of Seizures", xlab="Time Period", type="l", 
      ylim=c(0, 40), main="Placebo")  
axis(side=1, at= 0:1,labels =c("pre*", "week 4 of trt") )  

#then the progabide group 
plot(c(mean(seizure.trt$count[seizure.trt$fperiod == 0])/4, mean(seizure.trt$count[seizure.trt$fperiod == 4])) ~ c(0:1),
     xaxt="n", ylab="Avg Number of Seizures", xlab="Time Period", type="l", 
     ylim=c(0, 40), main="Progabide")  
axis(side=1, at= 0:1,labels =c("pre*", "week 4 of trt") )  




# Plot residuals
par(mfrow=c(1,1))
plot(residuals(AQuad20), xlab = "Observation", ylab = "", main = "Residuals from the Adaptive Quadrature with 20 pts")


# full model
Aquad20 <- glmer(count~new_x + (1|id), family="poisson", data= nlong, offset=(I(log(length))), nAGQ = 20)
pres.Aquad20=residuals(Aquad20,"pearson")

# null model (no random effects)
int_only<- glm(count ~ 1, family="poisson", data= nlong, offset=log(length))
pres=residuals(int_only,"pearson")

# null model with random effects
Aquad20_int_only=glmer(count ~ 1 + (1|id), family="poisson", data= nlong, offset=(I(log(length))), nAGQ = 20)
pres.Aquad20_int_only=residuals(Aquad20_int_only,"pearson")


# compare residuals for null and full models
plot(pres,ylab="Pearson residuals", main="Pearson residuals of Null and Full Model",col="blue",pch=4)
points(residuals(Aquad20,"pearson"),col="red",pch=1)
legend(200,25,c("Null Model", "Full Model"),col=c("blue","red"),pch=c(4,1))


# compare null mixed and full mixed models with global test
anova(Aquad20_int_only,Aquad20,"LRT")

# dispersion index for the three models
dispar=sum(pres**2)/(nrow(nlong)-10) # overdispersed
dispar.Aquad20_int_only=var(residuals(Aquad20_int_only,"pearson")) # still overdispersed but better
dispar.Aquad20=var(residuals(Aquad20,"pearson")) # still overdispersed but better yet



# Overdispersion (?)
var(residuals(AQuad20, "pearson"))

which(residuals(AQuad20)==max(residuals(AQuad20)))
long[124,]
