########################################################################### 
###################### PLOTTING THE TMLE RESULTS - VEG #################### 
########################################################################### 


#-------------------------------------------------------------------------------------------------------------------------------------
#PREPARATION
#-------------------------------------------------------------------------------------------------------------------------------------

# INSTALL AND LOAD PACKAGES
packages <- c("foreach","doParallel","boot","rmutil","mvtnorm","gam","sandwich","ggplot2", "SuperLearner","xgboost",
              "devtools","glmnet","tmle","data.table","rpart","ranger","nnet","arm","earth","e1071","tidyverse",
              "sl3", "tlverse", "tmle3", "beepr")
for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}
for (package in packages) {
  library(package, character.only=T)
}

# READ IN THE VEG TMLE OBJECT WE CREATED IN STEP #1
pree_veg <- readRDS(file="I:/nuMoM2B/Papers/Diet ML heterogeneity 2020/Data/tmle3_veg_pree_v4_cont_fixedcov_missind_0_new.rds")
dat <- read.csv(file="I:/nuMoM2B/Papers/Diet ML heterogeneity 2020/Data/numom_small_new_fixedcov_missind_2021.10.18.csv", header=TRUE, sep=",")
# Y = outcome (pree_acog), A = exposure (f_totdens80), W = covariates
#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------------------
# PREPARE FOR BOOTSTRAPPING
#-------------------------------------------------------------------------------------------------------------------------------------

# GET ITEMS OUT OF TMLE3 OBJECT:

# PROPENSITY SCORES
g.cf_task1 <- pree_veg$tmle_task$generate_counterfactual_task("cf1",data.frame(A=1))
pihat_0 <- pree_veg$likelihood$get_likelihood(g.cf_task1,"A")

# OUTCOME PREDICTIONS
Q.cf_task0 <- pree_veg$tmle_params[[1]]$cf_likelihood_control$cf_tasks[[1]]
Q.cf_task1 <- pree_veg$tmle_params[[1]]$cf_likelihood_treatment$cf_tasks[[1]]
mu0 <- as.numeric(pree_veg$likelihood$get_likelihood(Q.cf_task0,"Y", "validation"))
mu1 <- as.numeric(pree_veg$likelihood$get_likelihood(Q.cf_task1,"Y", "validation"))
muhat_0 <-  mu1*(dat$v_totdens80_0) + mu0*(1-dat$v_totdens80_0) 

# TAKE THE AIPW EFF 
X_ <- dat$v_totdens80_0
Y_ <- dat$pree_acog
aipw_EFF <- as.numeric((((2*X_-1)*(Y_ - muhat_0))/((2*X_-1)*pihat_0 + (1-X_)) + mu1 - mu0))

# ADD AIPW TO DATASET SO WE ARE RESAMPLING IT TOO
# I know I'm reading the data in a second time -- may change later
analysis <- read.csv(file="I:/nuMoM2B/Papers/Diet ML heterogeneity 2020/Data/numom_small_new_fixedcov_missind_2021.10.18.csv", 
                     sep=",", header=TRUE)
analysis <- cbind(analysis, aipw_EFF)

# We have to get rid of missing in BMI or we get an error. I am doing mean imputation:
#mean(analysis$bmi, na.rm=TRUE)
#analysis <- analysis %>% replace_na(list(bmi = 26.09507))


#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------------------
# BOOTSTRAPPING ------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------
bootNum = 50
res <- NULL
for(jj in 1:bootNum){
  set.seed(jj)
  boot_dat <- analysis %>% sample_n(.,nrow(analysis), replace=T)
  
  # HAVE TO ADD AS.NUMERIC TO GET RID OF ERROR
  M_ = as.numeric(boot_dat$bmi)
  y = boot_dat$aipw_EFF
  #M0_ = cbind(0, M_)
  x <- M_
  D <- data.frame(y,x)
  
  #set.seed(123)
  mm_numom <- seq(min(M_), max(M_), by = 1)
  print(jj)
  
  folds=10
  index<-split(1:nrow(D),1:folds)
  splt <- lapply(1:folds,function(ind) D[index[[ind]],])
  SL.library <- c("SL.glm","SL.step", "SL.earth", "SL.gam", "SL.mean", "SL.bayesglm")
  
  # HAND CODED ALGORITHMS
  m1 <- lapply(1:folds,function(ii) gam(y~s(x,5),family="gaussian",data=rbindlist(splt[-ii])))
  m2 <- lapply(1:folds,function(ii) gam(y~s(x,4),family="gaussian",data=rbindlist(splt[-ii])))
  m3 <- lapply(1:folds,function(ii) gam(y~s(x,3),family="gaussian",data=rbindlist(splt[-ii])))
  m4 <- lapply(1:folds,function(ii) mean(rbindlist(splt[-ii])$y))
  m5 <- lapply(1:folds,function(ii) bayesglm(y~x,data=rbindlist(splt[-ii]),family=gaussian))
  
  
  p1 <- lapply(1:folds,function(ii) predict(m1[[ii]],newdata=rbindlist(splt[ii]),type="response"))
  p2 <- lapply(1:folds,function(ii) predict(m2[[ii]],newdata=rbindlist(splt[ii]),type="response"))
  p3 <- lapply(1:folds,function(ii) predict(m3[[ii]],newdata=rbindlist(splt[ii]),type="response"))
  p4 <- lapply(1:folds,function(ii) m4[[ii]])
  p5 <- lapply(1:folds,function(ii) predict(m5[[ii]],newdata=rbindlist(splt[ii]),type="response"))
  
  for(i in 1:folds){
    splt[[i]] <- cbind(splt[[i]],p1[[i]],p2[[i]],p3[[i]],p4[[i]],p5[[i]])
  }
  
  X <- data.frame(do.call(rbind,splt))[,-2]
  names(X) <- c("y","gam1","gam2","gam3","mean","bayesglm")
  head(X)
  
  SL.r <- nnls(cbind(X[,2],X[,3],X[,4],X[,5],X[,6]),X[,1])$x
  alpha <- as.matrix(SL.r/sum(SL.r))
  round(alpha,3)
  
  m1_new <- gam(y~s(x,5),family="gaussian",data=D)
  m2_new <- gam(y~s(x,4),family="gaussian",data=D)
  m3_new <- gam(y~s(x,3),family="gaussian",data=D)
  m4_new <- mean(D$y)
  m5_new <- bayesglm(y~x,data=D,family=gaussian)
  
  p1_new <- predict(m1_new,newdata=data.frame(x=mm_numom),type="response")
  p2_new <- predict(m2_new,newdata=data.frame(x=mm_numom),type="response")
  p3_new <- predict(m3_new,newdata=data.frame(x=mm_numom),type="response")
  p4_new <- m4_new
  p5_new <- predict(m5_new,newdata=data.frame(x=mm_numom),type="response")
  
  SL.predict <- cbind(p1_new,p2_new,p3_new,p4_new,p5_new)%*%as.matrix(round(alpha,3))
  
  #identical(SL.predict,as.matrix(p1_new))
  
  res <- rbind(res,cbind(jj,SL.predict,mm_numom))
}

head(res)
tail(res)


d <- tibble(boot = res[,1],psi = res[,2], m = res[,3])
head(d)

d <- d %>% mutate(m_round = round(m, 0))
#d_test <- d %>% filter(boot==1 | boot==3)
d <- na.omit(d)
#psi is missing (NA) for some of the bootstrap resamples 
#is it okay to exclude them? 

d_mean <- aggregate(d$psi, list(d$m_round), mean)
colnames(d_mean) <- c("bmi", "rd_mean")


d_sd <- aggregate(d$psi, list(d$m_round), sd)
colnames(d_sd) <- c("bmi","sd_mean") 
d_sd <- d_sd %>% select("sd_mean")

plot_dat <- cbind(d_mean, d_sd)
plot_dat <- plot_dat %>% mutate(sd_plus = rd_mean + sd_mean,
                                sd_minus = rd_mean - sd_mean)

beep("mario")
#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------------------
# PLOTTING ------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------
ggplot(plot_dat) + 
  geom_smooth(aes(y = rd_mean, x = bmi),color="black", alpha = 1, size = 0.9, se=F) +
  geom_ribbon(aes(ymin = sd_minus, ymax =sd_plus, x = bmi), linetype = 2, size = 1.2, color="gray20", alpha = 0.1) +
  labs(x = "BMI (kg/m^2)",
       y = "Adjusted risk difference") +
  geom_hline(yintercept = 0,color="red",size=.9) +
  theme_bw()+
  scale_x_continuous(limits = c(15,45), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.15, 0.05), expand = c(0, 0)) +
  theme(plot.title=element_text(size=75),
        plot.subtitle=element_text(size=40),
        axis.title.x=element_text(size=60, vjust=-0.5),
        axis.title.y=element_text(size=60),
        axis.text.x=element_text(size=60, vjust=-0.5),
        axis.text.y=element_text(size=60),
        plot.margin = unit(c(1,1,1,1), "cm"))
#p + geom_ribbon(aes(ymin = plot_dat$sd_minus, ymax = plot_dat$sd_plus, x = plot_dat$heix_tot), colour="gray23", linetype = 2, alpha = 0.1)

ggsave("I:/nuMoM2B/Diet ML heterogeneity 2020/Plots/Plots/Sara Versions/VegPREE_BMI continuous_fixedcov_missind_0_new.png", width = 20, height = 20, dpi = 350)
#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------------------
# EDIT TO PLOT FOR AJE PUBLICATION ---------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------

## CODE TO TURN HYPHENS IN NEGATIVE NUMBERS INTO TRUE MINUS SIGNS
unicode_minus = function(x) sub('^-', '\U2212', format(x))

## INSTALL CURRENT VERSIONS OF CMAP AND EXPORT
devtools::install_github("CMAP-REPOS/cmapplot", build_vignettes=TRUE)
devtools::install_github("tomwenseleers/export", build_vignettes=TRUE)
install.packages("export")

## LOAD THE NEWLY INSTALLED PACKAGES
library(cmapplot)
library(export)

ggplot(plot_dat) + 
  geom_smooth(aes(y = rd_mean, x = bmi),color="black", alpha = 1, size=gg_lwd_convert(1), se=F) +
  geom_ribbon(aes(ymin = sd_minus, ymax =sd_plus, x = bmi), linetype =2, size=gg_lwd_convert(1), color="black", alpha = 0.1) +
  labs(x = "BMI",
       y = "Adjusted Risk Difference") +
  geom_hline(yintercept = 0,color="gray56", size=gg_lwd_convert(1)) +
  theme_bw()+
  scale_x_continuous(labels = unicode_minus, limits = c(15,45.1), expand = c(0, 0)) +
  scale_y_continuous(labels = unicode_minus, limits = c(-0.15, 0.05), expand = c(0, 0)) +
  theme(plot.title=element_text(size=14),
        plot.subtitle=element_text(size=12),
        # Hide panel borders and remove grid lines
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line color
        axis.line = element_line(color = "black", size=gg_lwd_convert(1)),
        # Change the appearance of the axes ticks and labels
        axis.title.x=element_text(size=14, margin = margin(t = 10), color = "black"),
        axis.title.y=element_text(size=14, margin = margin(r = 15), color = "black"),
        axis.text.x=element_text(size=14, color = "black"),
        axis.text.y=element_text(size=14, color = "black"),
        axis.ticks.length=unit(0.05,"inch"),
        axis.ticks=element_line(size=gg_lwd_convert(1), color="black"),
        plot.margin = unit(c(1,1,1,1), "cm"))
#p + geom_ribbon(aes(ymin = plot_dat$sd_minus, ymax = plot_dat$sd_plus, x = plot_dat$heix_tot), colour="gray23", linetype = 2, alpha = 0.1)

#ggsave("I:/nuMoM2B/Papers/Diet ML heterogeneity 2020/Plots/Plots/Sara Versions/VegPREE_BMI continuous_fixedcov_missind_0_new.png", width = 20, height = 20, dpi = 350)
ggsave("I:/nuMoM2B/Papers/Diet ML heterogeneity 2020/Plots/Plots/Sara Versions/VegPREE_BMI continuous_fixedcov_missind_0_new.png", width = 5, height = 5, dpi = 350)
ggsave("I:/nuMoM2B/Papers/Diet ML heterogeneity 2020/Plots/Plots/Sara Versions/VegPREE_BMI continuous_fixedcov_missind_0_new.pdf", width = 5, height = 5, dpi = 350)
graph2eps(file="I:/nuMoM2B/Papers/Diet ML heterogeneity 2020/Plots/Plots/Sara Versions/VegPREE_BMI continuous_fixedcov_missind_0_new.eps", width=5, height=5, fallback_resolution=600)

#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------









######### Read in the TMLE object we created
pree_veg <- readRDS(file="I:/nuMoM2B/Diet ML heterogeneity 2020/Data/tmle3_veg_pree_v4_cont_cov_missind_1.rds")
dat <- read.csv(file="I:/nuMoM2B/Diet ML heterogeneity 2020/Data/numom_small_new_fixedcov_missind_2021.09.01.csv", header=TRUE, sep=",")
# Y = outcome (pree_acog), A = exposure (f_totdens80), W = covariates
#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------------------
# PREPARE FOR BOOTSTRAPPING
#-------------------------------------------------------------------------------------------------------------------------------------
# Get stuff out of the tmle3 object:

# Propensity scores
g.cf_task1 <- pree_veg$tmle_task$generate_counterfactual_task("cf1",data.frame(A=1))
pihat_0 <- pree_veg$likelihood$get_likelihood(g.cf_task1,"A")

# Outcome predictions
Q.cf_task0 <- pree_veg$tmle_params[[1]]$cf_likelihood_control$cf_tasks[[1]]
Q.cf_task1 <- pree_veg$tmle_params[[1]]$cf_likelihood_treatment$cf_tasks[[1]]
mu0 <- as.numeric(pree_veg$likelihood$get_likelihood(Q.cf_task0,"Y", "validation"))
mu1 <- as.numeric(pree_veg$likelihood$get_likelihood(Q.cf_task1,"Y", "validation"))
muhat_0 <-  mu1*(dat$v_totdens80_1) + mu0*(1-dat$v_totdens80_1) 

# Take the AIPW EFF 
X_ <- dat$v_totdens80_1
Y_ <- dat$pree_acog
aipw_EFF <- as.numeric((((2*X_-1)*(Y_ - muhat_0))/((2*X_-1)*pihat_0 + (1-X_)) + mu1 - mu0))

# Add AIPW to data so we're resampling it too
# I know I'm reading the data in a second time -- may change later
analysis <- read.csv(file="I:/nuMoM2B/Diet ML heterogeneity 2020/Data/numom_small_new_fixedcov_missind_2021.09.01.csv", 
                     sep=",", header=TRUE)
analysis <- cbind(analysis, aipw_EFF)

# We have to get rid of missing in BMI or we get an error. I am doing mean imputation:
mean(analysis$bmi, na.rm=TRUE)
analysis <- analysis %>% replace_na(list(bmi = 26.09507))


#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------------------
# BOOTSTRAPPING ------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------
bootNum = 50
res <- NULL
for(jj in 1:bootNum){
  set.seed(jj)
  boot_dat <- analysis %>% sample_n(.,nrow(analysis), replace=T)
  
  # have to add as.numeric to get rid of error 
  M_ = as.numeric(boot_dat$bmi)
  y = boot_dat$aipw_EFF
  #M0_ = cbind(0, M_)
  x <- M_
  D <- data.frame(y,x)
  
  #set.seed(123)
  mm_numom <- seq(min(M_), max(M_), by = 1)
  print(jj)
  
  folds=10
  index<-split(1:nrow(D),1:folds)
  splt <- lapply(1:folds,function(ind) D[index[[ind]],])
  SL.library <- c("SL.glm","SL.step", "SL.earth", "SL.gam", "SL.mean", "SL.bayesglm")
  
  # hand coded algorithms
  m1 <- lapply(1:folds,function(ii) gam(y~s(x,5),family="gaussian",data=rbindlist(splt[-ii])))
  m2 <- lapply(1:folds,function(ii) gam(y~s(x,4),family="gaussian",data=rbindlist(splt[-ii])))
  m3 <- lapply(1:folds,function(ii) gam(y~s(x,3),family="gaussian",data=rbindlist(splt[-ii])))
  m4 <- lapply(1:folds,function(ii) mean(rbindlist(splt[-ii])$y))
  m5 <- lapply(1:folds,function(ii) bayesglm(y~x,data=rbindlist(splt[-ii]),family=gaussian))
  
  
  p1 <- lapply(1:folds,function(ii) predict(m1[[ii]],newdata=rbindlist(splt[ii]),type="response"))
  p2 <- lapply(1:folds,function(ii) predict(m2[[ii]],newdata=rbindlist(splt[ii]),type="response"))
  p3 <- lapply(1:folds,function(ii) predict(m3[[ii]],newdata=rbindlist(splt[ii]),type="response"))
  p4 <- lapply(1:folds,function(ii) m4[[ii]])
  p5 <- lapply(1:folds,function(ii) predict(m5[[ii]],newdata=rbindlist(splt[ii]),type="response"))
  
  for(i in 1:folds){
    splt[[i]] <- cbind(splt[[i]],p1[[i]],p2[[i]],p3[[i]],p4[[i]],p5[[i]])
  }
  
  X <- data.frame(do.call(rbind,splt))[,-2]
  names(X) <- c("y","gam1","gam2","gam3","mean","bayesglm")
  head(X)
  
  SL.r <- nnls(cbind(X[,2],X[,3],X[,4],X[,5],X[,6]),X[,1])$x
  alpha <- as.matrix(SL.r/sum(SL.r))
  round(alpha,3)
  
  m1_new <- gam(y~s(x,5),family="gaussian",data=D)
  m2_new <- gam(y~s(x,4),family="gaussian",data=D)
  m3_new <- gam(y~s(x,3),family="gaussian",data=D)
  m4_new <- mean(D$y)
  m5_new <- bayesglm(y~x,data=D,family=gaussian)
  
  p1_new <- predict(m1_new,newdata=data.frame(x=mm_numom),type="response")
  p2_new <- predict(m2_new,newdata=data.frame(x=mm_numom),type="response")
  p3_new <- predict(m3_new,newdata=data.frame(x=mm_numom),type="response")
  p4_new <- m4_new
  p5_new <- predict(m5_new,newdata=data.frame(x=mm_numom),type="response")
  
  SL.predict <- cbind(p1_new,p2_new,p3_new,p4_new,p5_new)%*%as.matrix(round(alpha,3))
  
  #identical(SL.predict,as.matrix(p1_new))
  
  res <- rbind(res,cbind(jj,SL.predict,mm_numom))
}

head(res)
tail(res)


d <- tibble(boot = res[,1],psi = res[,2], m = res[,3])
head(d)

d <- d %>% mutate(m_round = round(m, 0))
#d_test <- d %>% filter(boot==1 | boot==3)
d <- na.omit(d)
#psi is missing (NA) for some of the bootstrap resamples 
#is it okay to exclude them? 

d_mean <- aggregate(d$psi, list(d$m_round), mean)
colnames(d_mean) <- c("bmi", "rd_mean")


d_sd <- aggregate(d$psi, list(d$m_round), sd)
colnames(d_sd) <- c("bmi","sd_mean") 
d_sd <- d_sd %>% select("sd_mean")

plot_dat <- cbind(d_mean, d_sd)
plot_dat <- plot_dat %>% mutate(sd_plus = rd_mean + sd_mean,
                                sd_minus = rd_mean - sd_mean)

beep("mario")
#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------------------
# PLOTTING ------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------
ggplot(plot_dat) + 
  geom_smooth(aes(y = rd_mean, x = bmi),color="black", alpha = 1, size = 0.9, se=F) +
  geom_ribbon(aes(ymin = sd_minus, ymax =sd_plus, x = bmi), linetype = 2, size = 1.2, color="gray20", alpha = 0.1) +
  labs(x = "BMI (kg/m^2)",
       y = "Adjusted risk difference") +
  geom_hline(yintercept = 0,color="red",size=.9) +
  theme_bw()+
  scale_x_continuous(limits = c(15,45), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.15, 0.05), expand = c(0, 0)) +
  theme(plot.title=element_text(size=75),
        plot.subtitle=element_text(size=40),
        axis.title.x=element_text(size=60, vjust=-0.5),
        axis.title.y=element_text(size=60),
        axis.text.x=element_text(size=60, vjust=-0.5),
        axis.text.y=element_text(size=60),
        plot.margin = unit(c(1,1,1,1), "cm"))
#p + geom_ribbon(aes(ymin = plot_dat$sd_minus, ymax = plot_dat$sd_plus, x = plot_dat$heix_tot), colour="gray23", linetype = 2, alpha = 0.1)

ggsave("I:/nuMoM2B/Diet ML heterogeneity 2020/Plots/Plots/Sara Versions/VegPREE_BMI continuous_addcov_missind_1.png", width = 20, height = 20, dpi = 350)
#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------
