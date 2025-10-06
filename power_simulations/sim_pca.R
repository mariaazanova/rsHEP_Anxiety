library(MASS) # simulate data
library(fGarch) # skewed normal dist
library(dplyr) # for summary
library(ggplot2) # visualise 
library(factoextra) # percent explained for pca
set.seed(423) # fix rng


# input for simulations ####
num_sim <- 1000 
n_vars <-11
var_names <- c("HEP","anxiety","HR","HRV","bloodPressure","BMI","depressionScore","alphaPower","sex","age","time")
mu <- c(0.4, 5, 75, 90, 130, 26, 15, 0, 0, 60, 50) # means of dist for each var, can be rep(0,n_vars)
sds <- c(0.1, 4, 8, 10, 8, 1, 6, 1, 1, 10, 17) # sd of each var, can be rep(1,n_vars)
vars <- sds^2;
cor_HEP_anx <- 0.3
cor_data <-     c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, # HEP
                  cor_HEP_anx, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, # anxiety
                  0.05, 0.1, 1, 0, 0, 0, 0, 0, 0, 0, 0, # HR
                  -0.1, -0.3, -0.6, 1, 0, 0, 0, 0, 0, 0, 0, # HRV
                  0.05, 0.1, 0.3, -0.4, 1, 0, 0, 0, 0, 0, 0, # bloodPressure
                  0.05, 0.2, 0.3, -0.3, 0.5, 1, 0, 0, 0, 0, 0, # BMI
                  0.01, 0.5, 0.2, -0.2, 0.2, 0.2, 1, 0, 0, 0, 0, # depressionScore
                  0.1, 0.1, 0.01, -0.01, 0.01, 0.01, 0.05, 1, 0, 0, 0, # alphaPower
                  0.1, 0.3, 0.05, -0.01, 0.05, 0.1, 0.1, 0.05, 1, 0, 0, # sex, female is 1
                  0.1, 0.2, 0.2, -0.1, 0.2, 0.1, 0.1, -0.1, 0.1, 1, 0, # age
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1); # time

# if we want to try different strength of correlation among vars:
cor_data_stronger <- abs(cor_data)^(2/3)
cor_data_stronger[cor_data<0] <- cor_data_stronger[cor_data<0]*(-1)
cor_data_weaker <- cor_data^2
cor_data_weaker[cor_data<0] <- cor_data_weaker[cor_data<0]*(-1)
cor_data_weakest <- cor_data^4
cor_data_weakest[cor_data<0] <- cor_data_weakest[cor_data<0]*(-1)

make_skewed <- c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE) # which vars to skew, can do an integer array
make_skewed <- c(1:n_vars)[make_skewed]
make_skewed_xi <- c(0, 3, 1, 1, 0.5, 0.5, 1.4, 0, 0, 0.6, 1.4) # xi for qsnorm
# vary: n subs and/or HEP-anxiety association (through correlation or regression coefficient)
sim_n_subs <- c(40, 100, 300, 500, 1000, 1500, 3000)
sim_cor <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4) # vary correlation strength
sim_coef <- c(0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035) # vary regression coefficient size, if needed

# functions ####

find_cov_mat <- function(n_vars,vars,cor_data,lower_tri=TRUE){
  
  cor_matrix <- matrix(data = cor_data, nrow = n_vars, ncol = n_vars, byrow = TRUE);
  
  if (lower_tri){
    cor_matrix[upper.tri(cor_matrix)]<-t(cor_matrix)[upper.tri(cor_matrix)]
  } else {
    cor_matrix[lower.tri(cor_matrix)]<-t(cor_matrix)[lower.tri(cor_matrix)]
  }
  
  cov_matrix <- matrix(0, nrow=n_vars, ncol=n_vars)
  for (i in 1:n_vars){
    for (j in 1:n_vars){
      cov_matrix[i,j] <- cor_matrix[i,j]*sqrt(vars[i]*vars[j])
    }
  }
  return(cov_matrix)
}

simulate_data <- function(n_vars,n_subs,var_names,mu,cov_matrix,sim_reg=TRUE,skew=TRUE,make_skewed=c(),make_skewed_xi=c()){
  
  if (sim_reg){
    # if sim_reg=TRUE (default), simulate only IVs with mvrnorm, 
    # and generate HEP from them using reg coef later
    # add HEP column as a placeholder
    sim_data <- mvrnorm(n = n_subs, mu = mu[-1], Sigma = cov_matrix[-1,-1]);
    sim_data <- cbind(HEP = 0, sim_data)
    sim_data <- setNames(data.frame(sim_data),var_names);
  } else {
    # if sim_reg=FALSE, simulate all vars (including HEP) with mvrnorm through cor matrix
    sim_data <- mvrnorm(n = n_subs, mu = mu, Sigma = cov_matrix);
    sim_data <- setNames(data.frame(sim_data),var_names);
  }
  
  # make sex binary
  p_sex <- pnorm(sim_data$sex, mean = 0, sd = 1)
  sim_data$sex <- qbinom(p=p_sex, size=1, prob=0.5)
  
  if (skew){
    for (i in make_skewed){
      pvar<-pnorm(sim_data[,i],mean=mu[i], sd=sds[i])
      sknormvar <- qsnorm(pvar, mean=mu[i], sd=sds[i], xi=make_skewed_xi[i])
      sim_data[,i]<-sknormvar
    }
  }
  return(sim_data)
}

# to set up regression results, just supply sim_coef to the function instead of sim_cor
# overall there are 3 options for regression result: no controls, all controls, PCA
set_up_sim_res <- function(num_sim,sim_n_subs,sim_cor){
  
  nrows <- num_sim*3*length(sim_n_subs)*length(sim_cor) 
  cols <- c("sim_num","controls","n_sub","sim_value","result_coef", "pval", "sig","num_comp")
  sim_res <- setNames(data.frame(matrix(ncol = length(cols), nrow = nrows)),cols)
  
  sim_res["sim_num"] <- rep(1:num_sim,3*length(sim_n_subs)*length(sim_cor))
  sim_res["controls"] <- rep(rep(c('no','all','pca'),each=num_sim),length(sim_n_subs)*length(sim_cor))
  sim_res["n_sub"] <- rep(sim_n_subs, each=3*num_sim*length(sim_cor))
  sim_res["sim_value"] <- rep(rep(sim_cor, each=3*num_sim),length(sim_cor))
  
  return(sim_res)
}

# power simulations ####

correlation_data<-cor_data_stronger
correlation_data<-cor_data_weaker
correlation_data<-cor_data_weakest
correlation_data<-cor_data

sim_reg=FALSE # generating HEP through correlation matrix with other variables
sim_values<-sim_cor


sim_res<-set_up_sim_res(num_sim,sim_n_subs,sim_values)
for (sim_value in sim_values){
  
  if(!sim_reg){
    correlation_data[12]<-sim_value  
    cov_matrix<-find_cov_mat(n_vars,vars,correlation_data,lower_tri=TRUE)
  }
  
  for (n_subs in sim_n_subs){
    for (s in 1:num_sim){
      
      find_simulation <- (sim_res$sim_num==s)&(sim_res$n_sub==n_subs)&(sim_res$sim_value==sim_value)
      find_no_controls <- (sim_res$controls=='no')&find_simulation
      find_controls <- (sim_res$controls=='all')&find_simulation
      find_pca <- (sim_res$controls=='pca')&find_simulation
      
      sim_data <- simulate_data(n_vars,n_subs,var_names,mu,cov_matrix,sim_reg,skew=TRUE,make_skewed,make_skewed_xi)
      
      # simulate HEP here if doing a regression (doesnt include time-anxiety interaction, although we could add something)
      if(sim_reg)
      {
        HEP <- rnorm(n_subs, mean=mu[1]-0.2,sd=sds[1]) 
        # subtracted 0.2 from the mean used to generate HEP with correlation -- to make them comparable in size to correlation approach
        sim_data$HEP <- HEP + sim_value*sim_data$anxiety + rowSums(mapply(`*`, sim_data[,3:11], betas))
      }
      
      # add no_controls result
      lmm<-lm(HEP~anxiety, data=sim_data)
      sim_res[find_no_controls,"result_coef"] <- summary(lmm)$coefficients['anxiety',1]
      sim_res[find_no_controls,"pval"] <- summary(lmm)$coefficients['anxiety',4]
      # add all_controls 
      lmm<-lm(HEP~anxiety+HR+HRV+bloodPressure+BMI+depressionScore+alphaPower+sex+age+time*anxiety, data=sim_data)
      sim_res[find_controls,"result_coef"] <- summary(lmm)$coefficients['anxiety',1]
      sim_res[find_controls,"pval"] <- summary(lmm)$coefficients['anxiety',4]  
      
      # add PCA result 
      # pca doesnt take into account time-anxiety interaction
      res.pca <-prcomp(sim_data[3:11])
      num_comp <- sum(get_eigenvalue(res.pca)$variance.percent>1) # save this value? can also make it a var
      res_pca<-data.frame(res.pca$x[,1:num_comp]) 
      res_pca[c('HEP','anxiety')] <- sim_data[c('HEP','anxiety')]
      lmm<-lm(HEP~., data=res_pca)
      sim_res[find_pca,"result_coef"] <- summary(lmm)$coefficients['anxiety',1]
      sim_res[find_pca,"pval"] <- summary(lmm)$coefficients['anxiety',4]  
      sim_res[find_pca,"num_comp"] <- num_comp
    }
  }
}
sim_res["sig"] <- as.integer(sim_res["pval"]<0.05)

# save the result of simulations
save_path = "/data/p_02865/output/simulations/pca/"
write.csv(sim_res, "sim_res_reg.csv", row.names=FALSE)

# visualise ####

mean_data <- group_by(sim_res, controls, n_sub, sim_value) %>%
  summarise(sig = mean(sig), coef = mean(result_coef))

p1<-ggplot(mean_data[mean_data$controls=='no',], aes(x=sim_value, y=sig, group=factor(n_sub)))+
  geom_line(aes(color=factor(n_sub))) +
  geom_hline(yintercept=0.95, linetype="dashed", color = "red")+
  labs(x = "Simulated association value",y = "Power", color = "Subjects") +
  ggtitle('Power estimation for models with no controls')
p1

p2<-ggplot(mean_data[mean_data$controls=='all',], aes(x=sim_value, y=sig, group=factor(n_sub)))+
  geom_line(aes(color=factor(n_sub))) +
  geom_hline(yintercept=0.95, linetype="dashed", color = "red")+
  labs(x = "Simulated association value",y = "Power", color = "Subjects") +
  ggtitle('Power estimation for models with all controls')
p2

p3<-ggplot(mean_data[mean_data$controls=='pca',], aes(x=sim_value, y=sig, group=factor(n_sub)))+
  geom_line(aes(color=factor(n_sub))) +
  geom_hline(yintercept=0.95, linetype="dashed", color = "red")+
  labs(x = "Simulated association value",y = "Power", color = "Subjects") +
  ggtitle('Power estimation for models with PCA')
p3