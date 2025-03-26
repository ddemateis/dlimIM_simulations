library(Matrix)

#3 modifiers
set.seed(011425)
n_mod <- 3
corr_mat_raw <- matrix(runif(n_mod^2,0.5,0.7), ncol=n_mod) #correlation matrix 
corr_mat_raw[!upper.tri(corr_mat_raw)] <- 1 #set every element that is not the upper triangular to 1
Sigma <- corr_mat_raw * t(corr_mat_raw) #perform multiplication to create symmetric matrix
Sigma <- nearPD(Sigma)$mat #make sure it is positive definite

save(Sigma, file = "modifier_Sigma_3.rda")

#10 modifiers
set.seed(011425)
n_mod <- 10
corr_mat_raw <- matrix(runif(n_mod^2,0.5,0.7), ncol=n_mod) #correlation matrix 
corr_mat_raw[!upper.tri(corr_mat_raw)] <- 1
Sigma <- corr_mat_raw * t(corr_mat_raw) 
Sigma <- nearPD(Sigma)$mat

save(Sigma, file = "modifier_Sigma_10.rda")

#50 modifiers
set.seed(011425)
n_mod <- 50
corr_mat_raw <- matrix(runif(n_mod^2,0.5,0.7), ncol=n_mod) #correlation matrix 
corr_mat_raw[!upper.tri(corr_mat_raw)] <- 1
Sigma <- corr_mat_raw * t(corr_mat_raw) 
Sigma <- nearPD(Sigma)$mat

save(Sigma, file = "modifier_Sigma_50.rda")
