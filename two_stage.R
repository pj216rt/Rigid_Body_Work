#Two-stage fit of the rigid body model.
#
#Stage 1: for every frame, find the jaw's rotation R and translation T 
# from the 3 markers (Kabsch), and turn R into a rotation vector w.
#Stage 2: regress each motion curve on time, task, and the covariate with
# penalized splines (mgcv), plus random curves for subjects.

library(mgcv)
source("rotations.R")

marker_names <- c("incisor", "condyle_L", "condyle_R")
motion_names <- c("w1", "w2", "w3", "T1", "T2", "T3")

#load data
sim <- readRDS("sim_data.rds")
markers    <- sim$markers
covariates <- sim$covariates

# Stage 0: the jaw's own frame for each subject, from the rest frame