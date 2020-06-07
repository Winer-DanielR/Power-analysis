##############################################################
################# Comp 598 Final Project #####################
###################    Power Analysis  #######################
##############################################################
#https://cran.r-project.org/web/packages/simr/vignettes/fromscratch.html
#Covariates and parameters

x <- 1:10
g <- letters[1:3]

X <- expand.grid(x=x, g=g)
show(X)

#Specify some fixed and random paramenters

b <- c(2, -0.1) #fixed intercept and slope
V1 <- 0.5 #random intercept variance
V2 <- matrix(c(0.5,0.05,0.05,0.1), 2) #random intercept and slope variance-covariance matrix
s <- 1 #residual standard deviation

#Build a model object
model1 <- makeLmer(y ~ x + (1|g), fixef = b, VarCorr = V1, sigma = s, data = X)
print(model1)
  
model2 <- makeGlmer(z ~ x + (x|g), family = "poisson", fixef = b, VarCorr = V2, data = X)
print(model2)

#Start the power analysis
#Now we have "pilot" models, which can be used with simr.
### I need to define the models and also define the pilot data ###

powerSim(model1, nsim=20)

powerSim(model2, nsim=20)
