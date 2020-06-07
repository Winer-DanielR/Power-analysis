#Test examples from the simr tutorial package#
#This provides examples of some of the hypothesis tests that cna be specified in simr. The function dotest can be used to apply a test
#to an input model, which lets you check that the test works before running a power simulation

#Binomial GLMM with a categorical predictor
#glamer dataframe cbpp. An observation variable is added to allow for overdispersion. Note that the respose is specified using cbind-simr expects a binomial model

cbpp$obs <- 1:nrow(cbpp)
View(cbpp)
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1|herd) + (1|obs), data = cbpp, family = binomial)
summary(gm1)

#Note that period is a categorical variable with four levels, which enters the model as 3 dummy variables. To test all 3 dummy variables simultanously,
#you can use a likelihood ratio test
doTest(gm1, fixed("period", "lr"))

#If you were interested in the significance for the dummy variable period2 you could use a Z test. 
doTest(gm1, fixed("period2", "z"))

#Suppose your model also has a continous predictor. You can use fized to choose which fixed effect to apply tests to.
gm2 <- glmer(cbind(incidence, size - incidence) ~ period + size + (1|herd), data = cbpp, family = binomial)
doTest(gm2, fixed("size", "z"))

#Once you have chosen your tests,  you can run a power analysis by replacing dotest with powerSim. 
#Dont forget to specify an appropiate effect size.
fixef(gm2)["size"] <- 0.05
powerSim(gm2, fixed("size", "z"), nsim=50)
