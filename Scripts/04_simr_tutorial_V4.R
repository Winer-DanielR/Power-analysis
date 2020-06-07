#Tutorial from paper#
#Environmental data with a response variable z (bird abundance) measured at 10 levels of the continous fized effect variable x (study year)
#on three groups (study site), there is also a continous response variable y, which is not used in this tutorial.

#Start by fitting a model in lme4 to the dataset. In this case we have a random intercept model, where each group (g) has its own intercept
#but the groups share a common trend.
simdata
model <- glmer(z ~ x + (1|g), family = "poisson", data = simdata)
summary(model)

#In this case the effect size for x is -0.11, which is significant at the 0.01 level. A proper model would have a larger number of groups, and
#and would consider problems such as overdispersion. 

#If the effect size is real, would we have enough power to expect a positive result?

#Power increases with effect size.
#If we want to detect a slope of -0.05 you can change the fixed effect size from -0.11 to -0.05 as follows:
fixef(model)["x"]

fixef(model)["x"] <- -0.05

#Running the power analysis
powerSim(model)

#The power to reject the null hypothesis of zero trend in x is about 33%. Usually 80% power is considered adequate.

#Increase sample size
#A small pilot study often will not have enough power to detect a small effect, but a larger study might. In SIMR
#the extend function can be used to add rows to a dataframe. In this study the observations were 10 of x representing
#for example study years 1 - 10. What whould happen if we increse it to 20 years
model2 <- extend(model1, along = "x", n = 20)
powerSim(model2)

#Power analysis at a range of sample sizes
#When data collection is costly, the user might want to collect only as much data as needed. The Power Curve function
#can be used to explore trade-offs betweeen size and power.
#Minimum sample size required
#Can we reduce the number of samples while keeping the power above the usual 80%
pc2 <- powerCurve(model2)
print(pc2)
plot(pc2)

#If x were study year, we might be unwilling to wait for longer for our results. In this case increasing the number of 
#sites may be the best option. These two analyses start back with our original model1 which had 10 years.
#Adding more groups
model3 <- extend(model1, along = "g", n=15)
pc3 <- powerCurve(model3, along = "g")
plot(pc3)

#To reach 80% power, we need at least 11 sites.

#Increasing the size within groups
#We can replace the along argument to extend and powerCruve with the within argumen to increase the sample size
#within groups. Each group has only one obervation at each level of x and g. We can extend this to five observations
#per site per year.
model4 <- extend(model1, within = "x + g", n=5)
pc4 <- powerCurve(model4, within = "x + g", breaks = 1:5)
print(pc4)

