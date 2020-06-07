#Power analysis with simR. Introduction
#Determine power required/sample size for linear mixed effect models. The package lme4 is used for modelling.
#Simulating data. In this example is two groups one experiment the other control. Measurements were taking inmediate after intervention and then after
#later. Summary 5 classes and 10 participants per class.

subj <- factor(1:10)
class_id <- letters[1:5]
time <- 0:2
group <- c("control", "intervention")

subj_full <- rep(subj, 15)
class_full <- rep(rep(class_id, each=10),3)
time_full <- rep(time, each=50)
group_full <- rep(rep(group, each = 5), 15)

covars <- data.frame(id = subj_full, class = class_full, treat = group_full, time=factor(time_full))

covars

#Define the model. In this case y ~ treatment + time + treatment x time + (1|class/id) + e
#In this example parameter values were chosen arbitrarily but should be based on literature or experience as much as possible.

#Intercep and slopes for intervention, time1, time2, intervention:time1, intervention:time2
fixed <- c(5, 0, 0.1, 0.2, 1, 0.9)

#Random intercepts for participants clustered by class
rand <- list(0.5, 0.1)
rand

#Residual variance
res <- 2

#Create the model. The makeLmer function from the simr package allows us to combine all this information to create a fitted lmer model from scratch
model <- makeLmer(y ~ treat*time + (1|class/id), fixef = fixed, VarCorr = rand, sigma = res, data = covars)
model

#Power analysis. Once you fitted the lmer model, whether it was fitted to real data or created from scratch, you can use that to simulate new data
#and assess the required sample size.
#powerSim function allows us to estimate the power to detect a specific effect in the model. Here we are interested in the effect of the intervention
#fcompare function compares the model to another one that does not include the treatment variable. This allows us to specify the fixed effects in the
#alternative model. All random effects will be assumed to be the same as in the orginal model.
sim_treat <- powerSim(model, nsim = 100, test = fcompare(y~time))
sim_treat

#We can test for the effect of time in the same way
sim_time <- powerSim(model, nsim = 100, test = fcompare(y~treat))
sim_time

#Changing effect size. To compute for an effect size that is smaller than the one observed in a pilot study, or to study power for a range of effect
#sizes
model_large <- model
fixef(model_large)["treatintervention:time1"] <- 2

sim_treat_large <- powerSim(model_large, nsim = 100, test = fcompare(y~time))
sim_treat_large

#Changing the number of classes used in study. To study the effect an increase in sample size has on our ability to detect the effec of interest
#We can increase the number of levels for one of the factors in our model. This can be achieved with the extend function. For example increase
#the number of classes from 5 to 20.
model_ext_class <- extend(model, along = "class", n = 20)
model_ext_class

#We can visualize the effect that varying the number of classes has on the power to detect the intervention effect by asking simr to plot
#a power curve
sim_treat_class <- powerSim(model_ext_class, nsim = 100, test = fcompare(y~time))
sim_treat_class

p_curve_test <- powerCurve(model_ext_class, test = fcompare(y~time), along = "class")
plot(p_curve_test)

#Changing the number of participants per class
#Instead of increasing the number of classes we could increase the number of participants per class. This can be achieved with the within argument
#to extend. Here we are extending the number of sutdent for each treatment at each time potin within each class to 20.

model_ext_subj <- extend(model, within="class+treat+time", n=20)
model_ext_subj

#Once again we can obtain a power estimate at the increased sample size as well as a power curve taht combine estimates at different sample sizes
sim_treat_subj <- powerSim(model_ext_subj, nsim = 100, test=fcompare(y~time))
sim_treat_subj

p_curve_treat <- powerCurve(model_ext_subj, test = fcompare(y~time), within = "class+treat+time", breaks = c(5, 10, 15, 20))
plot(p_curve_treat)

#Changing numbers of classes and participants per class
#We may want to study models with changes to multiple sample size components. The relevant changes can be applied to the model sucessively
model_final <- extend(model, along = "class", n=8)
model_final <- extend(model_final, within = "class+treat+time", n=10)

sim_final <- powerSim(model_final, nsim=100, test = fcompare(y~time))
sim_final
