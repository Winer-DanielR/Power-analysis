#Tribulus analysis#
#Daniel Reyes Corral#


#The question here is what is the minimum power necesary to view significant differences on morphology betweeen three 
#groups of mericarps (small, large and two-spined). Morphology means a series of linear measurements on length, width,
#depth of mericarp, spine length, spine distance.
dataset <- read_csv("Data/Eco evo experiment.csv")
View(dataset)
head(dataset)

##############PCA##############
#Multidimensional variables on morphology I could use a PCA to see the patterns in the data to help explain the variance
myPr <- prcomp(dataset[,4:8], scale = TRUE)
myPr
summary(myPr)
plot(myPr, type = "l")
biplot(myPr, scale = 0)
#PCA on hardness measurements
hardPr <- prcomp(dataset[,12:15], scale = TRUE)
hardPr
summary(hardPr)
plot(hardPr, type = "l")
biplot(hardPr, scale = 0)

#Extracting PC scores
str(myPr)
myPr$x

str(hardPr)
hardPr$x
#Adding PC scores to dataset. I created two new sets to avoid confusion between PC scores from morphology and hardness
dataset2 <- cbind(dataset, myPr$x[,1:2]) #Morphology PCs
dataset3 <- cbind(dataset, hardPr$x[,1:2]) #Hardness PCs
head(dataset2)
head(dataset3)
#plot with ggplot

ggplot(dataset2, aes(PC1, PC2, col = treatment, fill = treatment)) + stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) + geom_point(shape = 21, col = "black")

ggplot(dataset3, aes(PC1, PC2, col = treatment, fill = treatment)) + stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) + geom_point(shape = 21, col = "black")

#Correlations between variables and PCs
cor(dataset[,4:8], dataset2[,16:17])
cor(dataset[,4:8], dataset3[,16:17])

###########DEFINING THE MODEL#############

#model 1 is using PC scores for morphology
model1 <- lmer(PC1 ~ treatment + lower_spine + (1|parcel) + (1|mericarp), data = dataset2) #Morphology PCA
summary(model1)
drop1(model1)

model2 <- lmer(PC1 ~ treatment + lower_spine + (1|parcel) + (1|mericarp), data=dataset3) #Hardness PCA
summary(model2)
drop1(model2)

## Checking assumptions model 1
plot(model1)
opar <- par(mfrow=c(2,2))
hist(dataset2$PC1)
plot(fitted(model1),resid(model1)) 
abline(h=0,lty=2,col="red")
qqnorm(resid(model1))
qqline(resid(model1), col="red")
hist(resid(model1)) 
par(opar)

#Checking assumptions model 2
plot(model2)
opar <- par(mfrow=c(2,2))
hist(dataset2$PC1)
plot(fitted(model2),resid(model2)) 
abline(h=0,lty=2,col="red")
qqnorm(resid(model2))
qqline(resid(model2), col="red")
hist(resid(model)) 
par(opar)

###########POWER ANALYSIS############

#Check fixed effect intercepts. I am interested in small mericarps
fixef(model1)["treatmenttwo spines"] <- 0.05
fixef(model2)["treatmenttwo spines"] <- 0.05

#Power analysis for both models
sim_morph <- powerSim(model1, nsim = 20)
sim_morph

sim_hard <- powerSim(model2, nsim = 20) 
sim_hard

#Power increases with effect size.
curve_morph <- powerCurve(model1, along = "parcel")
plot(curve_morph)

curve_hard <- powerCurve(model2, within = "parcel")
plot(curve_hard)

