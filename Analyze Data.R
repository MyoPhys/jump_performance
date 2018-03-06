## This code will analyze the scaling of jump velocity and power from an intraspecific size range of Cuban tree frogs
# The first required file is "Jump Performance Data v2.csv"
# The second required file is "Dissection Data.csv"

require(pspline)
require(caTools)
require(RColorBrewer)
require(scales)

# Load in the data file for jump performance
# Use file "Jump Performance Data v2.csv"
filename <- file.choose()
data <- (read.table(filename, sep = ",", skip=2))

# Define the variables
Name <- data$V1
BodyMass <- data$V2
SVL <- data$V3
Temp <- data$V4
Velocity <- data$V5
Power <- data$V6

# Find relationship between total hindlimb muscle mass and body length
# Use file "Dissection Data.csv"
filename <- file.choose()
dissect.data <- (read.table(filename, sep = ",", skip=2))
y <- dissect.data$V4 # Hind limb muscle mass [g]
x <- dissect.data$V3 # Snout-vent length [mm]
datax <- data.frame(x=x,y=y)
datax <- datax[order(datax$x),]

# Create a linear model of log-transformed data for seed values for non-linear regression
testmodel <- lm(log10(y)~log10(x), data=datax)
b_lr <- summary(testmodel)$coef[2,1]
a_lr <- summary(testmodel)$coef[1,1]

# Create a non-linear regression of muscle mass and body length using linear model for seed values
model_nlr = nls(y ~ a1 * x ^ a2, start = list(a1 = a_lr, a2 = b_lr),data=datax,
                control = nls.control(maxiter = 2000, warnOnly = TRUE))
a_nlr = coef(summary(model_nlr))[1, 1]
b_nlr = coef(summary(model_nlr))[2, 1]

#Apply this relationship to estimate total muscle mass of one hind limb for each individual
MMmodelInt <-  a_nlr
MMmodelSlope <- b_nlr
MuscleMass <- a_nlr*SVL^b_nlr
MuscleMass <- MuscleMass/1000 #in kg
#Assume 85% of mass is extensor (Peplowski and Marsh, 1997) and doulbe to include other limb
MuscleMass <- MuscleMass*0.85*2

# Convert jump power from body-mass-specific to muscle-mass-specific
MSpower <- Power*BodyMass/MuscleMass

# Divide up the variables into each temperature
mass30 <- BodyMass[which(Temp==30)]*1000
vel30 <- Velocity[which(Temp==30)]
power30 <- MSpower[which(Temp==30)]
mass20 <- BodyMass[which(Temp==20)]*1000
vel20 <- Velocity[which(Temp==20)]
power20 <- MSpower[which(Temp==20)]
mass10 <- BodyMass[which(Temp==10)]*1000
vel10 <- Velocity[which(Temp==10)]
power10 <- MSpower[which(Temp==10)]

## Start the plot Window
quartz()

## Plot the raw data with their linear regressions
#Plot the data points for 30C
plot(mass30,vel30, log ="xy", col ="firebrick", pch=16,ylim=c(1,4),xlim=c(1,35), ylab = "Takeoff Velocity (m/s)", xlab = "Body Mass (g)")

# Create the linear regression for jumps at 30C
model30 <- lm(log10(vel30)~log10(mass30))

#Plot the 95% confidence intervals for 30C
prd<-predict(model30,newdata=data.frame(x=log10(mass30)),interval = c("confidence"), 
             level = 0.95,type="response")
points(mass30,10^prd[,1],typ="l",col="firebrick")
colortray <- brewer.pal(4,"Set1")
testdata <- data.frame(mass30,lower = prd[,2], upper=prd[,3])
testdata <- testdata[order(mass30),]
polygon(c(testdata$mass30, rev(testdata$mass30)), c(10^testdata$lower, rev(10^testdata$upper)),
        col = alpha("firebrick",0.3), border = NA)

points(mass20,vel20,log ="xy", col = "black", pch=16)
# Create the linear regression for jumps at 20C
model20 <- lm(log10(vel20)~log10(mass20))

#Plot the 95% confidence intervals for 20C
prd<-predict(model20,newdata=data.frame(x=log10(mass20)),interval = c("confidence"), 
             level = 0.95,type="response")
points(mass20,10^prd[,1],typ="l",col="black")
colortray <- brewer.pal(4,"Set1")
testdata <- data.frame(mass20,lower = prd[,2], upper=prd[,3])
testdata <- testdata[order(mass20),]
polygon(c(testdata$mass20, rev(testdata$mass20)), c(10^testdata$lower, rev(10^testdata$upper)),
        col = alpha("black",0.3), border = NA)

points(mass10,vel10,log ="xy", col = "royalblue4", pch=16)
# Create the linear regression for jumps at 10C
model10 <- lm(log10(vel10)~log10(mass10))

#Plot the 95% confidence intervals for 10C
prd<-predict(model10,newdata=data.frame(x=log10(mass10)),interval = c("confidence"), 
             level = 0.95,type="response")
points(mass10,10^prd[,1],typ="l",col="royalblue4")
colortray <- brewer.pal(4,"Set1")
testdata <- data.frame(mass10,lower = prd[,2], upper=prd[,3])
testdata <- testdata[order(mass10),]
polygon(c(testdata$mass10, rev(testdata$mass10)), c(10^testdata$lower, rev(10^testdata$upper)),
        col = alpha("royalblue4",0.3), border = NA)

legend("bottomright", legend=c("10C", "20C", "30C"), col=c("royalblue4", "black", "firebrick"), pch=19, bty="n")

# Run the regressions to test effects of size and temperature
velmodel <- aov(log10(Velocity)~log10(BodyMass)*as.factor(Temp))

# create a dataframe containing only two temperature to get Q10 for each interval
compdata <- data.frame(Velocity,BodyMass,Temp)
temp1 <- 30
temp2 <- 20
compdata1 <- compdata[which(compdata$Temp==temp1|compdata$Temp==temp2),]
compmodel <- lm(log10(Velocity)~log10(BodyMass)+Temp,dat=compdata1)
VelHighQ10 <- 10^(10*compmodel$coefficients[3])

compdata <- data.frame(Velocity,BodyMass,Temp)
temp1 <- 20
temp2 <- 10
compdata1 <- compdata[which(compdata$Temp==temp1|compdata$Temp==temp2),]
compmodel <- lm(log10(Velocity)~log10(BodyMass)+Temp,dat=compdata1)
VelLowQ10 <- 10^(10*compmodel$coefficients[3])

## Plot for Power

## Start the plot Window
quartz()

## Plot the raw data with their linear regressions
#Plot the data points for 30C
plot(mass30,power30, log ="xy", col ="firebrick", pch=16,ylim=c(50,4000),xlim=c(1,35), ylab = "Peak Muscle Mass Specific Power (W/kg)", xlab = "Body Mass (g)")

# Create the linear regression for jumps at 30C
model30 <- lm(log10(power30)~log10(mass30))

#Plot the 95% confidence intervals for 30C
prd<-predict(model30,newdata=data.frame(x=log10(mass30)),interval = c("confidence"), 
             level = 0.95,type="response")
points(mass30,10^prd[,1],typ="l",col="firebrick")
colortray <- brewer.pal(4,"Set1")
testdata <- data.frame(mass30,lower = prd[,2], upper=prd[,3])
testdata <- testdata[order(mass30),]
polygon(c(testdata$mass30, rev(testdata$mass30)), c(10^testdata$lower, rev(10^testdata$upper)),
        col = alpha("firebrick",0.3), border = NA)

points(mass20,power20,log ="xy", col = "black", pch=16)
# Create the linear regression for jumps at 20C
model20 <- lm(log10(power20)~log10(mass20))

#Plot the 95% confidence intervals for 20C
prd<-predict(model20,newdata=data.frame(x=log10(mass20)),interval = c("confidence"), 
             level = 0.95,type="response")
points(mass20,10^prd[,1],typ="l",col="black")
colortray <- brewer.pal(4,"Set1")
testdata <- data.frame(mass20,lower = prd[,2], upper=prd[,3])
testdata <- testdata[order(mass20),]
polygon(c(testdata$mass20, rev(testdata$mass20)), c(10^testdata$lower, rev(10^testdata$upper)),
        col = alpha("black",0.3), border = NA)

points(mass10,power10,log ="xy", col = "royalblue4", pch=16)
# Create the linear regression for jumps at 10C
model10 <- lm(log10(power10)~log10(mass10))

#Plot the 95% confidence intervals for 10C
prd<-predict(model10,newdata=data.frame(x=log10(mass10)),interval = c("confidence"), 
             level = 0.95,type="response")
points(mass10,10^prd[,1],typ="l",col="royalblue4")
colortray <- brewer.pal(4,"Set1")
testdata <- data.frame(mass10,lower = prd[,2], upper=prd[,3])
testdata <- testdata[order(mass10),]
polygon(c(testdata$mass10, rev(testdata$mass10)), c(10^testdata$lower, rev(10^testdata$upper)),
        col = alpha("royalblue4",0.3), border = NA)

legend("bottomright", legend=c("10C", "20C", "30C"), col=c("royalblue4", "black", "firebrick"), pch=19, bty="n")


# Scaling Equation Coefficients for Peak MS Power
PPmodelInt <- log10(328)
PPmodelSlope <- 0    

# Define a range of body lengths
BodyLengths <- seq(from = 20, to=95, by=(95-20)/100)

# Predict PeakPower for range of Body Masses
PeakMSPowers <- 10^((log10(BodyLengths)*PPmodelSlope)+PPmodelInt)

# Predict Peak Powers for a range of BodyMasses with positive and negative 10C change
R1 <- PeakMSPowers
Temp1 <- 20
Temp2 <- 10
Q10 <- 3.53
R2 <- R1*exp(log(Q10)/(10/(Temp2-Temp1)))
PeakPowers10C <- R2

R1 <- PeakMSPowers
Temp1 <- 20
Temp2 <- 30
Q10 <- 1.51
R2 <- R1*exp(log(Q10)/(10/(Temp2-Temp1)))
PeakPowers30C <- R2

# Predict Body Mass based on Body Length
y <- BodyMass # Hind limb muscle mass [g]
x <- SVL # Snout-vent length [mm]
datax <- data.frame(x=x,y=y)
datax <- datax[order(datax$x),]

# Create a linear model of log-transformed data for seed values for non-linear regression
testmodel <- lm(log10(y)~log10(x), data=datax)
b_lr <- summary(testmodel)$coef[2,1]
a_lr <- summary(testmodel)$coef[1,1]

# Create a non-linear regression of muscle mass and body length using linear model for seed values
model_nlr = nls(y ~ a1 * x ^ a2, start = list(a1 = a_lr, a2 = b_lr),data=datax,
                control = nls.control(maxiter = 2000, warnOnly = TRUE))
a_nlr = coef(summary(model_nlr))[1, 1]
b_nlr = coef(summary(model_nlr))[2, 1]

# Scaling Equation Coefficients for Body Mass
BMmodelInt <- log10(a_nlr)
BMmodelSlope <- b_nlr

# Predict Body Masses
BodyMasses <- 10^((log10(BodyLengths)*BMmodelSlope)+BMmodelInt)
BodyMasses <- BodyMasses*1000

points(BodyMasses,PeakPowers30C,col="firebrick",type="l")
points(BodyMasses,PeakPowers10C,col="royalblue4",type="l")
points(BodyMasses,PeakMSPowers,col="black",type="l")


# Run the regressions to get effects of temperature and size
powermodel <- aov(log10(MSpower)~log10(BodyMass)+as.factor(Temp))

# create a dataframe containing only two temperature to get Q10 for each interval
compdata <- data.frame(Power,BodyMass,Temp)
temp1 <- 30
temp2 <- 20
compdata1 <- compdata[which(compdata$Temp==temp1|compdata$Temp==temp2),]
compmodel <- lm(log10(Power)~log10(BodyMass)+Temp,dat=compdata1)
PowerHighQ10 <- 10^(10*compmodel$coefficients[3])

compdata <- data.frame(Power,BodyMass,Temp)
temp1 <- 20
temp2 <- 10
compdata1 <- compdata[which(compdata$Temp==temp1|compdata$Temp==temp2),]
compmodel <- lm(log10(Power)~log10(BodyMass)+Temp,dat=compdata1)
PowerLowQ10 <- 10^(10*compmodel$coefficients[3])

VelLowQ10
VelHighQ10
PowerLowQ10
PowerHighQ10
summary(velmodel)
summary(powermodel)
velmodel$coefficients[2]
powermodel$coefficients[2]


