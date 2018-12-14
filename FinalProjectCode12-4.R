##### STT 420 Final Project #####
#####   Logistic Incedence  #####
##### Authors: Thomas Billman, Danielle Chaung, James Royal #####

##### Population Simulation #####
library(tidyverse)
#Population probabilities
gender  <- c(.52,.48) #f, m
age <- c(.19, .13, .26, .25, .16) #0-17, 18-24, 25-44, 45-64, 65+
race <- c(.817, .14, .016, .006) #W, B, Asian, Indian
mental <- c(.8,.2) #n, y
substance <- c(.906, .094) #n,y
perscription <- c(.8, .2) #n, y
hosp <- c(.927, .073) #n, y
rural <- c(.978, .022) #n, y

#Simulation
set.seed(420)
n <- 216430
gend <- rbernoulli(n, gender[2])
ages <- sample(1:5, n, replace = T, prob = age)
races <- sample(1:4, n, replace = T, prob = race)
mentals <- rbernoulli(n, mental[2])
substances <- rbernoulli(n, substance[2])
perscriptions <- rbernoulli(n, perscription[2])
hosps <- rbernoulli(n, hosp[2])
rurals <- rbernoulli(n,rural[2])

#Combination
dat <- data.frame(gend,
                  ages,
                  races,
                  mentals,
                  substances,
                  perscriptions,
                  hosps,
                  rurals)

##### Emergency Room Visits Simulation #####

#There were a total of 535 visits to New Hanover Medical Center in 2016
#so average p should be 535/n = .002471931
c <- 535

#gender parameter estimation
p <- .7 #proportion of cases that are male
pc <- round(p*c) #positive cases are the rounded product of these
nc <- c - pc #negative cases are the remainders
or <- pc*(n-sum(gend))/(sum(gend)*nc)
bgen <- log(or) # Beta for male gender 


#Age parameter estimation
  pops <- table(ages)/n
  risks <- c(0, .06, .55, .38, .01)
  nrisks <- round(risks * c)
  nrisks[5] <- nrisks[5] + 1 #Due to odd rounding
  
  agemat <- matrix(c(nrisks, table(ages)-nrisks), nrow = 2, byrow = T)
  agemat
  #Baseline defined as second group due to the first being degenerate
  
  bage1 <- -Inf #This approaches negative inifinity so that the p for this group goes to 0
  bage2 <- 0
  bage3 <- log(agemat[1,3]*agemat[2,2]/agemat[1,2]/agemat[2,3]) #Remaining logs odds ratios used to
  bage4 <- log(agemat[1,4]*agemat[2,2]/agemat[1,2]/agemat[2,4]) #To esimated beta factors
  bage5 <- log(agemat[1,5]*agemat[2,2]/agemat[1,2]/agemat[2,5])
  

#Race parameter estimation
  pops <- table(races)/n
  risks <- c(.85, .08, 0, .07)
  nrisks <- round(risks * c)
  
  racemat <- matrix(c(nrisks, table(races)-nrisks), nrow = 2, byrow = T)
  racemat
  #Baseline defined as second group due to the first being degenerate
  
  brace1 <- 0 #reference group
  brace2 <- log(racemat[1,2]*racemat[2,1]/racemat[1,1]/racemat[2,2]) #Remaining logs odds ratios used to
  brace3 <- log(racemat[1,3]*racemat[2,1]/racemat[1,1]/racemat[2,3])
  brace4 <- log(racemat[1,4]*racemat[2,1]/racemat[1,1]/racemat[2,4])


#Other parameter estimation
  #We already had odds ratios for these from previous studies
  bment <- log(3.94)
  bsubs <- log(5.24)
  bpers <- log(4.1)
  bhosp <- log(2.9)
  brurl <- log(.93)

#Concatenating into a vector
  betas <- c(bgen,
             bage1, bage2, bage3, bage4, bage5,
             brace1, brace2, brace3, brace4,
             bment,
             bsubs,
             bpers,
             bhosp,
             brurl)
#Testing with an intercept of 0 to see what the population expectation looks like
#This takes a few minutes
totbeta <- rep(0,n)

for(i in 1:n){
              print(i)
              totbeta[i] <- sum(gend[i]*betas[1],
                 betas[ages[i] + 1],
                 betas[races[i] + 6],
                 betas[11:15] %*% t(dat[i,4:8])
                 )
}  

probs <- 1/(1+exp(-totbeta))
ev <- sum(probs) # 129944.5
b0 <- log(535/ev)

tests <- NULL
for(b0 in -50:-2){
  ntotbeta <- totbeta + b0
  nprobs <- 1/(1+exp(-ntotbeta))
  nev <- sum(nprobs) # 129944.5
  tests <- c(tests,nev)
}
tests
#so b0 should be between -9 and -8?

tests <- NULL
vct <- -9 + .01*1:100
for(b0 in vct){
  ntotbeta <- totbeta + b0
  nprobs <- 1/(1+exp(-ntotbeta))
  nev <- sum(nprobs) # 129944.5
  tests <- c(tests,nev)
}
tests

# b0 of -8.69 gives an Expected value of 536.2 which is very close to 535

ntotbeta <- totbeta -8.69
nprobs <- 1/(1+exp(-ntotbeta))
dat <- cbind(dat,ntotbeta, nprobs)

set.seed(12042018)
resp <- sapply(1:dim(dat)[1],function(i){
  base::sample(c(1,0),1,replace = T, prob = c(nprobs[i], 1-nprobs[i]))
})

dat <- cbind(dat,resp)

dat$ages <- as.factor(dat$ages)
dat$races <- as.factor(dat$races)



mod <- glm(resp ~ gend + ages + races + mentals + 
           substances + perscriptions + hosps + rurals,
           data = dat, family = "binomial")
summary(mod)
