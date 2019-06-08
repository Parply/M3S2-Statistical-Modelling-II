#1

bp$female <- as.factor(bp$female)
contrasts(bp$female)
summary(bp)
plot(bp[,-3])
bp <- bp[which(bp$dose != 69.4),]
summary(bp)
plot(bp[,-3])
plot(bp[,-3],col=bp$female)
legend(7,4.3,c("male","female"),col=c("red","black"),pch=1)

model0 <- lm(response~.,data = bp)
model1 <- lm(response~dose,data = bp)
model2 <- lm(response~.^2,data = bp)

par(mfrow=c(2,2))
plot(model0)
plot(model1)
plot(model2)


#2


model <- glm(y~.,data=read,family="poisson")
plot(model)

X <- cbind(1,read$x,read$y)
y <- read$y
beta <- c(3,0,0) # initial guess
r <- 3
seqvec <- Vectorize(seq.defalt,vectorize.args = c("from","to"))
#devience function
D <- function(n){ # n represents inputed mu's and p are the probabilities and y is the count 
  a <- (n-y)*log(p) # term 1
  b <- apply(seqvec(from=n+1,to=n+r-1,by=1),2,prod)/apply(seqvec(from=y+1,to=y+r-1,by=1),2,prod) # term 2
  2*r*sum(a+b) # sum and calculate devience
}
nblnglm<-function(X,y,beta,r=3,terminate=1e-8){
  eta <- X%*%beta
  mu <- exp(eta)
  p <- mu/(r+mu)
  oldD <- D(as.vector(exp(as.numeric(X%*%beta))))
  control <- Inf
  while (control > terminate){
    eta <- X%*%beta
    mu <- exp(eta)
    p <- mu/(r+mu)
    detadmu <- 1/mu
    z <- eta + (y-mu)*detadmu
    w <- (1-p)*mu
    lmod <- lm(z~X+0,weights=w)
    beta <- as.vector(lmod$coeff)
    newD <- D(as.vector(exp(as.numeric(X%*%beta))))
    control <- abs(newD-oldD)/(abs(newD)+0.1)
    oldD <- newD
    
  }
  rtrn <- list(lmod,w,beta,newD,p)
  names(rtrn) <- c("lmod","w","beta","newD","p")
  return(rtrn)
}
rtrn <- nblnglm
lmod <- rtrn&lmod
w <- rtrn$w
beta <- rtrn$beta
newD <- rtrn$newD
p <- rtrn$p
#standard error
J <- t(X)%*%diag(as.vector(w))%*%X
inv.J <- solve(J)
beta.sd <- sqrt(as.vector(diag(inv.J)))

#deviance residuals
n <- as.vector(exp(X%*%beta))
a <- (n-y)*log(p)
b <- apply(seqvec(from=n+1,to=n+r-1,by=1),2,prod)/apply(seqvec(from=y+1,to=y+r-1,by=1),2,prod)
d <- sign(y-mu)*sqrt(2*(a+b))
summary(d)

#pvalues
z <- beta/beta.sd
pvalues <- 2*(1-pnorm(abs(z),lower.tail = TRUE))

#AIC
AIC <- -2*sum(log(choose(y+r-1,y)*((1-p)^r)*(p^y))) + 2*length(beta)

plot(lmod)

X <- cbind(1,read$x,read$y,read$x*read$y)
beta <- c(3,0,0,0) # initial guess
eta <- X%*%beta
mu <- exp(eta)
p <- mu/(r+mu)
oldD <- D(as.vector(exp(as.numeric(X%*%beta))))
control <- Inf
while (control > 1e-8){
  eta <- X%*%beta
  mu <- exp(eta)
  p <- mu/(r+mu)
  detadmu <- 1/mu
  z <- eta + (y-mu)*detadmu
  w <- (1-p)*mu
  lmod <- lm(z~X+0,weights=w)
  beta <- as.vector(lmod$coeff)
  newD <- D(as.vector(exp(as.numeric(X%*%beta))))
  control <- abs(newD-oldD)/(abs(newD)+0.1)
  print(newD)
  oldD <- newD
  
}

#standard error
J <- t(X)%*%diag(as.vector(w))%*%X
inv.J <- solve(J)
beta.sd <- sqrt(as.vector(diag(inv.J)))

#deviance residuals
n <- as.vector(exp(X%*%beta))
a <- (n-y)*log(p)
b <- apply(seqvec(from=n+1,to=n+r-1,by=1),2,prod)/apply(seqvec(from=y+1,to=y+r-1,by=1),2,prod)
d <- sign(y-mu)*sqrt(2*(a+b))
summary(d)

#pvalues
z <- beta/beta.sd
pvalues <- 2*(1-pnorm(abs(z),lower.tail = TRUE))

#AIC
AIC <- -2*sum(log(choose(y+r-1,y)*((1-p)^r)*(p^y))) + 2*length(beta)

plot(lmod)




X <- cbind(1,read$x) # create design matrix
y <- read$y # response vector
beta <- c(3,0) # initial guess
r <- 3 # set r
terminate <- 10e-8
eta <- X%*%beta # calculate initial eta
mu <- exp(eta) # apply inverse of log to eta to get mu
p <- mu/(r+mu) # calculate p from mu using its definition
oldD <- D(exp(as.numeric(X%*%beta))) # calculate initial deviance of initial guess
control <- Inf # set control
while (control > terminate){ # loop until deviance criterion is less than terminate
  eta <- X%*%beta # calculate eta
  mu <- exp(eta) # apply inverse log to get mu
  p <- mu/(r+mu) # calculate p from mu
  detadmu <- 1/mu # calculate differential of eta wrt mu
  z <- eta + (y-mu)*detadmu # calculate z
  w <- (1-p)*mu # calculate weights 
  lmod <- lm(z~X+0,weights=w) # fit model with weights to design matrix
  beta <- as.vector(lmod$coeff) # extract new betas
  newD <- D(exp(as.numeric(X%*%beta))) # calculate new deviance
  control <- abs(newD-oldD)/(abs(newD)+0.1) # calculate criterion value
  oldD <- newD # assign as oldD for next loop
  
}
sprintf("Final deviance: %s",round(newD,4)) # print deviance










