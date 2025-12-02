
#
# R script
#
# Markov chain Monte Carlo simulation of model:
#
# 	(n1, n) ~ JC69(t), and
#         t ~ Exponential(lambda)
#
# where n1 = # of changed sites, n = total # of sites and t is the branch lengths
#


jc69.likelihood = function(x, n1, n) {
	p = 0.25 + 0.75*exp(-4/3 * x)
	choose(n,n1) * p^(n-n1) * (1-p)^(n1)
}

prior = function(x,lambda) {
	lambda * exp(-lambda * x)
}

jc69.likeprior = function(x,n1,n,lambda=10) {
	jc69.likelihood(x,n1,n) * prior(x,lambda)
}

find.MLE = function(n1,n,min=0,max=1) {
	optimize(jc69.likelihood,interval=c(min,max),maximum=T,n1=n1,n=n)
}

marginal.likelihood = function(n1,n,lambda) {
	integrate(jc69.likeprior,lower=0,Inf,n1=n1,n=n,lambda=lambda)
}

make.plot = function(sample=NA,max=1) {
	if (is.na(sample)) {
		plot(0,0,ylim=c(0,9),xlim=c(0,max),type="n",ylab="Density",
			xlab="Branch Length")
	} else {
		hist(sample,ylim=c(0,9),probability=T,xlim=c(0,max),ylab="Density",
			xlab="Branch Length")
	}
}

add.likelihood = function(x,n1,n) {
	norm = integrate(jc69.likelihood,0,Inf,n1=n1,n=n)[[1]]
	lines(x,jc69.likelihood(x,n1,n)/norm,lty=2)
}

add.prior = function(x,lambda) {
	lines(x,prior(x,lambda),lty=3)
}

add.posterior = function(x,n1,n,lambda) {
	p.Y = marginal.likelihood(n1,n,lambda)[[1]]
	lines(x,jc69.likeprior(x,n1,n,lambda)/p.Y)
}

propose.sliding.window = function(theta,w=1) {
	abs(runif(1, theta-w/2,theta+w/2))
}

propose.proportion.scale = function(theta,eps=1) {
	theta * exp( eps * runif(1,-0.5,0.5) )
}

hastings.sliding.window = function(theta.star,theta,w=1) {
	1
}

hastings.proportion.scale = function(theta.star,theta,eps=1) {
	theta.star / theta
}

mcmc.sample.window = function(n1,n,lambda,initial=1, length=100,
                              proposal,
                              hastings,tuning=1) {
	theta = c(1:length)
	theta[1] = initial
	success = 0

	for(i in c(1:(length-1))) {
		theta.star = proposal(theta[i],tuning)
		alpha = jc69.likeprior(theta.star,n1,n,lambda) /
		        jc69.likeprior(theta[i],n1,n,lambda) *
			hastings(theta.star,theta[i],tuning)
		if( runif(1) < alpha ) {
			theta[i+1] = theta.star
			success = success + 1
		} else {
			theta[i+1] = theta[i]
		}
	}

	print(cat("Acceptance rate:",(success/(length-1))))
	theta
}

x.seq = seq(from=0,to=1,length=100)
