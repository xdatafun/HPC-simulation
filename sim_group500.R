
# lapply(c("runjags","coda","superdiag","mirt"),library, character.only=TRUE)

library(doParallel)
registerDoParallel(cores=16) 

N.group <- 500 # participants per group
G <- 2
N.total <- G*N.group # total number of people
K = 20 # number of items
# Baker 1998 a~unif(.2,1) d~normal(0,1) 
#reps <- 25

write.table(matrix(c(rep(c(sapply(seq(1:20),function(x) paste("a",x,sep="")),
                sapply(seq(1:20),function(x) paste("d",x,sep="")),
                sapply(seq(1:20),function(x) paste("S_X2_",x,sep="")),
                'BIC','LL'),2),
                c(sapply(seq(1:20),function(x) paste("a",x,sep="")),
                  sapply(seq(1:20),function(x) paste("d",x,sep=""))),
                c('rmse_bayes.a','rmse_bayes.dg1','rmse_bayes.dg2',
                  'bias_bayes.a','bias_bayes.dg1','bias_bayes.dg2',
                  'rmse_ml.ag1','rmse_ml.dg1','rmse_ml.ag2','rmse_ml.dg2',
                  'bias_ag1','bias_dg1','bias_ag2','bias_dg2')),
                nrow=1),
				paste("sim_item20_group",N.group,"-baker-7-2016.csv",sep=""),
				col.name=F,row.name=F,sep=",")
				

foreach(rep=1:24,
        .packages=c('mirt','runjags','superdiag','coda')) %dopar% {

# functions on bias and rmse
rmse <- function(pred,actual){
  error <- actual - pred
  sqrt(mean(error^2))
}

bias <- function(pred,actual){
  error <- actual - pred
  mean(error)
}



diff <- 0
infile <- paste("sim_item20_group",N.group,"-baker-7-2016.csv",sep="")

finals <- matrix(0,nrow=1,ncol=(8*K+4+14))


  ### data generation 

  # a * (theta + d)
  
  ####### MI data 
  a <- matrix(runif(K,0.2,1),ncol=1) # item discrimination 
  d <- matrix(rnorm(K,0,1),ncol=1) # item difficulty 

  Theta <- matrix(rnorm(N.group)) # participants' latent score
  items <- rep('dich', K) # item type
  data1 <- simdata(a, d, N.group, items, Theta=Theta) # generate dichotomous response data
  data2 <- simdata(a, (d+diff), N.group, items, Theta=Theta) 
  scar.data <- cbind(rbind(data1,data2),group=gl(2,N.group)) # add grouping
  
  #------------------------------------------------------------------
  modelstring = "
  model {
  for ( g in 1:G ) {
  for ( i in 1:N ) {
  for ( j in 1:J ) {
  # Likelihood:
  Y[i,j,g] ~ dbern(p[i,j,g])
  #  p[i,j,g] <- exp(a[j,g]*(theta[i,g]+d[j,g]))/(1+exp(a[j,g]*(theta[i,g]+d[j,g]))) # different a & d
  #  p[i,j,g] <- exp(a[j]*(theta[i,g]+d[j,g]))/(1+exp(a[j]*(theta[i,g]+d[j,g]))) # same a, different d
    p[i,j,g] <- exp(a[j]*(theta[i,g]+d[j]))/(1+exp(a[j]*(theta[i,g]+d[j]))) # same a, same d
  }
  # Prior on theta
  theta[i,g] ~ dnorm(0,1)
  } 
  }  
  
  # Priors on item parameters--same a, same d
  for ( j in 1:J) {
    a[j] ~ dnorm( 1, 1)T(0,)
    d[j] ~ dnorm( 0, 1)
  }
  #data# J, N, G, Y
  #monitor# a, d
  }
  
  " # close quote for modelstring
  
  
  #------------------------------------------------------------------------------
  # THE DATA.
  
  
  # Specify the data in a form that is compatible with BRugs model, as a list:
  J <- K
  N <- N.group
  data <- scar.data
  Y <- array(0, dim=c(N,J,G))
  
  for (g in 1:G){ # to re-shape Y into 3 indices
    Y[1:N,1:J,g]<-data[data[,(J+1)]==g,1:J]
  }
  
  
  #------------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # to generate multiple random starting values
  
  n.chains <- 2
  # initials for same a & d  
  initlist <- replicate(n.chains,list(a=rlnorm(J,0,1), 
                                      d=rnorm(J,0,1)),
                        simplify=FALSE)
  #------------------------------------------------------------------------------
  # RUN THE CHAINS.
  
  # Create, initialize, and adapt the model:
  # results.p  <- autorun.jags( model = modelstring, data = dataList,
  #                          monitor=c("theta", "a", "b"), method="parallel")
  
  results.scar <- autorun.jags( modelstring,
                                n.chains = n.chains,
                                inits=initlist
                      #          max.time='3m'
  ) 
  
  ### diag suggested by Jeff Gill ###
  full.out.list <- as.mcmc.list(results.scar$mcmc)
  sink("superdiag.modes.txt")
  superdiag(full.out.list,burnin=5000)
  sink()
  
  
  bayes <- summary(results.scar)
  bayes.a <- bayes[1:J,'Median']
  bayes.d <- bayes[(J+1):(2*J),'Median']
  
  ml.g1<-mirt(data[1:N.group,1:J],1,"2PL")
  ml.est.g1 <-sapply(coef(ml.g1)[1:J],function(x) x)
  #out[i,1:10]<-ml.est[1,] # extract a
  #out[i,11:20]<-ml.est[2,] # extract d 
  # ml.fit.g1 <- itemfit(ml.g1,Zh=F,S_X2=T) # extract item fit S_X2
  
  ml.g2<-mirt(data[(N.group+1):N.total,1:J],1,"2PL")
  ml.est.g2 <-sapply(coef(ml.g2)[1:J],function(x) x)
  # ml.fit.g2 <- itemfit(ml.g2,Zh=F,S_X2=T) # extract item fit S_X2
  
  

  
  quality <- c( rmse(bayes.a,a)
                ,rmse(bayes.d,d)
                ,rmse(bayes.d,(d+diff))
                ,bias(bayes.a,a)
                ,bias(bayes.d,d)
                ,bias(bayes.d,(d+diff))
                ,rmse(ml.est.g1[1,],a)
                ,rmse(ml.est.g1[2,],d)
                ,rmse(ml.est.g2[1,],a)
                ,rmse(ml.est.g2[2,],(d+diff))
                ,bias(ml.est.g1[1,],a)
                ,bias(ml.est.g1[2,],d)
                ,bias(ml.est.g2[1,],a)
                ,bias(ml.est.g2[2,],(d+diff)))
  
  names(quality) <- c('rmse_bayes.a','rmse_bayes.dg1','rmse_bayes.dg2',
                      'bias_bayes.a','bias_bayes.dg1','bias_bayes.dg2',
                      'rmse_ml.ag1','rmse_ml.dg1','rmse_ml.ag2','rmse_ml.dg2',
                      'bias_ag1','bias_dg1','bias_ag2','bias_dg2')
  
  #------------------------------------------------------------------------------
  # EXAMINE THE RESULTS.
  # to export a parameter and take a closer look 
  
  iter <- 200
  mcmcresults <- as.mcmc.list(results.scar, vars=c("a","d"))
  pp.par1 <- as.matrix(tail(mcmcresults[1],(iter-1)))
  pp.par2 <- as.matrix(tail(mcmcresults[2],(iter-1)))
  pp.par <- (pp.par1+pp.par2)/2
  
  
  colnames <- c(sapply(seq(1:J),function(x) paste("a",x,sep="")),
                sapply(seq(1:J),function(x) paste("d",x,sep="")),
                sapply(seq(1:J),function(x) paste("S_X2_",x,sep="")),
                'BIC','LL'
  ) # pc means proportion of correct; num -- number of people in raw score interval
  pp <- matrix(0,nrow=iter,ncol=length(colnames))
  colnames(pp) <- colnames
  # set.seed(123)
  for (i in 1:iter){ # group 1 replications 
    a.par <- matrix(pp.par[i,1:J],ncol=1) # one draw of a
    d.par <- matrix(pp.par[i,(J+1):(2*J)],ncol=1) # one draw of d
    temp <- simdata(a.par, d.par, N.total, itemtype='dich') # one draw of par sim one replication
    mod <- mirt(temp,1,"2PL") # fit the replication data
    est <-sapply(coef(mod)[1:J],function(x) x)
    pp[i,1:J]<-est[1,] # extract a
    pp[i,(J+1):(2*J)]<-est[2,] # extract d 
    # fit <- itemfit(mod,Zh=F,S_X2=T) # extract item fit S_X2
    pp[i,(2*J+1):(3*J)]<-itemfit(mod)[,'S_X2'] 
    pp[i,(3*J+1)] <- mod@Fit$BIC
    pp[i,(3*J+2)] <- mod@Fit$logLik
  }# i
  
  ml.org <- rep(0,dim(pp)[2])
  names(ml.org) <-colnames  
  ml.org[1:J] <- ml.est.g1[1,]
  ml.org[(J+1):(2*J)] <- ml.est.g1[2,]
  ml.org[(2*J+1):(3*J)] <- itemfit(ml.g1)[,'S_X2']
  ml.org[(3*J+1)] <- ml.g1@Fit$BIC
  ml.org[(3*J+2)] <- ml.g1@Fit$logLik

  
  out1 <- rep(0,dim(pp)[2])
  names(out1) <-colnames
  for (name in colnames(pp)){
    out1[name] <- sum(pp[,name]>ml.org[name])/iter
  }
  
  ml.org <- rep(0,dim(pp)[2])
  names(ml.org) <-colnames  
  ml.org[1:J] <- ml.est.g2[1,]
  ml.org[(J+1):(2*J)] <- ml.est.g2[2,]
  ml.org[(2*J+1):(3*J)] <- itemfit(ml.g2)[,'S_X2']
  ml.org[(3*J+1)] <- ml.g2@Fit$BIC
  ml.org[(3*J+2)] <- ml.g2@Fit$logLik
  
  
  out2 <- rep(0,dim(pp)[2])
  names(out2) <-colnames
  for (name in colnames(pp)){
    out2[name] <- sum(pp[,name]>ml.org[name])/iter
  }
  
  out3 <- rep(0,2*J)
  names(out3) <- colnames[1:(2*J)]
  bayes <- c(bayes.a,bayes.d)
  names(bayes) <- colnames[1:(2*J)]
  for (name in names(out3)){
    out3[name] <- sum(pp[,name]>bayes[name])/iter
  }
  
  finals[1,]<-c(out1,out2,out3,quality)
  


# boxplot(finals)

write.table(finals,infile,col.name=F,row.name=F,sep=",",append=T)

}
