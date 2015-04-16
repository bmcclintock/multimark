library(multimark)
setup <- processdata(bobcat)

iter <- 2200000
burnin <- 200000
adapt <- 100000
thin <- 20

miter <- 300000

bobcat.dot.same <- multimarkClosed(mms=setup,mod.p=~1,mod.delta=~1,parms="all",nchains=2,iter=iter,burnin=burnin,adapt=adapt,thin=thin,printlog=TRUE)
bobcat.time.same <- multimarkClosed(mms=setup,mod.p=~time,mod.delta=~1,parms="all",nchains=2,iter=iter,burnin=burnin,adapt=adapt,thin=thin,printlog=TRUE)
bobcat.c.same <- multimarkClosed(mms=setup,mod.p=~c,mod.delta=~1,parms="all",nchains=2,iter=iter,burnin=burnin,adapt=adapt,thin=thin,printlog=TRUE)
bobcat.h.same <- multimarkClosed(mms=setup,mod.p=~h,mod.delta=~1,parms="all",nchains=2,iter=iter,burnin=burnin,adapt=adapt,thin=thin,printlog=TRUE)
bobcat.c.time.same <- multimarkClosed(mms=setup,mod.p=~c+time,mod.delta=~1,parms="all",nchains=2,iter=iter,burnin=burnin,adapt=adapt,thin=thin,printlog=TRUE)
bobcat.c.h.same <- multimarkClosed(mms=setup,mod.p=~c+h,mod.delta=~1,parms="all",nchains=2,iter=iter,burnin=burnin,adapt=adapt,thin=thin,printlog=TRUE)
bobcat.time.h.same <- multimarkClosed(mms=setup,mod.p=~time+h,mod.delta=~1,parms="all",nchains=2,iter=iter,burnin=burnin,adapt=adapt,thin=thin,printlog=TRUE)
bobcat.c.time.h.same <- multimarkClosed(mms=setup,mod.p=~c+time+h,mod.delta=~1,parms="all",nchains=2,iter=iter,burnin=burnin,adapt=adapt,thin=thin,printlog=TRUE)
bobcat.dot <- multimarkClosed(mms=setup,mod.p=~1,parms="all",nchains=2,iter=iter,burnin=burnin,adapt=adapt,thin=thin,printlog=TRUE)
bobcat.time <- multimarkClosed(mms=setup,mod.p=~time,parms="all",nchains=2,iter=iter,burnin=burnin,adapt=adapt,thin=thin,printlog=TRUE)
bobcat.c <- multimarkClosed(mms=setup,mod.p=~c,parms="all",nchains=2,iter=iter,burnin=burnin,adapt=adapt,thin=thin,printlog=TRUE)
bobcat.h <- multimarkClosed(mms=setup,mod.p=~h,parms="all",nchains=2,iter=iter,burnin=burnin,adapt=adapt,thin=thin,printlog=TRUE)
bobcat.c.time <- multimarkClosed(mms=setup,mod.p=~c+time,parms="all",nchains=2,iter=iter,burnin=burnin,adapt=adapt,thin=thin,printlog=TRUE)
bobcat.c.h <- multimarkClosed(mms=setup,mod.p=~c+h,parms="all",nchains=2,iter=iter,burnin=burnin,adapt=adapt,thin=thin,printlog=TRUE)
bobcat.time.h <- multimarkClosed(mms=setup,mod.p=~time+h,parms="all",nchains=2,iter=iter,burnin=burnin,adapt=adapt,thin=thin,printlog=TRUE)
bobcat.c.time.h <- multimarkClosed(mms=setup,mod.p=~c+time+h,parms="all",nchains=2,iter=iter,burnin=burnin,adapt=adapt,thin=thin,printlog=TRUE)

modlist<-list(mod1=bobcat.dot.same,mod2=bobcat.time.same,mod3=bobcat.c.same,mod4=bobcat.h.same,mod5=bobcat.c.time.same,mod6=bobcat.c.h.same,mod7=bobcat.time.h.same,mod8=bobcat.c.time.h.same,mod9=bobcat.dot,mod10=bobcat.time,mod11=bobcat.c,mod12=bobcat.h,mod13=bobcat.c.time,mod14=bobcat.c.h,mod15=bobcat.time.h,mod16=bobcat.c.time.h)

stime=proc.time()
bobcat.M.all<-multimodelClosed(mms=setup,modlist=modlist,miter=miter,monparms=c("N","p","c"))
etime=proc.time()-stime