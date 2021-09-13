
context("CJS")

test_that("markCJS",{
  data<-simdataCJS(delta_1=1,delta_2=0)$Enc.Mat
  expect_error(markCJS(data,iter=10,burnin=0),NA)
  expect_error(markCJS(data,mod.phi=~age,iter=10,burnin=0,parameters=list(Phi=list(age.bins=c(0,1,4))),right=FALSE),NA)
})

test_that("multimarkCJS",{
  expect_error(multimarkCJS(Enc.Mat=bobcat,data.type="never",iter=10,burnin=0),NA)
  expect_error(multimarkCJS(Enc.Mat=bobcat,mod.phi=~age,iter=10,burnin=0,parameters=list(Phi=list(age.bins=c(0,1,7))),right=FALSE),NA)
  expect_error(getprobsCJS(multimarkCJS(Enc.Mat=bobcat,data.type="never",iter=10,burnin=0)),NA)
  
  setup<-processdata(bobcat)
  expect_error(test.dot<-multimarkCJS(mms=setup,parms="all",iter=10,burnin=0),NA)
  expect_error(multimodelCJS(modlist=list(mod1=test.dot,mod2=test.dot)),NA)
})
