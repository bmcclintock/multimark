
context("Closed")

test_that("markClosed",{
  expect_error(markClosed(Enc.Mat=simdataClosed(delta_1=1,delta_2=0)$Enc.Mat,iter=10,burnin=0,bin=5),NA)
})

test_that("multimarkClosed",{
  expect_error(multimarkClosed(Enc.Mat=bobcat,data.type="never",iter=10,burnin=0,bin=5),NA)
  expect_error(getprobsClosed(multimarkClosed(Enc.Mat=bobcat,data.type="never",iter=10,burnin=0,bin=5)),NA)
  
  setup<-processdata(bobcat)
  expect_error(test.dot<-multimarkClosed(mms=setup,parms="all",iter=10,burnin=0,bin=5),NA)
  expect_error(multimodelClosed(modlist=list(mod1=test.dot,mod2=test.dot)),NA)
})
