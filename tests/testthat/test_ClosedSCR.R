
context("ClosedSCR")

test_that("markClosedSCR",{
  sim.data<-simdataClosedSCR(delta_1=1,delta_2=0)
  Enc.Mat<-sim.data$Enc.Mat
  trapCoords<-sim.data$spatialInputs$trapCoords
  studyArea<-sim.data$spatialInputs$studyArea
  expect_error(markClosedSCR(Enc.Mat,trapCoords,studyArea,iter=10,burnin=0,bin=5),NA)
})

test_that("multimarkClosedSCR",{
  sim.data<-simdataClosedSCR(N=30,noccas=5,ntraps=4)
  Enc.Mat <- sim.data$Enc.Mat
  trapCoords <- sim.data$spatialInputs$trapCoords
  studyArea <- sim.data$spatialInputs$studyArea
  expect_error(multimarkClosedSCR(Enc.Mat,trapCoords,studyArea,iter=10,burnin=0,bin=5),NA)

  sim.data<-simdataClosedSCR()
  Enc.Mat<-sim.data$Enc.Mat
  trapCoords<-sim.data$spatialInputs$trapCoords
  studyArea<-sim.data$spatialInputs$studyArea
  expect_error(getprobsClosedSCR(multimarkClosedSCR(Enc.Mat,trapCoords,studyArea,iter=10,burnin=0,bin=5)),NA)
  expect_error(getdensityClosedSCR(multimarkClosedSCR(Enc.Mat,trapCoords,studyArea,iter=10,burnin=0,bin=5)),NA)
  
  setup<-processdataSCR(Enc.Mat,trapCoords,studyArea)
  expect_error(test.dot<-multimarkClosedSCR(mms=setup,parms="all",iter=10,burnin=0,bin=5),NA)
  expect_error(multimodelClosedSCR(modlist=list(mod1=test.dot,mod2=test.dot)),NA)
})
