# transition matrix for illness-death model
tmat <- trans.illdeath()
# data in wide format, for transition 1 this is dataset E1 of
# Therneau & Grambsch (2000)
tg <- data.frame(illt=c(1,1,6,6,8,9),ills=c(1,0,1,1,0,1),
                 dt=c(5,1,9,7,8,12),ds=c(1,1,1,1,1,1),
                 x1=c(1,1,1,0,0,0),x2=c(6:1))
# data in long format using msprep
tglong <- msprep(time=c(NA,"illt","dt"),status=c(NA,"ills","ds"),
                 data=tg,keep=c("x1","x2"),trans=tmat)
# expanded covariates
tglong <- expand.covs(tglong,c("x1","x2"))
# Cox model with different covariate
cx <- coxph(Surv(Tstart,Tstop,status)~x1.1+x2.2+strata(trans),
            data=tglong,method="breslow")

# new data, to check whether results are the same for transition 1 as
# those in appendix E.1 of Therneau & Grambsch (2000)
newdata1 <- data.frame(trans=1:3,x1.1=c(0,0,0),x2.2=c(0,1,0),strata=1:3)


msf1 <- msfit(cx,newdata1,trans=tmat, variance = FALSE)
expect_warning(msfit_subjects(cx, newdata1, trans = tmat))
msf_subjects1 <- suppressWarnings(msfit_subjects(cx, newdata1, trans = tmat))
                                  
                                  
#Compare probtrans outputs from msfit and msfit.subjects
pt_msfit <- probtrans(msf1, predt = 0, variance = FALSE)
pt_msfit_subjects <- probtrans_subjects(msf_subjects1, predt = 0)

#We want to obtain the same transition probabilities
expect_equal(pt_msfit[[1]], pt_msfit_subjects[[1]][[1]])
expect_equal(pt_msfit[[2]], pt_msfit_subjects[[1]][[2]])
expect_equal(pt_msfit[[3]], pt_msfit_subjects[[1]][[3]])





