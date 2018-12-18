test_that("pw is exact", {
  SG <- SGcreate(d=3, batchsize=50)
  y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
  # Make logtheta pretty small to avoid singularity
  logtheta <- c(-2,-2.1,-2.2)
  pw <- calculate_pw(SG=SG, y=y, logtheta=logtheta)
  
  n <- nrow(SG$design)
  R <- outer(1:n, 1:n,
             Vectorize(function(i,j) {
               CorrMatMatern32(SG$design[i,1],
                               SG$design[j,1],
                               logtheta=logtheta[1])*
                 CorrMatMatern32(SG$design[i,2],
                                 SG$design[j,2],
                                 logtheta=logtheta[2])*
                 CorrMatMatern32(SG$design[i,3],
                                 SG$design[j,3],
                                 logtheta=logtheta[3])
             })
  )
  Rinvy <- solve(R, y)
  
  # summary(pw-Rinvy)
  expect_equal(pw, Rinvy)
})

test_that("dpw matches numerical derivative, C and R versions match", {
  SG <- SGcreate(d=3, batchsize=50)
  y <- apply(SG$design, 1, function(x){x[1]+x[2]^2+rnorm(1,0,.01)})
  # Make logtheta pretty small to avoid singularity
  logtheta <- c(-2,-2.1,-2.2)
  
  # Check that R and C give same result
  pwR <- calculate_pw(SG=SG, y=y, logtheta=logtheta, useC=F)
  pwC <- calculate_pw(SG=SG, y=y, logtheta=logtheta, useC=T)
  pw_dpwR <- calculate_pw_and_dpw(SG=SG, y=y, logtheta=logtheta, useC=F)
  pw_dpwC <- calculate_pw_and_dpw(SG=SG, y=y, logtheta=logtheta, useC=T)
  expect_equal(pwR, pwC)
  expect_equal(pwR, pw_dpwR$pw)
  expect_equal(pw_dpwR$pw, pw_dpwC$pw)
  expect_equal(pw_dpwR$dpw, pw_dpwC$dpw)
  
  # Now just make sure that R is correct
  dpw <- pw_dpwR$dpw
  
  eps <- 1e-5
  for (i in 1:3) {
    lteps <- c(0,0,0)
    lteps[i] <- eps/2
    expect_equal(
      (calculate_pw(SG=SG, y=y, logtheta=logtheta+lteps)-calculate_pw(SG=SG, y=y, logtheta=logtheta-lteps))/eps,
      dpw[,i],
      tol=1e-4,
      info=paste("loop is", i, "useC is", useC)
    )
  }
  
  
})

