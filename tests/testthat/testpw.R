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

test_that("dpw matches numerical derivative", {
  pw_dpw <- calculate_pw_and_dpw(SG=SG, y=y, logtheta=logtheta)
  dpw <- pw_dpw$dpw
  numDeriv::grad(function(x)calculate_pw(SG=SG, y=y, logtheta=x), x=logtheta)
  eps <- 1e-5
  expect_equal(
    (calculate_pw(SG=SG, y=y, logtheta=logtheta+c(eps/2,0,0))-calculate_pw(SG=SG, y=y, logtheta=logtheta-c(eps/2,0,0)))/eps,
    dpw[,1]
  )
})