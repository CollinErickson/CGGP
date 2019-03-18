# How I would run laGP on Borehole function

# I normalize both because I was getting errors with darg or garg
# when using a function that was on a small scale.

d <- 8
f <- TestFunctions::borehole
n <- 1000
x <- lhs::maximinLHS(n, d)

np <- 150
xtest <- lhs::maximinLHS(np, d)
ytest <- f(xtest)


# 1. Using full GP model. Hard part is setting up d and g
y <- f(x)
sdy <- sd(y)
mny <- mean(y)
y <- (y-mny) / sdy
mod.agp <- laGP::newGPsep(X=x, Z=y,
                          d=laGP::darg(d=list(mle = TRUE, max = 100), X=x)$start,
                          g=laGP::garg(g=list(mle = TRUE), y=y)$start)
laGP::updateGPsep(mod.agp, x, y)
pred <- laGP::predGPsep(mod.agp, xtest, lite=T)
pred_var <- pred$s2 * sdy^2
pred_mean <- pred$mean * sdy + mny
plot(ytest, pred_mean); abline(a=0, b=1, col=2)


# 2. Using local approximate GP model. 
y <- f(x)
sdy <- sd(y)
mny <- mean(y)
y <- (y-mny) / sdy
pred <- laGP::aGPsep(X=x, Z=y, XX=xtest, method="alc")
pred_var <- pred$var * sdy^2
pred_mean <- pred$mean * sdy + mny
plot(ytest, pred_mean); abline(a=0, b=1, col=2)
