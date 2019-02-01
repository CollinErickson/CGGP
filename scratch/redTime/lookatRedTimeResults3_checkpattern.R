

pred.sgT2_1063 <- SGGPpred(sgT2_1063, xtest)
plot(ytest, pred.sgT2_1063$mean)
SGGPvalplot(sgT2_1063, xtest, ytest)
pairs(cbind(xtest, abs(ytest - pred.sgT2_1063$mean)[,1]))
lm(abs(ytest - pred.sgT2_1063$mean)[,1] ~ xtest)
pairs(cbind(abs(xtest-.5), abs(ytest - pred.sgT2_1063$mean)[,1]))
lm(abs(ytest - pred.sgT2_1063$mean)[,1] ~ abs(xtest-.5))


pred.sgT3_227 <- SGGPpred(sgT3_227, xtest)
pred.sgT3_455 <- SGGPpred(sgT3_455, xtest)
pred.sgT3_1061 <- SGGPpred(sgT3_1061, xtest)
pred.sgT3_1669 <- SGGPpred(sgT3_1669, xtest)

plot(ytest, pred.sgT3_227$mean)
points(ytest, pred.mlegp$mean, col=3)
summary(ytest - pred.sgT3_227$mean)
colMeans(ytest - pred.sgT3_227$mean)
plot(colMeans(ytest - pred.sgT3_227$mean))
points(colMeans(ytest - pred.mlegp$mean), col=3)
colMeansErrs <- rbind(data.frame(CME=colMeans(ytest - pred.sgT3_227$mean), clm=1:100, model="SGGP227"),
                      data.frame(CME=colMeans(ytest - pred.mlegp$mean), clm=1:100, model="mlegp100")
)
colMeansErrs <- rbind(colMeansErrs,
                      data.frame(CME=colMeans(ytest - pred.sgT3_455$mean), clm=1:100, model="SGGP455"),
                      data.frame(CME=colMeans(ytest - pred.sgT3_1061$mean), clm=1:100, model="SGGP1061"),
                      data.frame(CME=colMeans(ytest - pred.sgT3_1669$mean), clm=1:100, model="SGGP1669")
)
ggplot(colMeansErrs, aes(clm, CME, color=model)) + geom_point()
valstats(pred.sgT3_227$me, pred.sgT3_227$v, ytest)
valstats(pred.mlegp$me, pred.mlegp$v, ytest)
# Var1 is point number, Var2 is output dimension
allErrs <- rbind(
  cbind(reshape2::melt(unname(ytest - pred.sgT3_227$mean)), model="SGGP227"),
  cbind(reshape2::melt(unname(ytest - pred.mlegp$mean)), model="mlegp100")
)
ggplot(allErrs) + geom_point(aes(x=Var2, y=value, color=model, group=Var1))
ggplot(allErrs %>% filter(model=="SGGP227")) + geom_point(aes(x=Var2, y=value, color=Var1)) + 
  geom_smooth(aes(x=Var2, y=value, color=Var1), color="orange") + ggtitle("SGGP residuals")
ggplot(allErrs %>% filter(model=="mlegp100")) + geom_point(aes(x=Var2, y=value, color=Var1))
qplot(ytest[,1], pred.sgT3_227$m[,1]) + geom_abline(intercept=0, slope=1, color="pink")
qplot(ytest[1,], pred.sgT3_227$m[1,]) + geom_abline(intercept=0, slope=1, color="pink")
qplot(ytest[1,], ytest[1,]-pred.sgT3_227$m[1,]) + geom_abline(intercept=0, slope=1, color="pink")


sgtmp <- rlang::duplicate(sgT3_227)
sgtmp$Y <- sgtmp$Y[,1]
sgtmp$thetaMAP <- sgtmp$thetaMAP[,1]
sgtmp <- SGGPfit(sgtmp, sgtmp$Y)
sgtmp$thetaMAP
sgT3_227$thetaMAP[,1]
pred.sgtmp <- SGGPpred(sgtmp, xtest)
valstats(pred.sgtmp$me, pred.sgtmp$v, ytest[,1])
plot(ytest[,1], ytest[,1] -pred.sgtmp$me)

if (F) {
  # # mlegp is bad for this
  # mod.mlegp.ours <- mlegp::mlegp(X=sgT3_227$design, Z=sgT3_227$Y[,1:5])
  # # pred.mlegp.our <- predict(mod.mlegp.ours, xtest)
  # pred.mlegp.ours.raw <- lapply(1:mod.mlegp.ours$numGPs,
  #                          function(i)predict(mod.mlegp.ours[[i]], xtest, se.fit=T))
  # pred.mlegp.ours <- list(
  #   mean = do.call(cbind, lapply(pred.mlegp.ours.raw, function(x) x$fit)),
  #   var = do.call(cbind, lapply(pred.mlegp.ours.raw, function(x) x$se^2))
  # )
  # pred.mlegp.ours %>% str
  # pm <- reshape2::melt(pred.mlegp.ours$mean)
  # ggplot(pm, aes(x=Var2, y=value, color=Var1)) + geom_point()
}
if (F) { # Try laGP instead 
  gps <- lapply(1:100, function(i) {
    print(i)
    IGP::IGP(sgT3_227$design, sgT3_227$Y[,i], "laGP")
  })
  pred.gps <- lapply(1:100, function(i) {
    print(i)
    predict(gps[[i]], xtest)
  })
  pred.gps.mean <- NULL
  for (i in 1:100) {
    pred.gps.mean <- rbind(pred.gps.mean, cbind(data.frame(y=pred.gps[[i]], err=ytest[,i]-pred.gps[[i]], inrow=1:1000, outdim=i)))
  }
  pred.gps.mean %>% str
  pred.gps.mean %>% head
  ggplot(pred.gps.mean, aes(x=outdim, y=y, color=inrow)) + geom_point() + ggtitle("laGP predictions")
  ggplot(pred.gps.mean, aes(x=outdim, y=err, color=inrow)) + geom_point() + ggtitle("laGP residuals") + geom_smooth(color="orange")
  ggplot(pred.gps.mean %>% group_by(outdim) %>% summarise(err=mean(err)), aes(x=outdim, y=err)) +
    geom_point() + ggtitle("laGP residuals means")
  colmeansall3 <- rbind(cbind(pred.gps.mean %>% group_by(outdim) %>% summarise(CME=mean(err), model="laGP") %>% rename(clm=outdim)),
                        colMeansErrs
                        )
  ggplot(colmeansall3, aes(x=clm, y=CME, color=model)) + geom_point() + ggtitle("Mean error by column")
}