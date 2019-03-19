# Look at lookatRedTimeResults_Big1 to load needed objects.

t50 <- CGGPfit(rt.sggp.2695, Y=rt.sggp.2695$Y[,rep(50, 100)], Xs=rt.sggp.2695$Xs, Ys=rt.sggp.2695$Ys[,rep(50,100)], separateoutputparameterdimensions = F)
t50 %>% CGGPplotcorr()
t50 %>% CGGPplotprojection(outdims = c(10,50,90))
CGGPappend(t50, 500)$design_unevaluated %>% pairs
t50_1 <- CGGPfit(rt.sggp.2695, Y=rt.sggp.2695$Y[,50], Xs=rt.sggp.2695$Xs, Ys=rt.sggp.2695$Ys[,50], separateoutputparameterdimensions = F)
t100_1 <- CGGPfit(rt.sggp.2695, Y=rt.sggp.2695$Y[,100], Xs=rt.sggp.2695$Xs, Ys=rt.sggp.2695$Ys[,100], separateoutputparameterdimensions = F)


CGGPappend(rt.sggp.2695, 500, selectionmethod = "UCB")$design_unevaluated %>% pairs(main="UCB")
CGGPappend(rt.sggp.2695, 500, selectionmethod = "TS")$design_unevaluated %>% pairs(main="TS")
CGGPappend(rt.sggp.2695, 500, selectionmethod = "Greedy")$design_unevaluated %>% pairs(main="Greedy")


CGGPappend(t50_1, 500, selectionmethod = "UCB")$design_unevaluated %>% pairs(main="UCB, only out dim 50")
CGGPappend(t50_1, 500, selectionmethod = "TS")$design_unevaluated %>% pairs(main="TS, only out dim 50")
CGGPappend(t50_1, 500, selectionmethod = "Greedy")$design_unevaluated %>% pairs(main="Greedy, only out dim 50")


CGGPappend(t100_1, 500, selectionmethod = "UCB")$design_unevaluated %>% pairs(main="UCB, only out dim 100")
CGGPappend(t100_1, 500, selectionmethod = "TS")$design_unevaluated %>% pairs(main="TS, only out dim 100")
CGGPappend(t100_1, 500, selectionmethod = "Greedy")$design_unevaluated %>% pairs(main="Greedy, only out dim 100")


CGGPappend(rt.sggp.2695, 500, selectionmethod = "UCB")
CGGPappend(CGGPfit(rt.sggp.2695, Y=rt.sggp.2695$Y,
                   Xs=rt.sggp.2695$Xs, Ys=rt.sggp.2695$Ys,
                   separateoutputparameterdimensions = T),
           500, selectionmethod = "UCB") %>% pairs


ods <- c(1,2);CGGPappend(CGGPfit(rt.sggp.2695, Y=rt.sggp.2695$Y[,ods],
                                 Xs=rt.sggp.2695$Xs, Ys=rt.sggp.2695$Ys[,ods],
                                 set_thetaMAP_to=cbind(rt.sggp.2695$thetaMAP,rt.sggp.2695$thetaMAP),
                                 separateoutputparameterdimensions = T),
                         500, selectionmethod = "UCB")$design_unevaluated %>% pairs(main="1,2 UCB")
ods <- c(1,2);CGGPappend(CGGPfit(rt.sggp.2695, Y=rt.sggp.2695$Y[,ods],
                                 Xs=rt.sggp.2695$Xs, Ys=rt.sggp.2695$Ys[,ods],
                                 set_thetaMAP_to=rt.sggp.2695$thetaMAP,
                                 separateoutputparameterdimensions = F),
                         500, selectionmethod = "UCB")$design_unevaluated %>% pairs(main="1,2 UCB")

ods <- c(1,2);ml <- CGGPfit(rt.sggp.2695, Y=rt.sggp.2695$Y[,ods],
              Xs=rt.sggp.2695$Xs, Ys=rt.sggp.2695$Ys[,ods],
              set_thetaMAP_to=rt.sggp.2695$thetaMAP,
              separateoutputparameterdimensions = F)
CGGPappend(ml, 500, selectionmethod = "Greedy")$design_unevaluated %>% pairs(main="1,2 UCB")

mlrep <- ml
mlrep$thetaMAP <- rep(ml$thetaMAP[1:2], 9)
mlrep$thetaPostSamples <- do.call(rbind, lapply(1:9, function(i)ml$thetaPostSamples[1:2,]))
CGGPappend(ml, 500, selectionmethod = "UCB")$design_unevaluated %>% pairs(main="1,2 UCB")

mat3 <- CGGPfit(rt.sggp.2695, Y=rt.sggp.2695$Y[,ods],
                Xs=rt.sggp.2695$Xs, Ys=rt.sggp.2695$Ys[,ods],
                corr="m32",
                separateoutputparameterdimensions = F)
mat3 %>% CGGPplottheta()
CGGPappend(mat3, 500, selectionmethod = "Greedy")$design_unevaluated %>% pairs(main="1,2 Greedy")
CGGPappend(mat3, 500, selectionmethod = "UCB")$design_unevaluated %>% pairs(main="1,2 UCB")


cau <- CGGPfit(rt.sggp.2695, Y=rt.sggp.2695$Y[,ods],
                Xs=rt.sggp.2695$Xs, Ys=rt.sggp.2695$Ys[,ods],
                corr="cauchy",
                separateoutputparameterdimensions = F)
cau %>% CGGPplottheta()
CGGPappend(cau, 500, selectionmethod = "Greedy")$design_unevaluated %>% pairs(main="1,2 Greedy")
CGGPappend(cau, 500, selectionmethod = "UCB")$design_unevaluated %>% pairs(main="1,2 UCB")

CGGP_internal_neglogpost(ml$thetaMAP, ml, y=ml$y, Xs=ml$Xs, ys=ml$ys)
CGGP_internal_neglogpost(mat3$thetaMAP, mat3, y=mat3$y, Xs=mat3$Xs, ys=mat3$ys)
CGGP_internal_neglogpost(cau$thetaMAP, cau, y=cau$y, Xs=cau$Xs, ys=cau$ys)
nlp_ml <- sapply(1:ml$numPostSamples, function(i)CGGP_internal_neglogpost(ml$thetaPostSamples[,i], ml, y=ml$y, Xs=ml$Xs, ys=ml$ys))
nlp_mat3 <- sapply(1:mat3$numPostSamples, function(i)CGGP_internal_neglogpost(mat3$thetaPostSamples[,i], mat3, y=mat3$y, Xs=mat3$Xs, ys=mat3$ys))
nlp_cau <- sapply(1:cau$numPostSamples, function(i)CGGP_internal_neglogpost(cau$thetaPostSamples[,i], cau, y=cau$y, Xs=cau$Xs, ys=cau$ys))
nlp_ml %>% is.infinite %>% table
nlp_mat3 %>% is.infinite %>% table
nlp_cau %>% is.infinite %>% table

nlp_all <- sapply(1:rt.sggp.2695$numPostSamples, function(i)CGGP_internal_neglogpost(rt.sggp.2695$thetaPostSamples[,i], rt.sggp.2695, y=rt.sggp.2695$y, Xs=rt.sggp.2695$Xs, ys=rt.sggp.2695$ys))



ods <- c(1)
ml1 <- CGGPfit(rt.sggp.2695, Y=rt.sggp.2695$Y[,ods],
              Xs=rt.sggp.2695$Xs, Ys=rt.sggp.2695$Ys[,ods],
              separateoutputparameterdimensions = F)
CGGPappend(ml1, 500, selectionmethod = "Greedy")$design_unevaluated %>% pairs(main="1,2 Greedy")
CGGPappend(ml1, 500, selectionmethod = "UCB")$design_unevaluated %>% pairs(main="1,2 UCB")
nlp_ml1 <- sapply(1:ml1$numPostSamples, function(i)CGGP_internal_neglogpost(ml1$thetaPostSamples[,i], ml1, y=ml1$y, Xs=ml1$Xs, ys=ml1$ys))
nlp_ml1 %>% sort

ods <- c(1, 1)
ml11 <- CGGPfit(rt.sggp.2695, Y=rt.sggp.2695$Y[,ods],
               Xs=rt.sggp.2695$Xs, Ys=rt.sggp.2695$Ys[,ods],
               separateoutputparameterdimensions = F)
CGGPappend(ml11, 500, selectionmethod = "Greedy")$design_unevaluated %>% pairs(main="1,2 Greedy")
CGGPappend(ml11, 500, selectionmethod = "UCB")$design_unevaluated %>% pairs(main="1,2 UCB")
ml11 %>% CGGPplottheta()
nlp_ml11 <- sapply(1:ml11$numPostSamples, function(i)CGGP_internal_neglogpost(ml11$thetaPostSamples[,i], ml11, y=ml11$y, Xs=ml11$Xs, ys=ml11$ys))
nlp_ml11 %>% sort


# Check that grad is still correct
to <- sapply(1:ml$numPostSamples, function(i)CGGP_internal_gneglogpost(ml$thetaPostSamples[,i], ml, y=ml$y, Xs=ml$Xs, ys=ml$ys))
to2 <- sapply(1:ml$numPostSamples,
              function(i) {browser()
                numDeriv::grad(CGGP_internal_neglogpost,
                               x=ml$thetaPostSamples[,i], CGGP=ml, y=ml$y, Xs=ml$Xs, ys=ml$ys)
                }
              )



# Is it fixed?
rt.sggp.2695.refit <- CGGPfit(rt.sggp.2695, Y=rt.sggp.2695$Y,
                              Xs=rt.sggp.2695$Xs, Ys=rt.sggp.2695$Ys,
                              separateoutputparameterdimensions = F)
rt.sggp.2695 %>% CGGPplottheta()
rt.sggp.2695.refit %>% CGGPplottheta()
CGGPappend(rt.sggp.2695,       500, "UCB")$design_unevaluated %>% pairs(main="Old w/UCB")
CGGPappend(rt.sggp.2695.refit, 500, "UCB")$design_unevaluated %>% pairs(main="New w/UCB")
# Check if preds changed
CGGPvalstats(rt.sggp.2695      , x1000, y1000, bydim = F)
CGGPvalstats(rt.sggp.2695.refit, x1000, y1000, bydim = F)
# Check with good data
r2.sggp.8179.refit <- CGGPfit(r2.sggp.8179, Y=r2.sggp.8179$Y,
                              Xs=r2.sggp.8179$Xs, Ys=r2.sggp.8179$Ys,
                              separateoutputparameterdimensions = F)
CGGPvalstats(r2.sggp.8179      , x1000, y1000, bydim = F)
CGGPvalstats(r2.sggp.8179.refit, x1000, y1000, bydim = F)
r2.sggp.8179 %>% CGGPplottheta()
r2.sggp.8179.refit %>% CGGPplottheta()
# Check with smaller data to see how much parameter estimates change
r2.sggp.399.refit <- CGGPfit(r2.sggp.399, Y=r2.sggp.399$Y,
                              Xs=r2.sggp.399$Xs, Ys=r2.sggp.399$Ys,
                              separateoutputparameterdimensions = F)
CGGPvalstats(r2.sggp.399      , x1000, y1000, bydim = F)
CGGPvalstats(r2.sggp.399.refit, x1000, y1000, bydim = F)
r2.sggp.399 %>% CGGPplottheta()
r2.sggp.399.refit %>% CGGPplottheta()
