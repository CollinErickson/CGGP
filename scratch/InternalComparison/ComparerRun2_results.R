mod <- lm(score ~ corr + sel.method + factor(batchsize) + append.rimseperpoint + grid_size, e2$outcleandf)
mod

library(dplyr)
lm(score ~ corr + sel.method + factor(batchsize) + append.rimseperpoint + grid_size, e2$outcleandf %>% filter(nallotted==16384))

# Get LM coeffs for each nallotted
sapply(unique(e2$outcleandf$nallotted), function(nn) lm(score ~ corr + sel.method + factor(batchsize) + append.rimseperpoint + grid_size, e2$outcleandf %>% filter(nallotted==nn))$coeff) %>% rbind(nallotted=unique(e2$outcleandf$nallotted))

sapply(unique(e2$outcleandf$nallotted)[5], function(nn) lm(score ~ corr + sel.method + factor(batchsize) + append.rimseperpoint + grid_size, e2$outcleandf %>% filter(nallotted==nn))$coeff) %>% {stripchart(data.frame(t(.)), las=1)}

stripchart(elapsedtime ~ nallotted, e2$outcleandf)
# Some got really bad at last nallotted
stripchart(RMSE ~ nallotted, e2$outcleandf)
# All bad ones came from Gaussian
stripchart(RMSE ~ nallotted + corr, e2$outcleandf, las=1)
ggplot(data=e2$outcleandf, mapping=aes(x=nallotted,y=RMSE)) + geom_point() +
  facet_grid(corr ~ .)
ggplot(data=e2$outcleandf, mapping=aes(x=nallotted,y=RMSE, shape=factor(batchsize), color=sel.method)) + geom_point() +
  facet_grid(corr ~ .) + ggplot2::scale_y_log10()
ggplot(data=e2$outcleandf, mapping=aes(x=nallotted,y=score, shape=factor(batchsize), color=sel.method)) + geom_point() +
  facet_grid(corr ~ .)
ggplot(data=e2$outcleandf, mapping=aes(x=nallotted,y=score, shape=factor(batchsize), color=sel.method)) + geom_point() +
  facet_grid(corr ~ .)
ggplot(data=e2$outcleandf[e2$outcleandf$nallotted!=16384,], mapping=aes(x=nallotted,y=score, shape=factor(batchsize), color=sel.method)) + geom_point() +
  facet_grid(corr ~ .)
ggplot(data=e2$outcleandf, mapping=aes(x=nallotted,y=CRPscore, shape=factor(batchsize), color=sel.method)) + geom_point() +
  facet_grid(corr ~ .)
ggplot(data=e2$outcleandf[e2$outcleandf$nallotted!=16384,], mapping=aes(x=nallotted,y=CRPscore, shape=factor(batchsize), color=sel.method)) + geom_point() +
  facet_grid(corr ~ .)

sapply(unique(e2$outcleandf$nallotted), function(nn) lm(score ~ corr + sel.method + factor(batchsize) + append.rimseperpoint + grid_size, e2$outcleandf %>% filter(nallotted==nn, corr!="gaussian", corr!="m52"))$coeff) %>% rbind(nallotted=unique(e2$outcleandf$nallotted)) %>% round(4)
sapply(unique(e2$outcleandf$nallotted), function(nn) lm(log(RMSE) ~ corr + sel.method + factor(batchsize) + append.rimseperpoint + grid_size, e2$outcleandf %>% filter(nallotted==nn, corr!="gaussian", corr!="m52"))$coeff) %>% rbind(nallotted=unique(e2$outcleandf$nallotted)) %>% round(4)
sapply(unique(e2$outcleandf$nallotted), function(nn) lm(log(elapsedtime,2) ~ corr + sel.method + factor(batchsize) + append.rimseperpoint + grid_size, e2$outcleandf %>% filter(nallotted==nn, corr!="gaussian", corr!="m52"))$coeff) %>% rbind(nallotted=unique(e2$outcleandf$nallotted)) %>% round(4)
