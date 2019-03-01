# Check our supp options on redTime
# 6 HandlingSupp options
# A few different redTime objects.
# Can also compare laGP and MRFA on same data.
expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))

# Load test data
xlhs8039 <- unname(as.matrix(read.csv("../../../Desktop/redTimeData/LHS1L_n8039_s1226_matrix.csv")[,-1]))
ylhs8039 <- log(as.matrix(
  read.csv("../../../Desktop/redTimeData/LHS1L_n8039_s1226_all_output.csv")[,-1]))
xtest <- xlhs8039[1:1000,]
ytest <- ylhs8039[1:1000,]
# Load sup data
xsup <- as.matrix(unname(read.csv("../../../Desktop/redTimeData/LHS1L_n1000_s1225_Xmatrix.csv")[,-1]))
ysup <- log(unname(as.matrix(read.csv("../../../Desktop/redTimeData/LHS1L_n1000_s1225_all_output.csv")[,-1])))




run_func <- function(package, Ngrid, Nsupp, Supp, outdim) {
  # Need this first chunk inside when running parallel
  
  # Load test data
  xlhs8039 <- unname(as.matrix(read.csv("../../../Desktop/redTimeData/LHS1L_n8039_s1226_matrix.csv")[,-1]))
  ylhs8039 <- log(as.matrix(
    read.csv("../../../Desktop/redTimeData/LHS1L_n8039_s1226_all_output.csv")[,-1]))
  xtest <- xlhs8039[1:1000,]
  ytest <- ylhs8039[1:1000,]
  # Load sup data
  xsup <- as.matrix(unname(read.csv("../../../Desktop/redTimeData/LHS1L_n1000_s1225_Xmatrix.csv")[,-1]))
  ysup <- log(unname(as.matrix(read.csv("../../../Desktop/redTimeData/LHS1L_n1000_s1225_all_output.csv")[,-1])))
  
  
  
  # browser()
  # Get SGGP/grid data
  # if (Ngrid == 1319) {
  #   sgo <- readRDS("../../../Desktop/redTimeData/out_redTimeTest1_SGGP_after_fit_1319.rds")
  #   Y <- log(sgo$Y)
  # } else if (Ngrid == 8039) {
  #   sgo <- readRDS("../../../Desktop/redTimeData/out_redTimeTest1_SGGP_after_fit_8039.rds")
  #   Y <- log(sgo$Y)
  # } else if (Ngrid == 3119) {
  #   sgo <- readRDS("../../../Desktop/redTimeData/out_redTimeTest1_SGGP_after_fit_3119.rds")
  #   Y <- log(sgo$Y)
  # } else if (Ngrid == 1319) {
  #   sgo <- readRDS("../../../Desktop/redTimeData/out_redTimeTest1_SGGP_after_fit_1319.rds")
  #   Y <- log(sgo$Y)
  # } else if (Ngrid == 227) {
  #   sgo <- readRDS("../../../Desktop/redTimeData/out_T2_SGGP-227.rds")
  #   Y <- sgo$Y
  # } else if (Ngrid == 455) {
  #   sgo <- readRDS("../../../Desktop/redTimeData/out_T2_SGGP-455.rds")
  #   Y <- sgo$Y
  # } else if (Ngrid == 1063) {
  #   sgo <- readRDS("../../../Desktop/redTimeData/out_T2_SGGP-1063.rds")
  #   Y <- sgo$Y
  if (Ngrid == 227) {
    sgo <- readRDS("../../../Desktop/redTimeData/out_T4_SGGP-227.rds")
    Y <- sgo$Y
  } else if (Ngrid == 455) {
    sgo <- readRDS("../../../Desktop/redTimeData/out_T4_SGGP-455.rds")
    Y <- sgo$Y
  } else if (Ngrid == 1061) {
    sgo <- readRDS("../../../Desktop/redTimeData/out_T4_SGGP-1061.rds")
    Y <- sgo$Y
  } else if (Ngrid == 1517) {
    sgo <- readRDS("../../../Desktop/redTimeData/out_T4_SGGP-1517.rds")
    Y <- sgo$Y
  } else if (Ngrid == 2049) {
    sgo <- readRDS("../../../Desktop/redTimeData/out_T4_SGGP-2049.rds")
    Y <- sgo$Y
  } else if (Ngrid == 3265) {
    sgo <- readRDS("../../../Desktop/redTimeData/out_T4_SGGP-3265.rds")
    Y <- sgo$Y
  } else if (Ngrid == 0) {
    sgo <- list(Y <- numeric(0))
    Y <- numeric(0)
  } else {
    stop(paste("Not acceptable Ngrid", Ngrid))
  }
  if (length(Y) > 0) {Y <- Y[,outdim]}
  
  # Get LHS data
  if (Nsupp > 0) {
    supp_inds <- sample(1:nrow(xsup), Nsupp, replace=F)
    Xs <- xsup[supp_inds,]
    Ys <- ysup[supp_inds, outdim]
  } else {
    Xs <- NULL
    Ys <- NULL
  }
  
  
  # Run fitting and make predictions
  if (package == "SGGP") {
    require("SGGP")
    if (Ngrid==0) {
      sgo <- SGGPcreate(d=9,batchsize=0,Xs=Xs,Ys=Ys)
    } else {
      sgo$HandlingSuppData <- Supp
      sgo$CorrMat <- SGGP_internal_CorrMatCauchySQ
      sgo$thetaMAP <- rep(0,18)
      sgo$numpara <- 2
      sgo <- SGGPfit(sgo, Y=Y, Xs=Xs, Ys=Ys) # Doesn't need HandlingSuppData since it is set to sgo
    }
    pred <- SGGPpred(sgo, xtest)
  } else {
    x <- rbind(sgo$design, Xs)
    # browser()
    # y <- c(if (length(sgo$Y)>0) log(sgo$Y[,outdim]) else numeric(0), Ys) # sup was already log
    y <- c(Y, Ys)
    if (nrow(x) != Ngrid+Nsupp) {stop("wrong agp x length")}
    if (nrow(x) != length(y)) {stop("agp x and y don't match")}
    if (package == "aGP") {
      pred <- laGP::aGPsep(X=x, Z=y, XX=xtest, method="alc", end=min(50,nrow(x)-1))
    } else if (package == "laGP") {
      mod.agp <- laGP::newGPsep(X=x, Z=y, d=laGP::darg(d=list(mle = TRUE, max = 100), X=x)$start,
                                g=laGP::garg(g=list(mle = TRUE), y=y)$start)
      laGP::updateGPsep(mod.agp, x, y)
      pred <- laGP::predGPsep(mod.agp, xtest, lite=T)
      pred$var <- pred$s2
    } else if (package == "MRFA") {
      mrfa.try <- try({
        mod <- MRFA::MRFA_fit(X=x, Y=y, verbose=FALSE)
        pred <- list(mean=predict(mod, xtest, lambda = min(mod$lambda))$y_hat)
      })
      # browser()
      if (inherits(mrfa.try, "try-error")) {
        pred <- list(mean=rep(mean(y), nrow(xtest)))
      }
      pred$var <- rep(NaN, nrow(xtest))
    } else if (package == "mlegp") {
      mod <- mlegp::mlegp(X=x, Z=y)
      # browser()
      pred <- list(mean=predict(mod, xtest))
      pred$var=rep(NaN, nrow(xtest))
      pred
    } else if (package == "SVM") {
      mod <- e1071::svm(x, y)
      pred <- list(mean=predict(mod, xtest))
      pred$var=rep(NaN, nrow(xtest))
    } else {
      stop(paste("Bad package given to run_func:", package))
    }
  }
  # browser()
  out.stats <- SGGP::valstats(c(pred$mean), pred$var, ytest[,outdim])
  out.stats
}





pko <- expand.grid.df(
  rbind(
    data.frame(
      package="SGGP",
      Supp=c("Ignore", "Only", "Correct","Mixture", "MarginalValidation","FullValidation")),
    data.frame(package=c("aGP","laGP","SVM","mlegp","MRFA"), Supp=c("aGP","laGP","SVM","mlegp","MRFA"), stringsAsFactors=F)),
  # data.frame(Ngrid=c(0,227,455,1063,1319,3119, 8039)),
  data.frame(Ngrid=c(0,227,455,1061,1517,2049, 3265)),
  data.frame(Nsupp=9*c(0, 5, 10, 20)),
  data.frame(outdim=10)
)
pko$package <- as.character(pko$package)
pko$Supp <- as.character(pko$Supp)
# Remove when Ngrid and Nsupp both 0
pko <- pko[!(pko$Ngrid==0 & pko$Nsupp==0),]
# Remove big laGP
pko <- pko[!(pko$package=="laGP" & (pko$Ngrid+pko$Nsupp)>600),]
# Remove medium mlegp
pko <- pko[!(pko$package=="mlegp" & (pko$Ngrid+pko$Nsupp)>300),]
rownames(pko) <- NULL

r1 <- comparer::ffexp$new(
  # eval_func = function(...)run_func(...),
  eval_func = run_func,
  pko=pko,
  parallel=F,
  folder_path = "./scratch/redTime/redTimeSuppExp2/redTimeSuppExp2_od10/"
)
r1$run_all(run_order = "random")

# saveRDS(r1, "./scratch/redTime/redTimeSuppExp/redTimeSuppExp1_od10_object.rds")
# saveRDS(r1, "./scratch/redTime/redTimeSuppExp/redTimeSuppExp1_od50_object.rds")
# saveRDS(r1, "./scratch/redTime/redTimeSuppExp2/redTimeSuppExp2_od50_object.rds")
# saveRDS(r1, "./scratch/redTime/redTimeSuppExp2/redTimeSuppExp2_od10_object.rds")

if (F) {
  r1$parallel <- T
  r1$parallel_cores <- 4
  r1$save_output <- T
  r1$recover_parallel_temp_save(delete_after=F)
  r1$run_all(parallel_temp_save = T, delete_parallel_temp_save_after = F, write_start_files = T, write_error_files = T)
}

if (F) {
  r1$plot_run_times()
  # plyr::dlply(r1$outcleandf, "d")
  require('ggplot2')
  ecdf <- r1$outcleandf[r1$completed_runs,]# & r1$outcleandf$package!="MRFA",]
  ggplot(data=ecdf, mapping=aes(Ngrid, RMSE, color=as.factor(Nsupp))) + geom_point(size=2) + facet_grid(. ~ package+Supp, scales="free_y") + scale_y_log10()
  ggplot(data=ecdf[ecdf$package!="aGP",], mapping=aes(Ngrid, RMSE, color=as.factor(Nsupp))) + geom_point(size=2) + facet_grid(. ~ package+Supp, scales="free_y") + scale_y_log10()
  ggplot(data=ecdf[!(ecdf$package %in% c("mlegp","MRFA","SVM")),], mapping=aes(Ngrid, score, color=as.factor(Nsupp))) + geom_point(size=2) + facet_grid(. ~ package+Supp, scales="free_y")
  ggplot(data=ecdf, mapping=aes(n, CRPscore)) + geom_point() + facet_grid(f ~ package, scales="free_y") + scale_y_log10()
  ggplot(data=ecdf, mapping=aes(Ngrid, runtime)) + geom_point() + facet_grid(. ~ package, scales="free_y") + scale_y_log10()
  # saveRDS(excomp, "./scratch/ExternalComparison/ExComp1_completed.rds")
}
