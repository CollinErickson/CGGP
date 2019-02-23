# Check our supp options on redTime
# 6 HandlingSupp options
# A few different redTime objects.
# Can also compare laGP and MRFA on same data.


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
  
  # Get SGGP/grid data
  if (Ngrid == 1319) {
    sgo <- readRDS("../../../Desktop/redTimeData/out_redTimeTest1_SGGP_after_fit_1319.rds")
  } else if (Ngrid == 8039) {
    sgo <- readRDS("../../../Desktop/redTimeData/out_redTimeTest1_SGGP_after_fit_8039.rds")
  } else if (Ngrid == 3119) {
    sgo <- readRDS("../../../Desktop/redTimeData/out_redTimeTest1_SGGP_after_fit_3119.rds")
  } else if (Ngrid == 1319) {
    sgo <- readRDS("../../../Desktop/redTimeData/out_redTimeTest1_SGGP_after_fit_1319.rds")
  } else if (Ngrid == 0) {
    sgo <- list()
  } else {
    stop(paste("Not acceptable Ngrid", Ngrid))
  }
  
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
    sgo$HandlingSuppData <- Supp
    # browser()
    sgo$CorrMat <- SGGP_internal_CorrMatCauchySQ
    sgo$thetaMAP <- rep(0,18)
    sgo$numpara <- 2
    sgo <- SGGPfit(sgo, Y=log(sgo$Y[,outdim]), Xs=Xs, Ys=Ys)
    pred <- SGGPpred(sgo, xtest)
  } else if (package == "aGP") {
    x <- rbind(sgo$design, Xs)
    # browser()
    y <- log(c(sgo$Y[,outdim], Ys))
    if (nrow(x) != Ngrid+Nsupp) {stop("wrong agp x length")}
    # browser()
    if (nrow(x) != length(y)) {stop("agp x and y don't match")}
    pred <- laGP::aGPsep(X=x, Z=y, XX=xtest, method="alc")
  } else {
    stop(paste("Bad package given to run_func:", package))
  }
  # browser()
  out.stats <- SGGP::valstats(c(pred$mean), pred$var, ytest[,outdim])
  out.stats
}





pko <- expand.grid.df(
  rbind(
    expand.grid(
      package="SGGP",
      Supp=c("Ignore", "Only", "Correct","Mixture", "MarginalValidation","FullValidation")),
    data.frame(package="aGP", Supp="NA", stringsAsFactors=F)),
  data.frame(Ngrid=c(1319,3119, 8039)),
             data.frame(Nsupp=9*c(0, 5, 10, 20)),
                        data.frame(outdim=10)
)
pko$package <- as.character(pko$package)
pko$Supp <- as.character(pko$Supp)

r1 <- comparer::ffexp$new(eval_func = function(...)run_func(...),
                    pko=pko,
                    parallel=F,
                    folder_path = NULL
)
r1$run_all(run_order = "random")
