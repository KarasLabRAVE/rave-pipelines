# Put shared functions in `shared-*.R` so the pipeline is clean

get_default_cores <- function(round = TRUE) {
  re <- (raveio::raveio_getopt("max_worker") + 1) / 2
  if( round ) {
    re <- ceiling(re)
  }
  re
}

ridge <- function(xt, xtp1, lambda, intercept = FALSE, iw) {
  if (!identical(dim(xt), dim(xtp1))) {
    stop("Unmatched dimension")
  }
  nel <- ncol(xt)
  ## Coefficient matrix A
  ## each column is coefficients from a linear regression
  ## formula: xtp1 = xt*A + E
  A <- matrix(0, nel + intercept, nel)
  ## for each electrode
  for (i in seq_len(nel)) {
    y <- xtp1[, i]
    fit <- glmnet::glmnet(xt, y,
                  alpha = 0, lambda = lambda,
                  standardize = FALSE, intercept = intercept
    )
    # fit <- glmnet::cv.glmnet(xt, y,
    #               alpha = 0,
    #               standardize = FALSE, intercept = intercept
    # )

    if (intercept) {
      A[, i] <- as.numeric(coef(fit))
    } else {
      A[, i] <- coef(fit)[-1]
    }
  }
  AEigen <- eigen(A)
  e <- Mod(AEigen$values)
  me <- max(e)
  # if(me >= 1){
  #   print(paste0("Solution in timewindow ",iw," is not stable (",me,")"))
  # } else {
  #   print(paste0("Solution in timewindow ",iw," is stable (",me,")"))
  # }
  A
}

ridgeR2 <- function(xt, xtp1, A) {
  nel <- ncol(xt)
  ypredMat <- predictRidge(xt, A)

  R2 <- rep(0, nel)
  for (i in seq_len(nel)) {
    y <- xtp1[, i]
    ypred <- ypredMat[, i]
    sst <- sum((y - mean(y))^2)
    sse <- sum((ypred - y)^2)
    rsq <- 1 - sse / sst
    R2[i] <- rsq
  }
  R2
}

ridgesearchlambdadichomotomy <- function(xt, xtp1, intercept = FALSE, iw){
  if(!identical(dim(xt),dim(xtp1)))
    stop("Unmatched dimension")
  nel <- ncol(xt)
  ## Coefficient matrix A
  ## each column is coefficients from a linear regression
  ## formula: xtp1 = xt*A + E
  lambdamin <- 0.0001
  lambdamax <- 10


  Aa <- ridge(xt,xtp1,lambda=lambdamin,intercept=F, iw = iw)

  stableam <- TRUE

  nel <- ncol(Aa)
  e <- Mod(eigen(Aa)$values)
  me <- max(e)

  if (me >= 1) {
    stableam <- FALSE
  }

  #print(stableam)

  if(stableam){
    lambdaopt <- lambdamin
    Afin <- Aa
  }else{

    stablea <- stableam
    lambdaa <- lambdamin

    lambdab=lambdamax

    Ab<-ridge(xt,xtp1,lambda=lambdab,intercept=F, iw = iw)

    stableb <- TRUE

    nel <- ncol(Ab)
    e <- Mod(eigen(Ab)$values)
    me <- max(e)

    if (me >= 1) {
      stableb <- FALSE
    }

    #print(stableb)
    k <- 0
    while(k<20){
      lambdac <- (lambdaa + lambdab)*0.5

      Ac<-ridge(xt,xtp1,lambda=lambdac,intercept=F, iw = iw)

      stablec <- TRUE

      nel <- ncol(Ac)
      e <- Mod(eigen(Ac)$values)
      me <- max(e)

      if (me >= 1) {
        stablec <- FALSE
      }

      if(!stablec){
        lambdaa <- lambdac
        lambdaopt <- lambdab

      }else{
        lambdab <- lambdac
        lambdaopt <- lambdac
      }
      k <- k+1

      # print("ite")
      # print(k)
      # print(lambdac)
      # print(stablec)
      # print(lambdaopt)
    }
  }

  Afin <- ridge(xt,xtp1,lambda=lambdaopt,intercept=F, iw = iw)

  attr(Afin, "lambdaopt") <- lambdaopt
  Afin
}

fragilityRowNormalized <- function(A, nSearch = 100) {
  ## The adjacency matrix A here is a transpose of the
  ## adjacency matrix in the original paper
  nel <- ncol(A)
  e <- Mod(eigen(A)$values)
  me <- max(e)

  if (me >= 1) {
    #return(0)
  }


  fragcol <- matrix(0, nel, nel)
  fragNorm <- rep(0, nel)
  omvec <- seq(0, 1, length.out = nSearch + 1)[-1]

  b <- c(0, -1)
  ## for each electrode
  for (i in 1:nel) {
    ## indicate which electrode is disturbed (ith)
    ek <- rep(0, nel)
    ek[i] <- 1
    tek <- t(ek)
    minNorm <- 100000
    minPerturbColumn <- NA
    for (k in seq_len(nSearch)) {
      ## imaginary part
      om <- omvec[k]
      ## real part
      sigma <- sqrt(1 - om^2)
      ## target eigenvalue
      lambda <- complex(real = sigma, imaginary = om)
      ## A - (sigma + j* omega)*I
      mat <- A - lambda * diag(nel)
      imat <- t(solve(mat))

      argument <- tek %*% imat
      B <- rbind(Im(argument), Re(argument))
      ## B^T*(B*B^T)^-1*b
      invBtB <- solve(B %*% t(B))
      prov <- t(B) %*% invBtB %*% b

      sigma_hat <- ek %*% t(prov)

      ## validation
      if (FALSE) {
        A2 <- A + sigma_hat
        e2 <- eigen(A2)$values
        closestIndex <- which.min(abs(e2 - lambda))
        e2[closestIndex]
      }

      norm_sigma_hat <- norm(sigma_hat, type = "2")
      if (norm_sigma_hat < minNorm) {
        minPerturbColumn <- prov
        minNorm <- norm(prov, type = "2")
      }
    }

    fragcol[, i] <- minPerturbColumn
    fragNorm[i] <- minNorm
  }

  maxf <- max(fragNorm)
  fragNorm2 <- (maxf - fragNorm) / maxf

  return(fragNorm2)
}

fragilityRow <- function(A, nSearch = 100) {
  ## The adjacency matrix A here is a transpose of the
  ## adjacency matrix in the original paper
  nel <- ncol(A)
  e <- Mod(eigen(A)$values)
  me <- max(e)

  if (me >= 1) {
    #return(0)
  }


  fragcol <- matrix(0, nel, nel)
  fragNorm <- rep(0, nel)
  omvec <- seq(0, 1, length.out = nSearch + 1)[-1]

  b <- c(0, -1)
  ## for each electrode
  for (i in 1:nel) {
    ## indicate which electrode is disturbed (ith)
    ek <- rep(0, nel)
    ek[i] <- 1
    tek <- t(ek)
    minNorm <- 100000
    minPerturbColumn <- NA
    for (k in seq_len(nSearch)) {
      ## imaginary part
      om <- omvec[k]
      ## real part
      sigma <- sqrt(1 - om^2)
      ## target eigenvalue
      lambda <- complex(real = sigma, imaginary = om)
      ## A - (sigma + j* omega)*I
      mat <- A - lambda * diag(nel)
      imat <- t(solve(mat))

      argument <- tek %*% imat
      B <- rbind(Im(argument), Re(argument))
      ## B^T*(B*B^T)^-1*b
      invBtB <- solve(B %*% t(B))
      prov <- t(B) %*% invBtB %*% b

      sigma_hat <- ek %*% t(prov)

      ## validation
      if (FALSE) {
        A2 <- A + sigma_hat
        e2 <- eigen(A2)$values
        closestIndex <- which.min(abs(e2 - lambda))
        e2[closestIndex]
      }

      norm_sigma_hat <- norm(sigma_hat, type = "2")
      if (norm_sigma_hat < minNorm) {
        minPerturbColumn <- prov
        minNorm <- norm(prov, type = "2")
      }
    }

    fragcol[, i] <- minPerturbColumn
    fragNorm[i] <- minNorm
  }

  return(fragNorm)
}

calc_adj_frag <- function(repository, trial_num, t_window, t_step, soz, sozc, lambda = NULL, nSearch = 100, fs_new = NULL) {
  n_tps <- length(repository$voltage$dimnames$Time)
  n_elec <- length(repository$voltage$dimnames$Electrode)

  # Number of steps
  n_steps <- floor((n_tps - t_window) / t_step) + 1

  # slice of data
  arr <- filearray::filearray_load_or_create(
    filebase = tempfile(),
    dimension = c(n_tps, n_elec),
    type = "float", mode = "readwrite", partition_size = 1L,

    # if repository has changed, re-calculate
    repository_signature = repository$signature,
    t_step = t_step, t_window = t_window,
    trial_num = trial_num,

    on_missing = function(arr) {
      arr$set_header("ready", value = FALSE)
    }
  )

  # check if header `ready` is not TRUE
  if(!isTRUE(arr$get_header("ready", FALSE))) {

    loaded_electrodes <- repository$electrode_list
    raveio::lapply_async(repository$voltage$data_list, function(v) {
      e <- dimnames(v)$Electrode
      idx_e <- loaded_electrodes == e

      arr[,idx_e] <- v[, trial_num, 1, drop = TRUE, dimnames = NULL]

      return()
    })
  }

  signalScaling <- 10^floor(log10(max(arr[])))
  arr[] <- arr[]/signalScaling

  if(!is.null(fs_new)) {
    arr <- gsignal::resample(arr[],fs_new,round(repository$sample_rate))
    n_tps <- dim(arr)[1]
    n_steps <- floor((n_tps - t_window) / t_step) + 1
  }

  pt01EpochRaw <- t(arr[])

  goodChannels <- repository$electrode_list
  rownames(pt01EpochRaw) <- repository$electrode_table$Label[repository$electrode_table$Electrode %in% goodChannels]
  sozNames<-repository$electrode_table$Label[repository$electrode_table$Electrode %in% soz]
  sozIndex<-which(repository$electrode_list%in%soz==TRUE)

  ## Add time stamps to the columns
  times <- seq(repository$time_windows[[1]][1], repository$time_windows[[1]][2], length.out=ncol(pt01EpochRaw))
  times_with_sign <- ifelse(times >= 0, paste0("+", times), as.character(times))
  colnames(pt01EpochRaw)<-times_with_sign

  pt01EcoG<-pt01EpochRaw
  attr(pt01EcoG, "sozIndex") <- sozIndex
  attr(pt01EcoG, "sozNames") <- sozNames

  cl <- parallel::makeCluster(4, type = "SOCK")
  doSNOW::registerDoSNOW(cl)

  epoch <- EZFragility::Epoch(pt01EcoG)
  title <- repository$epoch_name
  fragtest<-EZFragility::calcAdjFrag(
    epoch = epoch, window = t_window,
    step = t_step, nSearch = nSearch, parallel = TRUE, progress = TRUE
  )

  ## stop the parallel backend
  parallel::stopCluster(cl)

  return(list(
    voltage = arr[],
    frag = fragtest$frag,
    frag_ranked = fragtest$frag_ranked,
    R2 = fragtest$R2,
    lambdas = fragtest$lambdas
  ))
}

calc_adj_frag_old <- function(repository, trial_num, t_window, t_step, lambda = NULL, nSearch = 100, fs_new = NULL) {

  n_tps <- length(repository$voltage$dimnames$Time)
  n_elec <- length(repository$voltage$dimnames$Electrode)

  # Number of steps
  n_steps <- floor((n_tps - t_window) / t_step) + 1

  # slice of data
  arr <- filearray::filearray_load_or_create(
    filebase = tempfile(),
    dimension = c(n_tps, n_elec),
    type = "float", mode = "readwrite", partition_size = 1L,

    # if repository has changed, re-calculate
    repository_signature = repository$signature,
    t_step = t_step, t_window = t_window,
    trial_num = trial_num,

    on_missing = function(arr) {
      arr$set_header("ready", value = FALSE)
    }
  )

  # check if header `ready` is not TRUE
  if(!isTRUE(arr$get_header("ready", FALSE))) {

    loaded_electrodes <- repository$electrode_list
    raveio::lapply_async(repository$voltage$data_list, function(v) {
      e <- dimnames(v)$Electrode
      idx_e <- loaded_electrodes == e

      arr[,idx_e] <- v[, trial_num, 1, drop = TRUE, dimnames = NULL]

      return()
    })
  }

  signalScaling <- 10^floor(log10(max(arr[])))
  arr[] <- arr[]/signalScaling

  if(!is.null(fs_new)) {
    arr <- gsignal::resample(arr[],fs_new,repository$sample_rate)
    n_tps <- dim(arr)[1]
    n_steps <- floor((n_tps - t_window) / t_step) + 1
  }

  ## create adjacency array (array of adj matrices for each time window)
  ## iw: The index of the window we are going to calculate fragility

  res <- raveio::lapply_async(seq_len(n_steps), function(iw) {
    ## Sample indices for the selected window
    si <- seq_len(t_window-1) + (iw-1)*t_step
    ## measurements at time point t
    xt <- arr[si,]
    ## measurements at time point t plus 1
    xtp1 <- arr[si + 1,]

    ## Coefficient matrix A (adjacency matrix)
    ## each column is coefficients from a linear regression
    ## formula: xtp1 = xt*A + E
    if (!is.null(lambda)) {
      #message(paste0("running with lambda = ", lambda))
      Ai <- ridge(xt, xtp1, intercept = F, lambda = lambda, iw = iw)
    } else {
      #message("running lambda search")
      Ai <- ridgesearchlambdadichomotomy(xt, xtp1, intercept = F, iw = iw)
    }

    R2 <- ridgeR2(xt,xtp1,Ai)

    return(list(Ai = Ai, R2 = R2))
  }, callback = function(iw) {
    sprintf("Generating Adjacency Matrices|Timewindow %s", iw)
  })

  A <- unlist(lapply(res, function(w){
    w$Ai
  }))
  dim(A) <- c(n_elec, n_elec, n_steps)
  dimnames(A) <- list(
    Electrode1 = repository$electrode_list,
    Electrode2 = repository$electrode_list,
    Step = seq_len(n_steps)
  )

  R2 <- unlist(lapply(res, function(w){
    w$R2
  }))
  dim(R2) <- c(n_elec, n_steps)
  dimnames(R2) <- list(
    Electrode = repository$electrode_list,
    Step = seq_len(n_steps)
  )

  lambdas <- rep(NA, length(res))

  for(i in seq_len(n_steps)){
    lambdas[i] <- attr(res[[i]]$Ai,"lambdaopt")
  }

  # calculate fragility
  f <- unlist(raveio::lapply_async(seq_len(n_steps), function(iw){
    fragilityRowNormalized(A[,,iw], nSearch = nSearch)
  }, callback = function(iw) {
    sprintf("Calculating Fragility|Timewindow %s", iw)
  }))
  dim(f) <- c(n_elec, n_steps)
  f_naked <- f
  dimnames(f) <- list(
    Electrode = repository$electrode_list,
    Step = seq_len(n_steps)
  )

  # ranked fragility map
  f_rank <- matrix(rank(f), nrow(f), ncol(f))
  attributes(f_rank) <- attributes(f)
  f_rank <- f_rank/max(f_rank)

  return(list(
    voltage = arr[],
    adj = A,
    frag = f,
    frag_ranked = f_rank,
    R2 = R2,
    lambdas = lambdas
  ))
}

frag_quantile <- function(repository, f, t_window, t_step, soz, sozc, fs_new = NULL){
  n_tps <- length(repository$voltage$dimnames$Time)
  n_elec <- length(repository$voltage$dimnames$Electrode)
  n_steps <- dim(f)[2]
  epoch_time_window <- repository$time_windows[[1]]
  fs <- round(repository$sample_rate,-1)

  if(!is.null(fs_new)){
    fs <- fs_new
  }

  if(any(repository$electrode_table$Label == "NoLabel")) {
    elec_names <- repository$electrode_table$Electrode[match(c(soz,sozc), repository$electrode_table$Electrode)]
    elec_names <- as.character(elec_names)
  } else {
    elec_names <- repository$electrode_table$Label[match(c(soz,sozc), repository$electrode_table$Electrode)]
  }

  # create fragility map with soz electrodes separated from sozc electrodes
  #sozNames <- repository$electrode_table$Label[repository$electrode_table$Electrode %in% soz]
  #sozcNames <- repository$electrode_table$Label[repository$electrode_table$Electrode %in% sozc]
  soz_i <- match(soz,repository$electrode_list)
  sozc_i <- match(sozc,repository$electrode_list)
  fmap <- f[c(soz_i,sozc_i),]
  stimes <- (seq_len(n_steps)-1)*t_step/fs+epoch_time_window[1]

  # raw fragility map
  fplot <- expand.grid(Time = stimes, Electrode = elec_names)
  fplot$Value <- c(t(fmap))

  # create separate heatmaps for soz and sozc for quantile calcs
  hmapsoz <- fmap[soz_i,]
  hmapsozc <- fmap[sozc_i,]

  #f90soz=quantile(hmapsoz, probs=c(0.9))
  #f90sozc=quantile(hmapsozc,probs=c(0.9))
  #interpretabilityratiosoz=f90soz/f90sozc

  quantilematrixsozsozc=matrix(0,20,length(stimes))
  cmeansoz=c(1:length(stimes))*0
  cmeansozc=c(1:length(stimes))*0
  csdsoz=c(1:length(stimes))*0
  csdsozc=c(1:length(stimes))*0

  for(i in 1:length(stimes)){

    colsoz=hmapsoz[,i]
    colsozc=hmapsozc[,i]

    meansoz=mean(colsoz)
    sdsoz=sd(colsoz)
    meansozc=mean(colsozc)
    sdsozc=sd(colsozc)

    cmeansoz[i]=meansoz
    cmeansozc[i]=meansozc
    csdsoz[i]=sdsoz
    csdsozc[i]=sdsozc

    f10colsoz<-quantile(colsoz,probs=c(0.1))
    f20colsoz<-quantile(colsoz,probs=c(0.2))
    f30colsoz<-quantile(colsoz,probs=c(0.3))
    f40colsoz<-quantile(colsoz,probs=c(0.4))
    f50colsoz<-quantile(colsoz,probs=c(0.5))
    f60colsoz<-quantile(colsoz,probs=c(0.6))
    f70colsoz<-quantile(colsoz,probs=c(0.7))
    f80colsoz<-quantile(colsoz,probs=c(0.8))
    f90colsoz<-quantile(colsoz,probs=c(0.9))
    f100colsoz<-quantile(colsoz,probs=c(1.0))

    f10colsozc<-quantile(colsozc,probs=c(0.1))
    f20colsozc<-quantile(colsozc,probs=c(0.2))
    f30colsozc<-quantile(colsozc,probs=c(0.3))
    f40colsozc<-quantile(colsozc,probs=c(0.4))
    f50colsozc<-quantile(colsozc,probs=c(0.5))
    f60colsozc<-quantile(colsozc,probs=c(0.6))
    f70colsozc<-quantile(colsozc,probs=c(0.7))
    f80colsozc<-quantile(colsozc,probs=c(0.8))
    f90colsozc<-quantile(colsozc,probs=c(0.9))
    f100colsozc<-quantile(colsozc,probs=c(1.0))

    quantilematrixsozsozc[1,i]=f10colsoz
    quantilematrixsozsozc[2,i]=f20colsoz
    quantilematrixsozsozc[3,i]=f30colsoz
    quantilematrixsozsozc[4,i]=f40colsoz
    quantilematrixsozsozc[5,i]=f50colsoz
    quantilematrixsozsozc[6,i]=f60colsoz
    quantilematrixsozsozc[7,i]=f70colsoz
    quantilematrixsozsozc[8,i]=f80colsoz
    quantilematrixsozsozc[9,i]=f90colsoz
    quantilematrixsozsozc[10,i]=f100colsoz
    quantilematrixsozsozc[11,i]=f10colsozc
    quantilematrixsozsozc[12,i]=f20colsozc
    quantilematrixsozsozc[13,i]=f30colsozc
    quantilematrixsozsozc[14,i]=f40colsozc
    quantilematrixsozsozc[15,i]=f50colsozc
    quantilematrixsozsozc[16,i]=f60colsozc
    quantilematrixsozsozc[17,i]=f70colsozc
    quantilematrixsozsozc[18,i]=f80colsozc
    quantilematrixsozsozc[19,i]=f90colsozc
    quantilematrixsozsozc[20,i]=f100colsozc

  }

  quantilesname<-c("SOZ(10th)","SOZ(20th)","SOZ(30th)","SOZ(40th)","SOZ(50th)",
                   "SOZ(60th)","SOZ(70th)","SOZ(80th)","SOZ(90th)","SOZ(100th)",
                   "SOZc(10th)","SOZc(20th)","SOZc(30th)","SOZc(40th)","SOZc(50th)",
                   "SOZc(60th)","SOZc(70th)","SOZc(80th)","SOZc(90th)","SOZc(100th)")
  quantileplot<- expand.grid(Time = stimes, Stats=quantilesname)
  quantileplot$Value <- c(t(quantilematrixsozsozc))

  dimnames(quantilematrixsozsozc) <- list(
    Quantile = quantilesname,
    Time = stimes
  )

  return(list(
    fplot = fplot,
    q_matrix = quantilematrixsozsozc,
    q_plot = quantileplot
  ))
}

mean_f_calc <- function(repository, f, soz, sozc) {

  mean_f_soz <- rep(0,dim(f)[2])
  mean_f_sozc <- rep(0,dim(f)[2])
  se_f_soz <- rep(0,dim(f)[2])
  se_f_sozc <- rep(0,dim(f)[2])

  # sozNames <- repository$electrode_table$Label[repository$electrode_table$Electrode %in% soz]
  # sozcNames <- repository$electrode_table$Label[repository$electrode_table$Electrode %in% sozc]
  soz_i <- match(soz,repository$electrode_list)
  sozc_i <- match(sozc,repository$electrode_list)

  for (i in seq_len(dim(f)[2])){
    mean_f_soz[i] <- mean(f[soz_i,i])
    se_f_soz[i] <- sd(f[soz_i,i])/sqrt(length(soz))
    mean_f_sozc[i] <- mean(f[sozc_i,i])
    se_f_sozc[i] <- sd(f[sozc_i,i]/sqrt(length(sozc)))
  }
  return(list(
    mean_f_soz = mean_f_soz,
    mean_f_sozc = mean_f_sozc,
    se_f_soz = se_f_soz,
    se_f_sozc = se_f_sozc
  ))
}

threshold_buckets <- function(mat, thresholds) {

  thresholds <- sort(thresholds)

  if (max(mat) > 1) {
    stop("Matrix values must be between 0 and 1!")
  }

  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      k <- 1
      while(mat[i,j] >= thresholds[k]) {
        k <- k + 1
      }
      mat[i,j] <- mean(thresholds[(k-1):k])
    }
  }

  mat <- (mat - min(mat))/(max(mat) - min(mat)) # normalize to between 0 and 1

  return(mat)
}

moving_average <- function(x, n) {
  x_ma <- stats::filter(x, rep(1 / n, n), sides = 2)
  x_ma[is.na(x_ma)] <- x[is.na(x_ma)] # replace NAs at beginning and end with original values
  return(as.vector(x_ma))
}

threshold_fragility <- function(repository, adj_frag_info, t_step, threshold_start, threshold_end, threshold = 0.5) {
  n_windows <- dim(adj_frag_info$adj)[3]

  # convert from input t_start and t_end to timewindow indices
  tw_start <- floor(which.min(abs(threshold_start-repository$voltage$dimnames$Time))/t_step)
  tw_end <- floor(which.min(abs(threshold_end-repository$voltage$dimnames$Time))/t_step)
  if (tw_end > n_windows) { tw_end <- n_windows }

  # subset fragility matrix to specified timewindows
  mat <- adj_frag_info$frag[,tw_start:tw_end]

  avg_f <- rowMeans(mat)
  elec <- which(avg_f > threshold)

  return(list(
    avg_f = avg_f,
    elecnames = attr(elec, "names")
  ))
}

voltage_resample <- function(v,fs_new,fs_old) {
  dumm<-t(electrodes_by_time)
  v_resample<-gsignal::resample(v,fs_new,fs_old)
  elec_resamp<-t(v_resample)
}
#
# ridgecv <- function(xt, xtp1, parallel=FALSE) {
#     if (!identical(dim(xt), dim(xtp1))) {
#         stop("Unmatched dimension")
#     }
#     nel <- ncol(xt)
#     x <- as.matrix(xt)
#     ## parallel computing backend
#     if (parallel){
#         library(doParallel)
#         library(foreach)
#         if(!getDoParRegistered()){
#             registerDoParallel(cores=parallel::detectCores())
#         }
#         A <- foreach(i = seq_len(nel), .packages=c("glmnet"), .export = "ridgecvTwoPass",.combine = cbind) %dopar%{
#             y <- xtp1[, i]
#             fitcoef <- ridgecvTwoPass(x,y)
#             fitcoef <- as.numeric(fitcoef)
#             fitcoef[-1]
#         }
#     }else{
#         ## Coefficient matrix A
#         ## each column is coefficients from a linear regression
#         ## formula: xtp1 = xt*A + E
#         A <- matrix(0, nel, nel)
#         ## for each electrode
#         for (i in seq_len(nel)) {
#             y <- xtp1[, i]
#             fitcoef <- tryCatch({
#                 ridgecvTwoPass(x,y)
#             }, error = function(e) {
#                 fit <- glmnet::glmnet(x, y,
#                     alpha = 0, lambda = 0,
#                     standardize = FALSE,
#                     intercept = FALSE
#                 )
#                 coef(fit)
#             })
#
#             fitcoef <- as.numeric(fitcoef)
#             A[, i] <- fitcoef[-1]
#         }
#         A
#     }
# }
#
# ridgecvTwoPass <- function(x, y, lambdaRange = 10^-rev(1:10)){
#     set.seed(1)
#     ## first pass: determine the scale of lambda
#     fitcv1 <- glmnet::cv.glmnet(x, y,
#         alpha = 0,
#         standardize = FALSE,
#         intercept = FALSE,
#         type.measure="mse",
#         lambda = lambdaRange,
#         thresh = 1e-10
#     )
#     # plot(fitcv1)
#     bestlambda <- fitcv1$lambda.min
#
#     ## if the best lambda is the smallest one, then we do not use ridge regression
#     if (bestlambda==min(lambdaRange)) {
#         fit <- glmnet::glmnet(x, y,
#             alpha = 0, lambda = 0,
#             standardize = FALSE,
#             intercept = FALSE
#         )
#         return(coef(fit))
#     }
#
#     ## second pass: determine the optimal lambda
#     lambdaRange2 <- seq(bestlambda/10, bestlambda*10, length.out = 1000)
#     fitcv2 <- glmnet::cv.glmnet(x, y,
#         alpha = 0,
#         standardize = FALSE,
#         intercept = FALSE,
#         type.measure="mse",
#         lambda = lambdaRange2,
#         thresh = 1e-10
#     )
#
#     coef(fitcv2, s = fitcv2$lambda.min)
# }

predictRidge <- function(xt, A) {
    ## the data matrix
    if (nrow(A) == ncol(A) + 1) {
        x <- cbind(1, as.matrix(xt))
    } else {
        x <- as.matrix(xt)
    }
    x %*% A
}
