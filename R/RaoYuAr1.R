#' @title Small Area Estimation using Hierarchical Bayesian under Rao-Yu Model
#' @description This function is implemented to variable of interest \code{ydi}
#' @param formula Formula that describe the fitted model
#' @param area Number of areas (domain) of the data
#' @param period Number of periods (subdomains) for each area of the data
#' @param vardir Sampling variances of direct estimations
#' @param iter.update Number of updates with default \code{3}
#' @param iter.mcmc Number of total iterations per chain with default \code{2000}
#' @param thin Thinning rate, must be a positive integer with default \code{1}
#' @param burn.in Number of iterations to discard at the beginning with default \code{1000}
#' @param tau.e Variance of area-by-time effect of variable interest with default \code{1}
#' @param tau.v Variance of random area effect of variable interest with default \code{1}
#' @param data The data frame
#' @export RaoYuAr1
#' @import stringr
#' @import coda
#' @import rjags
#' @import stats
#' @import grDevices
#' @import graphics
#' @return This function returns a list of the following objects:
#'    \item{Est}{A vector with the values of Small Area mean Estimates using Hierarchical bayesian method}
#'    \item{refVar}{Estimated random effect variances}
#'    \item{coefficient}{A dataframe with the estimated model coefficient}
#'    \item{alpha}{Parameter dispersion of Generalized Poisson distribution}
#'    \item{plot}{Trace, Density, Autocorrelation Function Plot of MCMC samples}
#'@examples
#' ##For data without any non-sampled area
#' data(dataAr1)     # Load dataset
#' formula = ydi ~ xdi1 + xdi2
#' area = max(dataAr1[, "area"])
#' period = max(dataAr1[,"period"])
#' vardir = dataAr1[,"vardir"]
#'
#' result <- RaoYuAr1(formula, area, period, vardir, data = dataAr1)
#' result$Est
#' result$refVar
#' result$coefficient
#' result$plot
#' ## For data with non-sampled area use dataAr1Ns
#' @export
RaoYuAr1<-function( formula, area, period, vardir, iter.update=3, iter.mcmc=2000,
                    thin = 2, burn.in =1000, tau.e = 1, tau.v=1, data){

  result <- list(Est = NA, refVar = NA, coefficient = NA, plot = NA)
  formuladata <- model.frame(formula, data, na.action = NULL)
  formuladata <- data.frame(vardir, formuladata)
  m=area
  t=period
  n_formuladata <- nrow(formuladata)

  if (any(is.na(formuladata[,-c(1:2)]))){ stop("Auxiliary Variables contains NA values.")
  } else if (iter.update < 3) {stop("the number of iteration updates at least 3 times")
  } else if (n_formuladata!=m*t) stop("length of variabel must be multiply of area and period")
  for (i in 1:n_formuladata) {
    if(!any(is.na(formuladata[i,2])) & any(is.na(formuladata[i, 1])))
      stop("any vardir contain NA when direct estimate are available")
  }

  vardirr = as.matrix(formuladata[,'vardir'])
  ydir = as.matrix(formuladata[,2])
  xdir = as.matrix(formuladata[,-c(1:2)])
  if(ncol(formuladata)>3){
    aux = ncol(formuladata[,-c(1:2)])
  }else{
    aux=1
  }
  nvar=aux+1
  y = matrix(0, nrow = m, ncol = t)
  vardir = matrix(0, nrow = m, ncol = t)
  k=0
  for (i in 1:m){
    for(j in 1:t){
      k = k+1
      vardir[i,j] = vardirr[k]
      y[i,j] = ydir[k,1]
    }
  }
  if (!any(is.na(formuladata[,2]))){
    x = list()
    for (i in 1:aux) {
      x[[i]] <- matrix(0, nrow = m, ncol = t)
    }
    kx=0
    for (i in 1:m){
      for(j in 1:t){
        kx = kx+1
        for (h in 1:aux){
          x[[h]][i,j] <- xdir[kx,h]
        }
      }
    }
    x_aux = c()
    for (r in 1:aux) {
      x_aux =  c(x_aux,x[[r]])
    }
    dim(x_aux) = c(m,t,aux)

    mu.b = rep(0, nvar)
    tau.b = rep(1, nvar)
    tau.va=tau.vb=1
    tau.ea=tau.eb=1
    a.var=1
    rho.a=-1
    rho.b=1

    for(iter in 1:iter.update){
      dat <- list("m"= m, "t"= t, "x" = x_aux, "vardir" = vardir, "y" = y, "nvar"=nvar, "aux"=aux,
                  "mu.b" = mu.b, "tau.b"=tau.b, "rho.a" = rho.a,"rho.b" = rho.b,"tau.va" = tau.va,
                  "tau.vb" = tau.vb, "tau.ea" = tau.ea,"tau.eb" = tau.eb)
      inits <- list(eps = matrix(0,m,t), b = mu.b, tau.v=1, tau.e = 1)
      cat("model {
        for (i in 1:m) {
				    v[i]~dnorm(0,tau.v)
				    for (j in 1:t){
				      y[i,j]~dnorm(mu[i,j],tau[i,j])
				      mu[i,j]<-b[1] + sum(b[2:nvar]*(x[i, j, 1:aux])) + v[i] + u[i,j]
				      eps[i,j]~dnorm(0,tau.e)
				    }
				    u[i,1]<-eps[i,1]
				    for(j in 2:t){
				      u[i,j]<-rho*u[i,j-1]+eps[i,j]
				    }
        }

				#Priors
				for (k in 1:nvar){
				    b[k] ~ dnorm(mu.b[k],tau.b[k])
				}
				rho ~ dunif(rho.a,rho.b)
        tau.v ~ dgamma(tau.va,tau.vb)
        tau.e ~ dgamma(tau.ea,tau.eb)
        a.var<-1/tau.v
        tau <-1/vardir
      }", file = "rao_yu.txt")

      jags.m <- jags.model( file = "rao_yu.txt", data=dat, inits=inits, n.chains=1, n.adapt=500 )
      file.remove("rao_yu.txt")
      params <- c("mu", "a.var", "rho", "b", "tau.v", "tau.e")
      samps <- coda.samples( jags.m, params, n.iter = iter.mcmc, thin = thin)
      samps1 <- window(samps, start = burn.in + 1, end = iter.mcmc)
      result_samps=summary(samps1)
      a.var=result_samps$statistics[1]
      beta=result_samps$statistics[2:(nvar+1),1:2]
      for (i in 1:nvar){
        mu.b[i]  = beta[i,1]
        tau.b[i] = 1/(beta[i,2]^2)
      }
      rho.a=result_samps$statistics[m*t+nvar+2,1]-sqrt(t)*result_samps$statistics[m*t+nvar+2,2]
      rho.b=result_samps$statistics[m*t+nvar+2,1]+sqrt(t)*result_samps$statistics[m*t+nvar+2,2]
      tau.ea = result_samps$statistics[m*t+nvar+3,1]^2/result_samps$statistics[m*t+nvar+3,1]^2
      tau.eb = result_samps$statistics[m*t+nvar+3,1]/result_samps$statistics[m*t+nvar+3,1]^2
      tau.va = result_samps$statistics[m*t+nvar+4,1]^2/result_samps$statistics[m*t+nvar+4,1]^2
      tau.vb = result_samps$statistics[m*t+nvar+4,1]/result_samps$statistics[m*t+nvar+4,1]^2
    }
    result_samps = summary(samps1)
    b.varnames <- list()
    for (i in 1:(nvar)) {
      idx.b.varnames <- as.character(i-1)
      b.varnames[i] <-str_replace_all(paste("b[",idx.b.varnames,"]"),pattern=" ", replacement="")
    }

    result_mcmc <- samps1[,c(2:(nvar+1))]
    colnames(result_mcmc[[1]]) <- b.varnames
    a.var=result_samps$statistics[1]
    beta=result_samps$statistics[2:(nvar+1),1:2]
    rownames(beta) <- b.varnames
    rho <- result_samps$statistics[1+nvar+(m*t),1:2]
    mu=result_samps$statistics[(nvar+2):(1+nvar+(m*t)),1:2]
    Estimation=data.frame(mu)
    coef = rbind(beta, rho)
    rownames(coef) <- c(b.varnames, "rho")

    Quantiles <- as.data.frame(result_samps$quantiles)
    q_mu <- Quantiles[(nvar+2):(1+nvar+(m*t)),]
    q_beta <- (Quantiles[2:(nvar+1),])
    q_rho <- Quantiles[2+nvar+(m*t),]
    q_coef <- rbind(q_beta, q_rho)
    rownames(q_coef) <- c(b.varnames, "rho")
    coef <- cbind(coef, q_coef)
    Estimation <- data.frame(Estimation,q_mu)
    mean_est = matrix(0, t, m)
    sd_est = matrix(0, t, m)
    est1 = matrix(0, t, m)
    est2 = matrix(0, t, m)
    est3 = matrix(0, t, m)
    est4 = matrix(0, t, m)
    est5 = matrix(0, t, m)

    k=0
    for (i in 1:t) {
      for (j in 1:m) {
        k=k+1
        mean_est[i,j] = Estimation[k,1]
        sd_est[i,j] = Estimation[k,2]
        est1[i,j] = Estimation[k,3]
        est2[i,j] = Estimation[k,4]
        est3[i,j] = Estimation[k,5]
        est4[i,j] = Estimation[k,6]
        est5[i,j] = Estimation[k,7]
      }
    }
    mean_est = t(mean_est)
    sd_est = t(sd_est)
    est1 = t(est1)
    est2 = t(est2)
    est3 = t(est3)
    est4 = t(est4)
    est5 = t(est5)
    Mean_est = c()
    Sd_est = c()
    Est1 = c()
    Est2 = c()
    Est3 = c()
    Est4 = c()
    Est5 = c()

    k = 0
    for (i in 1:m) {
      for (j in 1:t) {
        k=k+1
        Mean_est[k] = mean_est[i,j]
        Sd_est[k] = sd_est[i,j]
        Est1[k] = est1[i,j]
        Est2[k] = est2[i,j]
        Est3[k] = est3[i,j]
        Est4[k] = est4[i,j]
        Est5[k] = est5[i,j]
      }
    }
    Estimation <- data.frame(Mean_est, Sd_est, Est1, Est2, Est3, Est4, Est5)
    colnames(Estimation) <- c("MEAN","SD","2.5%","25%","50%","75%","97.5%")
  }else
  {
    Y = as.matrix(na.omit(y))
    V = na.omit(vardir)
    rowNA <- c()
    a=0
    for (i in 1:m) {
      if(is.na(y[i,1])){
        a = a+1
        rowNA[a] <- i
      }
    }

    x = list()
    for (i in 1:aux) {
      x[[i]] <- matrix(0, nrow = m, ncol = t)
    }
    k=0
    for (j in 1:m) {
      for(l in 1:t){
        k = k+1
        for (i in 1:aux) {
          x[[i]][j,l] <- xdir[k,i]
        }
      }
    }
    x_aux = c()
    for (r in 1:aux) {
      x_aux =  c(x_aux,x[[r]])
    }
    dim(x_aux) = c(m,t,aux)
    NTS = length(rowNA)
    NS = m-NTS
    x_auxS <- x_aux[-rowNA,,]
    dim(x_auxS) = c(NS,t,aux)
    x_auxTS <- x_aux[rowNA,,]
    dim(x_auxTS) = c(NTS,t,aux)

    mu.b = rep(0, nvar)
    tau.b = rep(1, nvar)
    tau.va=tau.vb=1
    tau.ea=tau.eb=1
    a.var=1
    rho.a=-1
    rho.b=1
    for (iter in 1:iter.update) {
      dat <- list("NS"=NS,"NTS"=NTS,"t"=t, "nvar"=nvar, "aux"=aux, "y"=Y,
                  "xS"=x_auxS, "xTS"=x_auxTS, "vardir"=V, "mu.b"=mu.b,"tau.b"=tau.b,
                  "tau.va"=tau.va, "tau.vb"=tau.vb, "tau.ea"=tau.ea, "tau.eb"=tau.eb,
                  "rho.a"=rho.a,"rho.b"=rho.b)  # names list of numbers
      inits <- list(eps = matrix(0,NS,t), b = mu.b)
      cat("model {
					for (i in 1:NS) {
					v[i]~dnorm(0,tau.v)
					for (j in 1:t){
		  			y[i,j]~dnorm(mu[i,j],tau[i,j])
				  	mu[i,j]<-b[1] + sum(b[2:nvar]*(xS[i, j, 1:aux])) + v[i] + u[i,j]
					  eps[i,j]~dnorm(0,tau.e)
					}
					  u[i,1]<-eps[i,1]
					for(j in 2:t){
					  u[i,j]<-rho*u[i,j-1]+eps[i,j]
					  }
					}

					for (i in 1:NTS) {
					  vT[i]~dnorm(0,tau.v)
					  for (j in 1:t){
					    muT[i,j]<-mu.b[1] + sum(mu.b[2:nvar]*(xTS[i, j, 1:aux])) + v[i] + u[i,j]
					    epsT[i,j]~dnorm(0,tau.e)
					  }
					  uT[i,1]<-epsT[i,1]
					  for(j in 2:t){
					    uT[i,j]<-rhoT*uT[i,j-1]+epsT[i,j]
					  }
					}

					#priors
					for (k in 1:nvar){
					    b[k] ~ dnorm(mu.b[k],tau.b[k])
					}
					  rhoT=rho
					  rho~dunif(rho.a,rho.b)
					  tau.e~dgamma(tau.ea,tau.eb)
					  tau.v~dgamma(tau.va,tau.vb)
					  a.var <- 1/tau.v
					  tau <-1/vardir
			  }", file="rao_yu.txt")
      jags.m <- jags.model( file = "rao_yu.txt", data=dat, inits=inits, n.chains=1, n.adapt=500 )
      file.remove("rao_yu.txt")
      params <- c("mu","muT", "a.var", "rho", "b", "tau.v", "tau.e")
      samps <- coda.samples( jags.m, params, n.iter = iter.mcmc, thin = thin)
      samps1 <- window(samps, start = burn.in, end = iter.mcmc)
      result_samps=summary(samps1)
      a.var=result_samps$statistics[1]
      beta=result_samps$statistics[2:(nvar+1),1:2]
      for (i in 1:nvar){
        mu.b[i]  = beta[i,1]
        tau.b[i] = 1/(beta[i,2]^2)
      }

      rho.a=result_samps$statistics[m*t+nvar+2,1]-sqrt(t)*result_samps$statistics[m*t+nvar+2,2]
      rho.b=result_samps$statistics[m*t+nvar+2,1]+sqrt(t)*result_samps$statistics[m*t+nvar+2,2]
      tau.ea = result_samps$statistics[m*t+nvar+3,1]^2/result_samps$statistics[m*t+nvar+3,1]^2
      tau.eb = result_samps$statistics[m*t+nvar+3,1]/result_samps$statistics[m*t+nvar+3,1]^2
      tau.va = result_samps$statistics[m*t+nvar+4,1]^2/result_samps$statistics[m*t+nvar+4,1]^2
      tau.vb = result_samps$statistics[m*t+nvar+4,1]/result_samps$statistics[m*t+nvar+4,1]^2
    }
    result_samps = summary(samps1)
    b.varnames <- list()
    for (i in 1:(nvar)) {
      idx.b.varnames <- as.character(i-1)
      b.varnames[i] <-str_replace_all(paste("b[",idx.b.varnames,"]"),pattern=" ", replacement="")
    }
    result_mcmc <- samps1[,c(2:(nvar+1))]
    colnames(result_mcmc[[1]]) <- b.varnames
    a.var=result_samps$statistics[1]
    beta=result_samps$statistics[2:(nvar+1),1:2]
    rownames(beta) <- b.varnames
    mu_stat=result_samps$statistics[(nvar+2):(NS*t+nvar+1),1:2]
    mu_quant = result_samps$quantiles[(nvar+2):(NS*t+nvar+1),]
    mu = cbind(mu_stat, mu_quant )
    muT_stat=result_samps$statistics[(NS*t+nvar+2):(m*t+nvar+1),1:2]
    muT_quant=result_samps$quant[(NS*t+nvar+2):(m*t+nvar+1),]
    muT = cbind(muT_stat, muT_quant)

    Mu_mean_est = matrix(0, t, m)
    Mu_sd_est = matrix(0, t, m)
    Mu_est1 = matrix(0, t, m)
    Mu_est2 = matrix(0, t, m)
    Mu_est3 = matrix(0, t, m)
    Mu_est4 = matrix(0, t, m)
    Mu_est5 = matrix(0, t, m)

    k=0
    for (i in 1:t) {
      for (j in 1:(m-length(rowNA))) {
        k=k+1
        Mu_mean_est[i,j] = mu[k,1]
        Mu_sd_est[i,j] = mu[k,2]
        Mu_est1[i,j] = mu[k,3]
        Mu_est2[i,j] = mu[k,4]
        Mu_est3[i,j] = mu[k,5]
        Mu_est4[i,j] = mu[k,6]
        Mu_est5[i,j] = mu[k,7]
      }
    }
    Mu_mean_est = t(Mu_mean_est)
    Mu_sd_est = t(Mu_sd_est)
    Mu_est1 = t(Mu_est1)
    Mu_est2 = t(Mu_est2)
    Mu_est3 = t(Mu_est3)
    Mu_est4 = t(Mu_est4)
    Mu_est5 = t(Mu_est5)
    Mu_Mean_est = c()
    Mu_Sd_est = c()
    Mu_Est1 = c()
    Mu_Est2 = c()
    Mu_Est3 = c()
    Mu_Est4 = c()
    Mu_Est5 = c()

    k = 0
    for (i in 1:(m-length(rowNA))) {
      for (j in 1:t) {
        k=k+1
        Mu_Mean_est[k] = Mu_mean_est[i,j]
        Mu_Sd_est[k] = Mu_sd_est[i,j]
        Mu_Est1[k] = Mu_est1[i,j]
        Mu_Est2[k] = Mu_est2[i,j]
        Mu_Est3[k] = Mu_est3[i,j]
        Mu_Est4[k] = Mu_est4[i,j]
        Mu_Est5[k] = Mu_est5[i,j]
      }
    }
    mu <- data.frame(Mu_Mean_est, Mu_Sd_est, Mu_Est1, Mu_Est2, Mu_Est3, Mu_Est4, Mu_Est5)
    colnames(mu) <- c("MEAN", "SD", "2.5%", "25%", "50%", "75%", "97.5%")

    MuT_mean_est = matrix(0, t, m)
    MuT_sd_est = matrix(0, t, m)
    MuT_est1 = matrix(0, t, m)
    MuT_est2 = matrix(0, t, m)
    MuT_est3 = matrix(0, t, m)
    MuT_est4 = matrix(0, t, m)
    MuT_est5 = matrix(0, t, m)

    k=0
    for (i in 1:t) {
      for (j in 1:(length(rowNA))) {
        k=k+1
        MuT_mean_est[i,j] = muT[k,1]
        MuT_sd_est[i,j] = muT[k,2]
        MuT_est1[i,j] = muT[k,3]
        MuT_est2[i,j] = muT[k,4]
        MuT_est3[i,j] = muT[k,5]
        MuT_est4[i,j] = muT[k,6]
        MuT_est5[i,j] = muT[k,7]
      }
    }
    MuT_mean_est = t(MuT_mean_est)
    MuT_sd_est = t(MuT_sd_est)
    MuT_est1 = t(MuT_est1)
    MuT_est2 = t(MuT_est2)
    MuT_est3 = t(MuT_est3)
    MuT_est4 = t(MuT_est4)
    MuT_est5 = t(MuT_est5)
    MuT_Mean_est = c()
    MuT_Sd_est = c()
    MuT_Est1 = c()
    MuT_Est2 = c()
    MuT_Est3 = c()
    MuT_Est4 = c()
    MuT_Est5 = c()

    k = 0
    for (i in 1:(length(rowNA))) {
      for (j in 1:t) {
        k=k+1
        MuT_Mean_est[k] = MuT_mean_est[i,j]
        MuT_Sd_est[k] = MuT_sd_est[i,j]
        MuT_Est1[k] = MuT_est1[i,j]
        MuT_Est2[k] = MuT_est2[i,j]
        MuT_Est3[k] = MuT_est3[i,j]
        MuT_Est4[k] = MuT_est4[i,j]
        MuT_Est5[k] = MuT_est5[i,j]
      }
    }
    muT <- data.frame(MuT_Mean_est, MuT_Sd_est, MuT_Est1, MuT_Est2, MuT_Est3, MuT_Est4, MuT_Est5)
    colnames(muT) <- c("MEAN", "SD", "2.5%", "25%", "50%", "75%", "97.5%")

    rho=result_samps$statistics[(m*t+nvar+2),1:2]
    coef <- rbind(beta,rho)
    result_all = rbind(mu,muT)
    idx.mu = 1
    idx.muT = 1
    idx = 0
    for(i in 1:m){
      for(j in 1:t){
        idx=idx+1
        if(data[idx,2] %in% rowNA){
          result_all[idx,]<-muT[idx.muT,]
          idx.muT = idx.muT+1
        }else{
          result_all[idx,]<-mu[idx.mu,]
          idx.mu = idx.mu+1
        }
      }
    }
    Estimation = data.frame(result_all)

    q_beta <- result_samps$quantiles[2:(nvar+1),]
    q_rho <- result_samps$quantiles[(m*t+nvar+2),]
    q_coef <- rbind(q_beta,q_rho)
    rownames(q_coef) <- c(b.varnames,"rho")
    coef <- cbind(coef,q_coef)

  }
  result$Est = Estimation
  result$refVar = a.var
  result$coefficient = coef
  result$plot = list(graphics.off(), par(mar = c(2, 2, 2, 2)),
                     autocorr.plot(result_mcmc, col = "brown2", lwd = 2),
                     plot(result_mcmc, col = "brown2", lwd = 2))
  return(result)
}
