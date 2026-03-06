#----- Functions for multivariate Gaussian mixture serology model -----#
library(emdbook)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(matrixStats)
library(stringr)
library(lhs)
library(mvtnorm)
library(matrixcalc)


#----- Generate infection status matrix
inf_matrix <- function(N_pathogen, pres = rep(1, N_pathogen)) {
  # list of possible outcomes for each pathogen
  combos <- list()
  for (c in 1:N_pathogen) {
    combos[[c]] <- c(0, 1)
  }

  # matrix of all possible infection status combinations
  m <- expand.grid(combos)
  colnames(m) <- letters[1:N_pathogen]

  # remove positives of absent pathogens
  if (sum(pres) < N_pathogen) {
    for (abs in which(pres == 0)) {
      m <- m[m[, abs] == 0, ]
    }
  }

  return(m)
}


#----- Extract prevalence estimates
extract_sero <- function(chains, data, pathogens) {
  sero <- data.frame(pathogen = pathogens, med = NA, ciL = NA, ciU = NA)

  for (p in 1:data$nP) {
    sero[p, 2:4] <- quantile(
      chains[, paste('seroAll[', paste(p, ']', sep = ''), sep = '')],
      c(0.5, 0.025, 0.975)
    )
  }
  sero <- sero[!sero$med == 0, ]

  return(sero)
}


#----- Extract prevalence estimates by location
extract_seroLoc <- function(chains, data, pathogens, loc) {
  sero <- data.frame(pathogen = NA, loc = NA, med = NA, ciL = NA, ciU = NA)
  ind <- 1
  for (p in 1:data$nPp) {
    for (l in 1:data$nL) {
      lp <- paste(l, p, sep = ',')
      sero[ind, 3:5] <- quantile(
        chains[, paste('seroLoc[', paste(lp, ']', sep = ''), sep = '')],
        c(0.5, 0.025, 0.975)
      )
      sero[ind, 1:2] <- c(pathogens[p], loc[l])
      ind <- ind + 1
    }
  }

  return(sero)
}


#----- Extract prevalence estimates by age
extract_seroAge <- function(chains, data, pathogens, ageG) {
  sero <- data.frame(pathogen = NA, age = NA, med = NA, ciL = NA, ciU = NA)
  ind <- 1
  for (p in 1:data$nPp) {
    for (a in 1:data$nA) {
      ap <- paste(a, p, sep = ',')
      sero[ind, 3:5] <- quantile(
        chains[, paste('seroAge[', paste(ap, ']', sep = ''), sep = '')],
        c(0.5, 0.025, 0.975)
      )
      sero[ind, 1:2] <- c(pathogens[p], ageG[a])
      ind <- ind + 1
    }
  }

  return(sero)
}


#----- Extract prevalence estimates by location & age
extract_seroLocAge <- function(chains, data, pathogens, loc, ageG) {
  sero <- data.frame(
    loc = NA,
    age = NA,
    pathogen = NA,
    med = NA,
    ciL = NA,
    ciU = NA
  )
  ind <- 1
  for (p in 1:data$nPp) {
    for (l in 1:data$nL) {
      for (a in 1:data$nA) {
        k <- paste(l, a, sep = ',')
        sero[ind, 1:3] <- c(loc[l], ageG[a], pathogens[p])
        sero[ind, 4:6] <- quantile(
          chains[, paste(
            'sero[',
            paste(k, paste(p, ']', sep = ''), sep = ','),
            sep = ''
          )],
          c(0.5, 0.025, 0.975)
        )
        ind <- ind + 1
      }
    }
  }

  return(sero)
}


#----- Extract cross-reactivity & correlation estimates
extract_phi <- function(chains, data, pathogens) {
  phi <- data.frame(pos = NA, neg = NA, med = NA, ciL = NA, ciU = NA)
  ind <- 1
  for (p in 1:data$nPp) {
    for (p2 in 1:data$nP) {
      if (!p == p2) {
        phi[ind, 1:2] <- c(pathogens[p], pathogens[p2])
        y <- str_replace_all(toString(c(p, p2)), " ", "")
        phi[ind, 3:5] <- quantile(
          chains[, paste('CR[', paste(y, ']', sep = ''), sep = '')],
          c(0.5, 0.025, 0.975)
        )
        ind <- ind + 1
      }
    }
  }

  rho <- data.frame(pars = c('rho00'), med = NA, ciL = NA, ciU = NA)
  rho[1, 2:4] <- quantile(chains[, paste('rho00')], c(0.5, 0.025, 0.975))

  return(list(phi = phi, rho = rho))
}


#----- Extract gaussian means
extract_mu <- function(chains, data, pathogens) {
  # label combination positives
  pos <- rep('neg', data$nC)
  for (c in 1:data$nC) {
    np <- sum(data$infM[c, ])
    if (np == 1) {
      pos[c] <- pathogens[which(data$infM[c, ] == 1)]
    } else if (np > 1) {
      pos[c] <- paste(
        pathogens[which(data$infM[c, ] == 1)],
        sep = '&',
        collapse = "&"
      )
    }
  }

  # all gaussian means
  mus0 <- data.frame(
    pg = rep(NA, length(which(data$infM == 0))),
    pos = NA,
    med = NA,
    ciL = NA,
    ciU = NA
  )
  mus1 <- data.frame(
    pg = rep(NA, length(which(data$infM == 1))),
    pos = NA,
    med = NA,
    ciL = NA,
    ciU = NA
  )
  ix0 <- 1
  ix1 <- 1
  for (c in 1:data$nC) {
    for (p in 1:data$nP) {
      if (data$infM[c, p] == 0) {
        mus0$pg[ix0] <- pathogens[p]
        mus0$pos[ix0] <- pos[c]
        y <- paste(c, p, sep = ',')
        mus0[ix0, 3:5] <- quantile(
          chains[, paste(paste('mu[', y, sep = ''), ']', sep = '')],
          c(0.5, 0.025, 0.975)
        )
        ix0 <- ix0 + 1
      } else {
        mus1$pg[ix1] <- pathogens[p]
        mus1$pos[ix1] <- pos[c]
        y <- paste(c, p, sep = ',')
        mus1[ix1, 3:5] <- quantile(
          chains[, paste(paste('mu[', y, sep = ''), ']', sep = '')],
          c(0.5, 0.025, 0.975)
        )
        ix1 <- ix1 + 1
      }
    }
  }
  colnames(mus0)[1] <- 'antigen'
  colnames(mus1)[1] <- 'antigen'

  return(list(mus0 = mus0, mus1 = mus1))
}


#----- Extract gaussian sds
extract_sds <- function(chains, data) {
  sig <- data.frame(par = c('sd0', 'sd1'), med = NA, ciL = NA, ciU = NA)
  sig[1, 2:4] <- quantile(chains$sd0, c(0.5, 0.025, 0.975))
  sig[2, 2:4] <- quantile(chains$sd1, c(0.5, 0.025, 0.975))

  return(sig)
}


#----- Extract gaussian sds for ELISA
extract_sdsELISA <- function(chains, data, pathogens) {
  sig <- data.frame(
    par = paste('sdE', seq(1, data$nC), sep = '_'),
    med = NA,
    ciL = NA,
    ciU = NA
  )
  sdE <- paste('sde[', paste(seq(1, data$nC), ']', sep = ''), sep = '')
  for (i in 1:data$nC) {
    sig[i, 2:4] <- quantile(chains[, sdE[i]], c(0.5, 0.025, 0.975))
  }

  # label
  sig$status <- 'negative'
  for (c in 2:nrow(sig)) {
    sig$status[c] <- paste(
      pathogens[data$wpos[c, data$wpos[c, ] > 0]],
      collapse = '+'
    )
  }

  return(sig)
}


#----- Extract covariance matrices per iteration
extract_covM <- function(chains, data) {
  iter <- length(chains$lp__)
  covM <- list()
  for (i in 1:iter) {
    covM[[i]] <- list()
    for (c in 1:data$nC) {
      x <- matrix(NA, ncol = data$nP, nrow = data$nP)
      for (p in 1:data$nP) {
        for (p2 in 1:data$nP) {
          y <- paste(paste(c, p, sep = ','), p2, sep = ',')
          x[p, p2] <- chains[
            i,
            paste(paste('covM[', y, sep = ''), ']', sep = '')
          ]
        }
      }
      covM[[i]][[c]] <- x
    }
  }
  return(covM)
}


#----- Plot gaussian distribution fits
plot_dists <- function(chains, data, pathogens) {
  iter <- length(chains$lp__)
  covM <- extract_covM(chains, data)

  # simulate multivariate gaussians per combination
  yy <- yp <- list()
  for (p in 1:data$nP) {
    yp[[p]] <- matrix(NA, nrow = 512, ncol = iter)
  }
  ypN <- ypP <- yp
  for (i in 1:iter) {
    for (c in 1:data$nC) {
      # simulate gaussion for combination c, iteration i
      g <- paste('mu[', c, sep = '')
      nn <- ceiling(sum(
        data$N * chains[, paste('theta[', paste(c, ']', sep = ''), sep = '')][i]
      ))
      muu <- vector()
      for (p in 1:data$nP) {
        muu[p] <- chains[i, paste(paste(g, p, sep = ','), ']', sep = '')]
      }
      yy[[c]] <- as.data.frame(rmvnorm(nn, mean = muu, sigma = covM[[i]][[c]]))
      yy[[c]]$C <- c
    }

    # density distributions per pathogen
    yc <- do.call('rbind', yy)
    for (p in 1:data$nP) {
      yp[[p]][, i] <- density(yc[, p], bw = 0.01, from = -2.5, to = 5.5)$y

      if (data$pres[p] == 1) {
        z <- which(data$infM[, p] == 1)
        pw <- vector()
        for (s in 1:length(z)) {
          pw[s] <- chains[
            i,
            paste(paste('theta[', paste(z[s]), sep = ''), ']', sep = '')
          ]
        }
        propP <- sum(data$N * pw) / data$N
        ypP[[p]][, i] <- density(
          yc[yc$C %in% which(data$infM[, p] == 1), p],
          bw = 0.01,
          from = -2.5,
          to = 5.5
        )$y *
          propP
        ypN[[p]][, i] <- density(
          yc[yc$C %in% which(data$infM[, p] == 0), p],
          bw = 0.01,
          from = -2.5,
          to = 5.5
        )$y *
          (1 - propP)
      } else {
        ypN[[p]][, i] <- density(
          yc[yc$C %in% which(data$infM[, p] == 0), p],
          bw = 0.01,
          from = -2.5,
          to = 5.5
        )$y
      }
    }
  }

  # quantiles of density distributions
  dpq <- dpqP <- dpqN <- list()
  titer <- density(yc[, 1], bw = 0.01, from = -2.5, to = 5.5)$x
  for (p in 1:data$nP) {
    dpq[[p]] <- as.data.frame(rowQuantiles(
      yp[[p]],
      probs = c(0.5, 0.025, 0.975)
    ))
    dpqN[[p]] <- as.data.frame(rowQuantiles(
      ypN[[p]],
      probs = c(0.5, 0.025, 0.975)
    ))
    dpqP[[p]] <- as.data.frame(rowQuantiles(
      ypP[[p]],
      probs = c(0.5, 0.025, 0.975)
    ))
    dpq[[p]]$titer <- dpqN[[p]]$titer <- dpqP[[p]]$titer <- titer
    dpq[[p]]$pathogen <- dpqN[[p]]$pathogen <- dpqP[[p]]$pathogen <- pathogens[
      p
    ]
  }
  dpq <- do.call('rbind', dpq)
  dpqN <- do.call('rbind', dpqN)
  dpqP <- do.call('rbind', dpqP)
  colnames(dpq)[1:3] <- colnames(dpqN)[1:3] <- colnames(dpqP)[1:3] <- c(
    'med',
    'ciL',
    'ciU'
  )

  # compile data for plotting
  dta <- as.data.frame(data$y)
  colnames(dta) <- pathogens
  dta <- tidyr::gather(dta, key = 'pathogen', value = 't')

  # overall fit
  fitD <- ggplot() +
    geom_histogram(
      data = dta,
      aes(t, y = ..density..),
      bins = 150,
      fill = 'grey80',
      col = 'grey70'
    ) +
    theme_minimal() +
    theme(text = element_text(size = 18)) + #xlim(-1.5,3.5)+
    geom_line(data = dpq, aes(titer, med), col = 'springgreen4') +
    facet_wrap(~pathogen, scales = 'free_y') +
    xlab('titer') +
    geom_ribbon(
      data = dpq,
      aes(x = titer, y = med, ymin = ciL, ymax = ciU),
      fill = 'springgreen3',
      alpha = 0.4
    )

  # pos-neg fit
  fitDPN <- ggplot() +
    geom_histogram(
      data = dta,
      aes(t, y = ..density..),
      bins = 150,
      fill = 'grey80',
      col = 'grey70'
    ) +
    theme_minimal() +
    theme(text = element_text(size = 18)) + #xlim(-1.5,3.5)+
    geom_line(data = dpqN, aes(titer, med), col = 'mediumblue') +
    geom_line(data = dpqP, aes(titer, med), col = 'violetred') +
    facet_wrap(~pathogen, scales = 'free_y') +
    xlab('titer') +
    geom_ribbon(
      data = dpqN,
      aes(x = titer, y = med, ymin = ciL, ymax = ciU),
      fill = 'mediumblue',
      alpha = 0.3
    ) +
    geom_ribbon(
      data = dpqP,
      aes(x = titer, y = med, ymin = ciL, ymax = ciU),
      fill = 'violetred',
      alpha = 0.5
    )

  # return plots
  return(list(fit = fitD, fitPN = fitDPN))
}


#----- Plot gaussian distribution fits for location model
plot_distsLoc <- function(chains, data, pathogens) {
  iter <- length(chains$lp__)
  covM <- extract_covM(chains, data)

  # simulate multivariate gaussians per combination
  yy <- yp <- list()
  for (p in 1:data$nP) {
    yp[[p]] <- matrix(NA, nrow = 512, ncol = iter)
  }
  ypN <- ypP <- yp
  for (i in 1:iter) {
    for (c in 1:data$nC) {
      # simulate gaussion for combination c, iteration i
      g <- paste('mu[', c, sep = '')
      nn <- 0
      for (l in 1:data$nL) {
        lc <- paste(l, c, sep = ',')
        nn <- nn +
          round(
            data$NperL[l] *
              chains[, paste('theta[', paste(lc, ']', sep = ''), sep = '')][i]
          )
      }

      muu <- vector()
      for (p in 1:data$nP) {
        muu[p] <- chains[i, paste(paste(g, p, sep = ','), ']', sep = '')]
      }
      if (nn == 0) {
        nn <- 1
      }
      yy[[c]] <- as.data.frame(rmvnorm(nn, mean = muu, sigma = covM[[i]][[c]]))
      yy[[c]]$C <- c
    }

    # density distributions per pathogen
    yc <- do.call('rbind', yy)
    for (p in 1:data$nP) {
      yp[[p]][, i] <- density(yc[, p], bw = 0.01, from = -2.5, to = 5.5)$y

      if (data$pres[p] == 1) {
        z <- which(data$infM[, p] == 1)
        pw <- rep(0, data$nL)
        for (l in 1:data$nL) {
          for (s in 1:length(z)) {
            lz <- paste(l, z[s], sep = ',')
            pw[l] <- pw[l] +
              chains[i, paste(paste('theta[', lz, sep = ''), ']', sep = '')]
          }
        }

        propP <- sum(data$NperL * pw) / data$N
        ypP[[p]][, i] <- density(
          yc[yc$C %in% which(data$infM[, p] == 1), p],
          bw = 0.01,
          from = -2.5,
          to = 5.5
        )$y *
          propP
        ypN[[p]][, i] <- density(
          yc[yc$C %in% which(data$infM[, p] == 0), p],
          bw = 0.01,
          from = -2.5,
          to = 5.5
        )$y *
          (1 - propP)
      } else {
        ypN[[p]][, i] <- density(
          yc[yc$C %in% which(data$infM[, p] == 0), p],
          bw = 0.01,
          from = -2.5,
          to = 5.5
        )$y
      }
    }
  }

  # quantiles of density distributions
  dpq <- dpqP <- dpqN <- list()
  titer <- density(yc[, 1], bw = 0.01, from = -2.5, to = 5.5)$x
  for (p in 1:data$nP) {
    dpq[[p]] <- as.data.frame(rowQuantiles(
      yp[[p]],
      probs = c(0.5, 0.025, 0.975)
    ))
    dpqN[[p]] <- as.data.frame(rowQuantiles(
      ypN[[p]],
      probs = c(0.5, 0.025, 0.975)
    ))
    dpqP[[p]] <- as.data.frame(rowQuantiles(
      ypP[[p]],
      probs = c(0.5, 0.025, 0.975)
    ))
    dpq[[p]]$titer <- dpqN[[p]]$titer <- dpqP[[p]]$titer <- titer
    dpq[[p]]$pathogen <- dpqN[[p]]$pathogen <- dpqP[[p]]$pathogen <- pathogens[
      p
    ]
  }
  dpq <- do.call('rbind', dpq)
  dpqN <- do.call('rbind', dpqN)
  dpqP <- do.call('rbind', dpqP)
  colnames(dpq)[1:3] <- colnames(dpqN)[1:3] <- colnames(dpqP)[1:3] <- c(
    'med',
    'ciL',
    'ciU'
  )

  # compile data for plotting
  dta <- as.data.frame(data$y)
  colnames(dta) <- pathogens
  dta <- tidyr::gather(dta, key = 'pathogen', value = 't')

  # overall fit
  fitD <- ggplot() +
    geom_histogram(
      data = dta,
      aes(t, y = ..density..),
      bins = 150,
      fill = 'grey80',
      col = 'grey70'
    ) +
    theme_minimal() +
    theme(text = element_text(size = 18)) + #xlim(-1.5,3.5)+
    geom_line(data = dpq, aes(titer, med), col = 'springgreen4') +
    facet_wrap(~pathogen, scales = 'free_y') +
    xlab('titer') +
    geom_ribbon(
      data = dpq,
      aes(x = titer, y = med, ymin = ciL, ymax = ciU),
      fill = 'springgreen3',
      alpha = 0.4
    )

  # pos-neg fit
  fitDPN <- ggplot() +
    geom_histogram(
      data = dta,
      aes(t, y = ..density..),
      bins = 150,
      fill = 'grey80',
      col = 'grey70'
    ) +
    theme_minimal() +
    theme(text = element_text(size = 18)) + #xlim(-1.5,3.5)+
    geom_line(data = dpqN, aes(titer, med), col = 'mediumblue') +
    geom_line(data = dpqP, aes(titer, med), col = 'violetred') +
    facet_wrap(~pathogen, scales = 'free_y') +
    xlab('titer') +
    geom_ribbon(
      data = dpqN,
      aes(x = titer, y = med, ymin = ciL, ymax = ciU),
      fill = 'mediumblue',
      alpha = 0.3
    ) +
    geom_ribbon(
      data = dpqP,
      aes(x = titer, y = med, ymin = ciL, ymax = ciU),
      fill = 'violetred',
      alpha = 0.5
    )

  # return plots
  return(list(fit = fitD, fitPN = fitDPN))
}


#----- Plot gaussian distribution fits for location & age model
plot_distsLocAge <- function(chains, data, pathogens, NperLA) {
  iter <- length(chains$lp__)
  covM <- extract_covM(chains, data)

  # simulate multivariate gaussians per combination
  yy <- yp <- list()
  for (p in 1:data$nP) {
    yp[[p]] <- matrix(NA, nrow = 512, ncol = iter)
  }
  ypN <- ypP <- yp
  for (i in 1:iter) {
    for (c in 1:data$nC) {
      # simulate gaussion for combination c, iteration i
      g <- paste('mu[', c, sep = '')
      nn <- 0
      for (l in 1:data$nL) {
        for (a in 1:data$nA) {
          la <- paste(l, a, sep = ',')
          nn <- nn +
            round(
              NperLA[a, l] *
                chains[, paste(
                  'theta[',
                  paste(paste(la, c, sep = ','), ']', sep = ''),
                  sep = ''
                )][i]
            )
        }
      }

      muu <- vector()
      for (p in 1:data$nP) {
        muu[p] <- chains[i, paste(paste(g, p, sep = ','), ']', sep = '')]
      }
      if (nn == 0) {
        nn <- 1
      }
      yy[[c]] <- as.data.frame(rmvnorm(nn, mean = muu, sigma = covM[[i]][[c]]))
      yy[[c]]$C <- c
    }

    # density distributions per pathogen
    yc <- do.call('rbind', yy)
    for (p in 1:data$nP) {
      yp[[p]][, i] <- density(yc[, p], bw = 0.01, from = -2.5, to = 5.5)$y

      if (data$pres[p] == 1) {
        z <- which(data$infM[, p] == 1)
        pw <- matrix(0, nrow = data$nA, ncol = data$nL)
        for (l in 1:data$nL) {
          for (a in 1:data$nA) {
            for (s in 1:length(z)) {
              la <- paste(l, a, sep = ',')
              pw[a, l] <- pw[a, l] +
                chains[
                  i,
                  paste(
                    paste('theta[', paste(la, z[s], sep = ','), sep = ''),
                    ']',
                    sep = ''
                  )
                ]
            }
          }
        }

        propP <- sum(NperLA * pw) / data$N
        ypP[[p]][, i] <- density(
          yc[yc$C %in% which(data$infM[, p] == 1), p],
          bw = 0.01,
          from = -2.5,
          to = 5.5
        )$y *
          propP
        ypN[[p]][, i] <- density(
          yc[yc$C %in% which(data$infM[, p] == 0), p],
          bw = 0.01,
          from = -2.5,
          to = 5.5
        )$y *
          (1 - propP)
      } else {
        ypN[[p]][, i] <- density(
          yc[yc$C %in% which(data$infM[, p] == 0), p],
          bw = 0.01,
          from = -2.5,
          to = 5.5
        )$y
      }
    }
  }

  # quantiles of density distributions
  dpq <- dpqP <- dpqN <- list()
  titer <- density(yc[, 1], bw = 0.01, from = -2.5, to = 5.5)$x
  for (p in 1:data$nP) {
    dpq[[p]] <- as.data.frame(rowQuantiles(
      yp[[p]],
      probs = c(0.5, 0.025, 0.975)
    ))
    dpqN[[p]] <- as.data.frame(rowQuantiles(
      ypN[[p]],
      probs = c(0.5, 0.025, 0.975)
    ))
    dpqP[[p]] <- as.data.frame(rowQuantiles(
      ypP[[p]],
      probs = c(0.5, 0.025, 0.975)
    ))
    dpq[[p]]$titer <- dpqN[[p]]$titer <- dpqP[[p]]$titer <- titer
    dpq[[p]]$pathogen <- dpqN[[p]]$pathogen <- dpqP[[p]]$pathogen <- pathogens[
      p
    ]
  }
  dpq <- do.call('rbind', dpq)
  dpqN <- do.call('rbind', dpqN)
  dpqP <- do.call('rbind', dpqP)
  colnames(dpq)[1:3] <- colnames(dpqN)[1:3] <- colnames(dpqP)[1:3] <- c(
    'med',
    'ciL',
    'ciU'
  )

  # compile data for plotting
  dta <- as.data.frame(data$y)
  colnames(dta) <- pathogens
  dta <- tidyr::gather(dta, key = 'pathogen', value = 't')

  # overall fit
  fitD <- ggplot() +
    geom_histogram(
      data = dta,
      aes(t, y = ..density..),
      bins = 150,
      fill = 'grey80',
      col = 'grey70'
    ) +
    theme_minimal() +
    theme(text = element_text(size = 18)) + #xlim(-1.5,3.5)+
    geom_line(data = dpq, aes(titer, med), col = 'springgreen4') +
    facet_wrap(~pathogen, scales = 'free_y') +
    xlab('titer') +
    geom_ribbon(
      data = dpq,
      aes(x = titer, y = med, ymin = ciL, ymax = ciU),
      fill = 'springgreen3',
      alpha = 0.4
    )

  # pos-neg fit
  fitDPN <- ggplot() +
    geom_histogram(
      data = dta,
      aes(t, y = ..density..),
      bins = 150,
      fill = 'grey80',
      col = 'grey70'
    ) +
    theme_minimal() +
    theme(text = element_text(size = 18)) + #xlim(-1.5,3.5)+
    geom_line(data = dpqN, aes(titer, med), col = 'mediumblue') +
    geom_line(data = dpqP, aes(titer, med), col = 'violetred') +
    facet_wrap(~pathogen, scales = 'free_y') +
    xlab('titer') +
    geom_ribbon(
      data = dpqN,
      aes(x = titer, y = med, ymin = ciL, ymax = ciU),
      fill = 'mediumblue',
      alpha = 0.3
    ) +
    geom_ribbon(
      data = dpqP,
      aes(x = titer, y = med, ymin = ciL, ymax = ciU),
      fill = 'violetred',
      alpha = 0.5
    )

  # return plots
  return(list(fit = fitD, fitPN = fitDPN))
}


#----- Plot posterior estimates of params
get_posterior = function(
  chains,
  data,
  pathogens,
  real_pars = NULL,
  real_CR = NULL
) {
  phi <- extract_phi(chains, data, pathogens)
  sero <- extract_sero(chains, data, pathogens)
  mu <- extract_mu(chains, data, pathogens)
  sds <- extract_sds(chains, data)

  mu0 = mu$mus0 %>%
    filter(pos == "neg") %>%
    mutate(
      par_label = " Mean negative\ntiter (Mu 0)",
      par = paste0("mu0 ", pathogens),
      idx = pathogens
    ) %>%
    select(par, par_label, med, ciL, ciU, idx)

  mu1 = mu$mus1 %>%
    filter(pos == antigen) %>%
    mutate(
      par_label = 'Mean positive\ntiter (Mu 1)',
      par = paste0("mu1 ", present),
      idx = present
    ) %>%
    select(par, par_label, med, ciL, ciU, idx)

  sdss = sds %>%
    mutate(
      par_label = "Std. dev.",
      par = c("Sd0", "Sd1"),
      idx = c("Sd0", "Sd1")
    ) %>%
    select(par, par_label, med, ciL, ciU, idx)

  phis = phi$phi %>%
    mutate(
      par_label = "Cross-reactivity (Phi)",
      par = paste0("Phi ", pos, " to ", neg),
      idx = paste0(pos, " to\n", neg)
    ) %>%
    select(par, par_label, med, ciL, ciU, idx)

  seros = sero %>%
    mutate(
      par_label = "Prevalence",
      par = paste0("Prev ", present),
      idx = present
    ) %>%
    select(par, par_label, med, ciL, ciU, idx)

  if (!is.null(real_pars)) {
    mu0 = mu0 %>% mutate(real = real_pars$true[real_pars$pars == "mu0"])
    mu1 = mu1 %>%
      mutate(
        real = real_pars$true[
          real_pars$pars == "mu1" & real_pars$pathogen %in% c("A", "B")
        ] +
          real_pars$true[
            real_pars$pars == "mu0" & real_pars$pathogen %in% c("A", "B")
          ],
      )
    sdss = sdss %>%
      mutate(
        real = real_pars$true[
          real_pars$pars %in% c("sd0", "sd1") & real_pars$pathogen %in% c("A")
        ]
      )
    colnames(real_cr$CR) = rownames(real_cr$CR) = pathogens
    phis = phis %>%
      mutate(real = diag(real_cr$CR[phi$phi$pos, phi$phi$neg]))
    seros = seros %>%
      mutate(real = real_pars$true[real_pars$pars == "sero"][1:2])
  }

  tab = rbind(mu0, mu1, sdss, phis, seros)

  p = tab %>%
    mutate(idx = as.factor(idx)) %>%
    ggplot() +
    geom_pointrange(aes(x = idx, y = med, ymin = ciL, ymax = ciU)) +
    facet_wrap(~par_label, scales = "free") +
    xlab("Index") +
    ylab("") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  tab = tab %>% select(-c(par_label, idx))

  if (!is.null(real_pars)) {
    p = p +
      geom_point(
        aes(x = idx, y = real),
        col = 2,
        size = 3,
        position = position_nudge(x = 0.2)
      )
    tab = tab %>% relocate(par, real)
  }

  return(list(table = tab, plot = p))
}


#----- BREAD data formatting
format_data <- function(
  df,
  present,
  nonpres,
  ageG = 0, # 1 if age groups included
  locID = 0, # 1 if location IDs included
  log_tranform = T
) {
  # create list for inputs
  data <- list()

  # antibody measurement data (on log scale)
  if (log_tranform) {
    data$y <- cbind(log(df[, c(present, nonpres)]))
  } else {
    data$y <- cbind((df[, c(present, nonpres)]))
  }

  # index for present/absent pathogens
  data$pres <- c(rep(1, length(present)), rep(0, length(nonpres)))

  # N individuals in study population
  data$N <- nrow(data$y)

  # N pathogens
  data$nP <- ncol(data$y)

  # N present pathogens
  data$nPp <- sum(data$pres)

  # Covar selection
  if (ageG == 1 & locID == 1) {
    # both age and location

    # Age group index
    data$ageG <- df$ageG
    # Location index
    data$loc <- df$locID
    # N age groups
    data$nA <- length(unique(df$ageG)) #7
    # N locations
    data$nL <- length(unique(df$locID)) #5
    # N individuals per location
    data$NperL <- as.vector(table(df$locID))
    # N individuals per location & age
    data$NperLA <- table(df$ageG, df$locID)
    # Proportion of study population by age group per location
    data$ageProp <- as.matrix(table(df$locID, df$ageG) / data$NperL)
  } else if (ageG == 0 & locID == 0) {
    # no strata

    data$ageG <- rep(1, data$N)
    data$loc <- rep(1, data$N)
    data$nA <- 1
    data$nL <- 1
    data$NperL <- as.matrix(data$N)
    data$NperLA <- as.matrix(data$N)
    data$ageProp <- as.matrix(1)
  } else if (ageG == 1 & locID == 0) {
    # age group only

    data$ageG <- df$ageG
    data$loc <- rep(1, data$N)
    data$nA <- length(unique(df$ageG))
    data$nL <- 1
    data$NperL <- as.matrix(data$N)
    data$NperLA <- as.matrix(table(df$ageG))
    data$ageProp <- as.matrix(table(df$ageG) / data$N)
  } else if (ageG == 0 & locID == 1) {
    # location only

    data$ageG <- rep(1, data$N)
    data$loc <- df$locID
    data$nA <- 1
    data$nL <- length(unique(data$locID))
    data$NperL <- as.vector(table(df$locID))
    data$NperLA <- as.matrix(table(df$ageG))
    data$ageProp <- as.matrix(1)
  }
  # Matrix of infection status combinations
  data$infM <- inf_matrix(data$nP, pres = data$pres)

  # N possible infection statuses
  data$nC <- nrow(data$infM) # N status combinations

  # N positive pathogens per infection status
  npos <- rowSums(data$infM)

  # Compute indices for which matrix cells are negative (wneg) or positive (wpos)
  wpos <- matrix(0, ncol = data$nP, nrow = data$nC)
  wneg <- matrix(0, ncol = data$nP, nrow = data$nC)
  for (c in 1:nrow(data$infM)) {
    for (p in 1:data$nP) {
      if (npos[c] > 0) {
        wpos[c, 1:npos[c]] <- which(data$infM[c, ] == 1)
      }
      if (npos[c] < data$nP) {
        wneg[c, 1:(data$nP - npos[c])] <- which(data$infM[c, ] == 0)
      }
    }
  }
  data$npos <- npos
  data$wpos <- wpos
  data$wneg <- wneg

  return(data)
}


#######   Get simulate data ########
MultiSeroSimulate = function(nP, nPp, N, nsims) {
  # Parameters
  # nsims <- 50 # N simulations
  # nP <- 3 # N pathogens
  # nPp <- 2 # N present pathogens (has to be <= nP)
  # N <- 1500 # N individuals

  pres <- c(rep(1, nPp), rep(0, nP - nPp)) # Pathogens A and B circulating, C not circulating

  # Recalculate number of parameters
  # - sero: nPp values (one per present pathogen)
  # - mu0: 1 value (shared across all pathogens)
  # - mu1: nPp values (one per present pathogen)
  # - sd0: 1 value (shared across all pathogens)
  # - sd1: 1 value (shared across all pathogens)
  # - phi: nPp * (nP - 1) cross-reactivity values
  # - corr00: 1 value
  # - corr11: 1 value

  npars <- nPp + 1 + nPp + 1 + 1 + (nPp * (nP - 1)) + 2
  # That's: 2 + 1 + 2 + 1 + 1 + 4 + 2 = 13 parameters

  params <- c(
    rep('sero', nPp), # 2 seroprevalence values
    'mu0', # 1 baseline mean (shared)
    rep('mu1', nPp), # 2 infection increases
    'sd0', # 1 baseline SD (shared)
    'sd1', # 1 infection SD (shared)
    rep('phi', nPp * (nP - 1)), # 4 cross-reactivity values (2 pathogens × 2 others each)
    'corr00', # 1 correlation
    'corr11' # 1 correlation
  )

  # Sample pars from Latin Hypercube
  pars <- randomLHS(nsims, npars)
  pars <- round(pars, 2)

  # Transform pars to desired scales
  parsT <- pars

  # Column indices
  col_sero <- 1:nPp
  col_mu0 <- (nPp + 1):(nPp + nP)
  col_mu1 <- (nPp + nP + 1):(nPp * 2 + nP)
  col_sd0 <- nPp * 2 + nP + 1
  col_sd1 <- nPp * 2 + nP + 2
  col_phi <- (nPp * 2 + nP + 3):(nPp * 2 + nP + 2 + nPp * (nP - 1))
  col_corr00 <- npars - 1
  col_corr11 <- npars

  # Transform to desired ranges
  parsT[, col_mu0] <- qunif(pars[, col_mu0], min = 0, max = 2)
  parsT[, col_mu1] <- qunif(pars[, col_mu1], min = 4, max = 6)
  parsT[, col_sd0] <- qunif(pars[, col_sd0], min = 0.75, max = 1.25)
  parsT[, col_sd1] <- qunif(pars[, col_sd1], min = 0.75, max = 1.25)
  parsT[, col_phi] <- qunif(pars[, col_phi], min = 0, max = 0.30)
  parsT[, col_sero] <- qunif(pars[, col_sero], min = 0.3, max = 0.70)
  parsT[, col_corr11] <- 0

  colnames(parsT) <- params

  #--- Simulate
  sims <- list()
  truepars <- truepars2 <- list()

  for (i in 1:nrow(parsT)) {
    # Extract seroprevalence (for present pathogens only)
    sero_present <- parsT[i, col_sero]

    # Create full sero vector (0 for non-circulating pathogen)
    sero <- numeric(nP)
    sero[pres == 1] <- sero_present

    # Extract shared baseline parameters (replicate for all pathogens)
    mu0 <- parsT[i, col_mu0]
    sd0 <- rep(parsT[i, col_sd0], nP)
    sd1_value <- parsT[i, col_sd1]
    sd1 <- rep(sd1_value, nP)

    # Extract mu1 (for present pathogens, expand to all)
    mu1_present <- parsT[i, col_mu1]
    mu1 <- numeric(nP)
    mu1[pres == 1] <- mu1_present

    # Build cross-reactivity matrix
    CR <- matrix(0, ncol = nP, nrow = nP)
    phi_vals <- parsT[i, col_phi]

    ind <- 1
    for (p in which(pres == 1)) {
      # Only for present pathogens
      for (p2 in 1:nP) {
        if (p != p2) {
          CR[p, p2] <- phi_vals[ind]
          ind <- ind + 1
        }
      }
    }

    # Extract correlations
    corr00 <- parsT[i, col_corr00]
    corr11 <- parsT[i, col_corr11]

    # Simulate
    sims[[i]] <- sim_multisero(
      N = N,
      nP = nP,
      sero = sero,
      pres = pres,
      mu0 = mu0,
      mu1 = mu1,
      sd0 = sd0,
      sd1 = sd1,
      CR = CR,
      corr00 = corr00,
      corr11 = corr11
    )

    # True params
    truepars[[i]] <- data.frame(
      pars = c(
        rep('sero', nP),
        rep('mu0', nP),
        rep('mu1', nP),
        rep('sd0', nP),
        rep('sd1', nP)
      ),
      pathogen = LETTERS[seq(1, nP)],
      present = rep(pres, 5),
      true = c(sero, mu0, mu1, sd0, sd1)
    )
    truepars2[[i]] <- list(CR = CR, r00 = corr00, r11 = corr11)
  }

  sims = lapply(sims, function(x) {
    temp = x$titers
    colnames(temp)[1:nP] = LETTERS[seq(1, nP)]
    x$titers = temp
    return(x)
  })

  return(list(sim_data = sims, truepars = truepars, trueCR = truepars2))
}


#----- Generate infection status matrix
inf_matrix <- function(N_pathogen, pres = rep(1, N_pathogen)) {
  # list of possible outcomes for each pathogen
  combos <- list()
  for (c in 1:N_pathogen) {
    combos[[c]] <- c(0, 1)
  }

  # matrix of all possible infection status combinations
  m <- expand.grid(combos)
  colnames(m) <- letters[1:N_pathogen]

  # remove positives of absent pathogens
  if (sum(pres) < N_pathogen) {
    for (abs in which(pres == 0)) {
      m <- m[m[, abs] == 0, ]
    }
  }

  return(m)
}


#----- Simulate multivariate gaussian mixtures
sim_multisero <- function(
  N,
  nP,
  sero,
  pres,
  mu0,
  mu1,
  sd0,
  sd1,
  CR,
  corr00,
  corr11
) {
  #--- infection status matrix
  infM <- inf_matrix(nP, pres = pres)

  #--- gaussian weights
  w <- infM
  for (p in 1:nP) {
    for (c in 1:nrow(infM)) {
      if (infM[c, p] == 0) {
        w[c, p] <- 1 - sero[p]
      } else {
        w[c, p] <- sero[p]
      }
    }
  }
  ws <- rowProds(as.matrix(w))
  Nws <- round(N * ws)

  #--- gaussian means & covariance matrices
  mus <- infM
  sigs <- infM
  for (c in 1:nrow(infM)) {
    if (sum(infM[c, ]) == 0) {
      # neg to all

      mus[c, ] <- mu0
      sigs[c, ] <- sd0
    } else if (sum(infM[c, ]) == 1) {
      # pos to just 1

      # which are pos/neg
      wp <- which(infM[c, ] == 1)
      wn <- which(infM[c, ] == 0)

      # pos
      mus[c, wp] <- mu0[wp] + mu1[wp]
      sigs[c, wp] <- sd1[wp]

      # neg
      for (n in wn) {
        mus[c, n] <- mu0[n] + mu1[wp] * CR[wp, n]
      }
      for (n in wn) {
        sigs[c, n] <- sqrt(
          sd0[n]^2 + (sd1[wp] * CR[wp, n])^2
        )
      }
    } else if (sum(infM[c, ]) > 1) {
      # pos to > 1

      # which are pos/neg
      wp <- which(infM[c, ] == 1)
      npos <- sum(infM[c, ])
      wn <- which(infM[c, ] == 0)

      # pos
      for (p in wp) {
        mus[c, p] <- mu0[p] + mu1[p]
      }
      for (p in wp) {
        sigs[c, p] <- sd1[p]
      }

      # neg
      for (n in wn) {
        mus[c, n] <- mu0[n]
        for (p in wp) {
          mus[c, n] <- mus[c, n] +
            mu1[p] * CR[p, n]
        }
        tmp <- rep(0, ncol(infM))
        tmp2 <- rep(0, ((npos * (npos - 1)) / 2))
        for (p in wp) {
          tmp[p] <- (sd1[p] * CR[p, n])^2
        }
        iy <- 1
        for (j in 1:(npos - 1)) {
          for (k in (j + 1):npos) {
            tmp2[iy] <- 2 *
              CR[wp[j], n] *
              CR[wp[k], n] *
              corr11 *
              sd1[wp[j]] *
              sd1[wp[k]]
            iy <- iy + 1
          }
        }
        sigs[c, n] <- sqrt(
          sd0[n]^2 + sum(tmp) + sum(tmp2)
        )
      }
    }
  }

  #--- covariance matrices
  covs <- list()
  Qcovs <- list()
  for (c in 1:nrow(infM)) {
    covs[[c]] <- matrix(0, ncol = nP, nrow = nP)
  }
  for (c in 1:nrow(infM)) {
    for (p in 1:nP) {
      for (p2 in p:nP) {
        if (p == p2) {
          # gaussian variance

          covs[[c]][p, p2] <- sigs[c, p]^2
        } else if (sum(infM[c, ]) == 0) {
          # negative to all pathogens

          covs[[c]][p, p2] <- covs[[c]][
            p2,
            p
          ] <- corr00 * sigs[c, p] * sigs[c, p2]
        } else if (infM[c, p] + infM[c, p2] == 0) {
          # negative to both where others are pos

          npos <- sum(infM[c, ])
          wp <- which(infM[c, ] == 1)

          if (npos == 1) {
            cv00 <- corr00 *
              sd0[p] *
              sd0[p2] +
              CR[wp, p] *
                CR[wp, p2] *
                sd1[wp]^2
          } else {
            vrs <- rep(0, npos)
            cvs <- rep(0, (npos^2 - npos))
            for (x in 1:npos) {
              vrs[x] = CR[wp[x], p] *
                CR[wp[x], p2] *
                sd1[wp[x]]^2
            }
            ind <- 1
            for (j in 1:(npos - 1)) {
              for (k in (j +
                1):npos) {
                if (j == k) {
                  NULL
                } else {
                  cvs[
                    ind
                  ] <- CR[
                    wp[
                      j
                    ],
                    p
                  ] *
                    CR[
                      wp[
                        k
                      ],
                      p2
                    ] *
                    corr11 *
                    sd1[wp[
                      j
                    ]] *
                    sd1[wp[
                      k
                    ]]
                  cvs[
                    ind +
                      1
                  ] <- CR[
                    wp[
                      k
                    ],
                    p
                  ] *
                    CR[
                      wp[
                        j
                      ],
                      p2
                    ] *
                    corr11 *
                    sd1[wp[
                      j
                    ]] *
                    sd1[wp[
                      k
                    ]]
                  ind <- ind +
                    2
                }
              }
            }
            cv00 <- corr00 *
              sd0[p] *
              sd0[p2] +
              sum(vrs) +
              sum(cvs)
          }
          covs[[c]][p, p2] <- covs[[c]][
            p2,
            p
          ] <- cv00
        } else if (infM[c, p] + infM[c, p2] == 1) {
          # negative & positive

          npos <- sum(infM[c, ])
          wp <- which(infM[c, ] == 1)
          pos <- c(p, p2)[which(
            infM[c, c(p, p2)] == 1
          )]
          neg <- c(p, p2)[which(
            infM[c, c(p, p2)] == 0
          )]
          wp <- wp[!wp == pos]

          tmp <- rep(0, npos)
          tmp[1] <- CR[pos, neg] * sd1[pos]^2
          if (npos > 1) {
            for (m in 1:(npos - 1)) {
              tmp[m + 1] <- CR[
                wp[m],
                neg
              ] *
                corr11 *
                sd1[pos] *
                sd1[wp[m]]
            }
          }

          covs[[c]][p, p2] <- covs[[c]][
            p2,
            p
          ] <- sum(tmp)
        } else if (infM[c, p] + infM[c, p2] == 2) {
          # positive to both

          covs[[c]][p, p2] <- covs[[c]][
            p2,
            p
          ] <- corr11 * sigs[c, p] * sigs[c, p2]
        }
      }
    }
  }
  for (c in 1:nrow(infM)) {
    Qcovs[[c]] <- is.positive.semi.definite(covs[[c]])
  }

  #--- simulate multivariate gaussians
  gs <- list()
  for (c in 1:nrow(infM)) {
    if (Nws[c] > 0) {
      gs[[c]] <- as.data.frame(rmvnorm(
        Nws[c],
        as.numeric(mus[c, ]),
        covs[[c]]
      ))
      gs[[c]]$status <- c
      gs[[c]][, LETTERS[seq(1, nP)]] <- NA
      for (p in 1:nP) {
        gs[[c]][, nP + 1 + p] <- infM[c, p]
      }
      gs[[c]]$stat <- paste(
        '{',
        paste(infM[c, ], collapse = ','),
        '}',
        sep = ''
      )
    }
  }
  gs <- do.call('rbind', gs)

  #--- return data
  return(list(
    titers = gs,
    mus = mus,
    sigs = sigs,
    covs = covs,
    Qcovs = Qcovs
  ))
}


relabel_chains <- function(chains, pathogens, present) {
  # general pars
  label_map <- c(
    "sd0" = "Std. dev. (negative)",
    "sd1" = "Std. dev. (positive)",
    "rho00" = "Rho",
    "lp__" = "Likelihood"
  )
  # mus - pathogen specific
  label_map <- c(
    label_map,
    setNames(
      paste0("Seroprev: ", pathogens),
      paste0("seroAll[", seq_along(pathogens), "]")
    ),
    setNames(
      paste0("Mu0: ", pathogens),
      paste0("mu0[", seq_along(pathogens), "]")
    ),
    setNames(paste0("Mu1:", present), paste0("mu1[", seq_along(present), "]"))
  )
  # phi labels - present * pathogens
  phi_map <- c()
  phi_idx <- 1
  for (p in seq_along(present)) {
    for (p2 in seq_along(pathogens)) {
      if (p2 != p) {
        phi_name <- paste0("phi[", phi_idx, "]")
        phi_map[phi_name] <- paste0("Phi: ", present[p], " to\n", pathogens[p2])
        phi_idx <- phi_idx + 1
      }
    }
  }

  label_map <- c(
    label_map,
    phi_map[intersect(names(phi_map), colnames(chains))]
  )

  pars <- names(label_map)
  trace_df <- chains |>
    select(.chain, .iteration, .draw, all_of(pars)) |>
    rename_with(~ unname(label_map[.x]), .cols = all_of(pars))

  return(trace_df)
}
