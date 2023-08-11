#' Offspring Distribution
#'
#' This function generates distribution of wild-type, mutant, and overall
#' offspring with a given R0, mutation rate, and a maximum number of offspring
#' you want to calculate for.
#'
#' @param R0 Reproduction number, a positive number
#' @param mu Mutation rate, a positive number
#' @param maxOff Maximum offspring number being calculated
#' @return A data frame of offspring size, type, and probabilities
#' @examples
#' pmf.alloff(1.2, 0.4, 5)
#' pmf.alloff(3.5, 1.6, 8)
#' @export
pmf.alloff <- function(R0, mu, maxOff) {
  k_overall <- 1
  k_wt <- exp(-mu)
  k_mu <- (1 - exp(-mu))
  p <- 1 / (R0 + 1)
  pmf <- data.frame(
    x = 0:maxOff,
    overall = dnbinom(0:maxOff, k_overall, p),
    wildtype = dnbinom(0:maxOff, k_wt, p),
    mutant = dnbinom(0:maxOff, k_mu, p)
  )
  pmf <- melt(pmf, id.vars = "x")
  names(pmf) <- c("size", "type", "prob")
  return(pmf)
}

#' Wild-type Establishment Probability
#'
#' This function calculates the probability of establishment of wild-type
#' offspring of a lineage starting with certain parameters.
#'
#' @param n Initial population size, a positive integer
#' @param R0 Reproduction number, a positive number
#' @param mu Mutation rate, a positive number
#' @return A number representing the probability
#' @examples
#' P_wt_est(2, 1.2, 0.4)
#' P_wt_est(1, 3.5, 1.6)
#' @export
P_wt_est <- function(n, R0, mu) {
  lhs <- function(p) {
    p
  }
  rhs <- function(p, R0, mu) {
    rhs <- 1 / (1 + (R0 * exp(-mu) * (1 - p)) / (exp(-mu)))^exp(-mu)
    return(rhs)
  }
  uni <- try(uniroot(function(p) rhs(p, R0, mu) - lhs(p), c(0, 0.9), extendInt = "yes"), silent = TRUE)
  uni <- if (inherits(uni, "try-error")) 1 else uni$root
  if (uni < 1 & uni > 0) {
    ext_1 <- uni
  } else {
    ext_1 <- 1
  }
  est_n <- 1 - ext_1^n
  return(est_n)
}

#' Final Size Probability
#'
#' This function calculates the probability of getting a specific final size y
#' with initial population size of 1 and some certain parameters.
#'
#' @param y Final size, a positive integer
#' @param R0 Reproduction number, a positive number
#' @param mu Mutation rate, a positive number
#' @return A number representing the probability
#' @examples
#' FS(3, 1.2, 0.4)
#' FS(6, 3.5, 1.6)
#' @export
FS <- function(y, R0, mu) {
  k <- exp(-mu)
  meanR0 <- R0 * k
  if (y == 1) {
    prob <- 1 / (1 + meanR0 / k)^k
  } else if (y > 1) {
    prod <- 1
    for (j in 0:(y - 2)) {
      prod <- prod * (j / k + y)
    }
    prob <- prod / factorial(y) * (k / (meanR0 + k))^(k * y) * (meanR0 * k / (meanR0 + k))^(y - 1)
  } else {
    prob <- 0
  }
  return(prob)
}

#' Final Size Distribution
#'
#' This function calculates the final size distribution with initial population
#' size of 1 and some certain parameters.
#'
#' @param R0 Reproduction number, a positive number
#' @param mu Mutation rate, a positive number
#' @param maxFS Maximum number of final sizes being calculated
#' @return A data frame of final sizes and corresponding probabilities
#' @examples
#' pmf.FS(1.2, 0.4, 60)
#' pmf.FS(3.5, 1.6, 120)
#' @export
pmf.FS <- function(R0, mu, maxFS) {
  prob <- rep(0, maxFS)
  for (i in 1:maxFS) {
    prob[i] <- FS(i, R0, mu)
  }
  names(prob) <- c(1:maxFS)
  df <- data.frame(FinalSize = names(prob), prob)
  df$FinalSize <- as.numeric(df$FinalSize)
  return(df)
}

#' Final Size Probability with Arbitrary Initial Population Size
#'
#' This function calculates the probability of getting a specific final size y
#' with initial population size of N and some certain parameters.
#'
#' @param n Initial population size, a positive integer
#' @param y Final size, a positive integer
#' @param R0 Reproduction number, a positive number
#' @param mu Mutation rate, a positive number
#' @return A number representing the probability
#' @examples
#' FSn(2, 3, 1.2, 0.4)
#' FSn(3, 6, 3.5, 1.6)
#' @export
FSn <- function(n, y, R0, mu) {
  prob <- 0
  if (y >= n) {
    if (n == 1) {
      prob <- FS(y, R0, mu)
    } else if (n > 1) {
      for (i in 1:(y - n + 1)) {
        prob <- prob + FS(i, R0, mu) * FSn(n - 1, y - i, R0, mu)
      }
    }
  }
  return(prob)
}

#' Final Size Distribution with Arbitrary Initial Population Size
#'
#' This function calculates the final size distribution with initial population
#' size of N and some certain parameters.
#'
#' @param n Initial population size, a positive integer
#' @param R0 Reproduction number, a positive number
#' @param mu Mutation rate, a positive number
#' @param maxFS Maximum number of final sizes being calculated
#' @return A data frame of final sizes and corresponding probabilities
#' @examples
#' pmf.FSn(2, 1.2, 0.4, 60)
#' pmf.FSn(3, 3.5, 1.6, 120)
#' @export
pmf.FSn <- function(n, R0, mu, maxFS) {
  df_n1 <- pmf.FS(R0, mu, maxFS)
  prob <- rep(0, maxFS)
  if (n == 1) {
    df <- df_n1
  } else if (n > 1) {
    df_n_one_less <- pmf.FSn(n - 1, R0, mu, maxFS)
    for (i in (1:(maxFS - n))) {
      for (j in ((n - 1):(maxFS - i))) {
        prob[i + j] <- prob[i + j] + df_n1[i, 2] * df_n_one_less[j, 2]
      }
    }
    names(prob) <- seq(1, maxFS)
    df <- data.frame(FinalSize = names(prob), prob)
  }
  pi <- max((1 - P_wt_est(n, R0, mu)), sum(df$prob))
  df$prob_normal <- df$prob / pi
  df$FinalSize <- as.numeric(df$FinalSize)
  return(df)
}

#' Mutant Lineages Generated Probability
#'
#' This function calculates the probability of generating a specific number of
#' mutant lineages.
#'
#' @param n Initial population size, a positive integer
#' @param m Number of mutant lineages generated, a positive integer
#' @param R0 Reproduction number, a positive number
#' @param mu Mutation rate, a positive number
#' @param maxFS Maximum number of final sizes being calculated
#' @return A number representing the probability
#' @examples
#' P_mu_n(2, 1, 1.2, 0.4, 60)
#' P_mu_n(3, 4, 3.5, 1.6, 120)
#' @export
P_mu_n <- function(n, m, R0, mu, maxFS) {
  prob <- 0
  df <- pmf.FSn(n, R0, mu, maxFS)
  for (i in 1:maxFS) {
    k <- (1 - exp(-mu)) * i
    p <- 1 / (R0 + 1)
    NBmu <- dnbinom(m, k, p)
    prob <- prob + df[i, 3] * NBmu
  }
  return(prob)
}

#' Mutant Lineages Generated Distribution
#'
#' This function calculates the mutant lineages generated distribution.
#'
#' @param n Initial population size, a positive integer
#' @param R0 Reproduction number, a positive number
#' @param mu Mutation rate, a positive number
#' @param maxMuGen Maximum number of mutant lineages being calculated
#' @param maxFS Maximum number of final sizes being calculated
#' @return A data frame of mutant lineages generated and corresponding probabilities
#' @examples
#' pmf.P_mu_n(2, 1.2, 0.4, 60, 60)
#' pmf.P_mu_n(3, 3.5, 1.6, 120, 130)
#' @export
pmf.P_mu_n <- function(n, R0, mu, maxMuGen, maxFS) {
  prob <- rep(0, maxMuGen + 1)
  df_FS <- pmf.FSn(n, R0, mu, maxFS)
  for (i in 0:maxMuGen) {
    prob_mugen <- 0
    for (j in 1:maxFS) {
      k <- (1 - exp(-mu)) * j
      p <- 1 / (R0 + 1)
      NBmu <- dnbinom(i, k, p)
      prob_mugen <- prob_mugen + df_FS[j, 3] * NBmu
    }
    prob[i + 1] <- prob[i + 1] + prob_mugen
  }
  names(prob) <- c(0:maxMuGen)
  df <- data.frame(MuLinGen = names(prob), prob)
  df$MuLinGen <- as.numeric(df$MuLinGen)
  return(df)
}

#' Mutant Lineages Established Probability
#'
#' This function calculates the probability of a specific number of mutant
#' lineages established.
#'
#' @param n Initial population size, a positive integer
#' @param m Number of mutant lineages established, a positive integer
#' @param R0 Reproduction number, a positive number
#' @param mu Mutation rate, a positive number
#' @param maxMuGen Maximum number of mutant lineages being calculated
#' @param maxFS Maximum number of final sizes being calculated
#' @return A number representing the probability
#' @examples
#' P_mu_est(2, 1, 1.2, 0.4, 60, 60)
#' P_mu_est(3, 4, 3.5, 1.6, 120, 130)
#' @export
P_mu_est <- function(n, m, R0, mu, maxMuGen, maxFS) {
  p_est <- 1 - 1 / R0
  prob <- rep(0, 1)
  pmf_mu_gen <- pmf.P_mu_n(n, R0, mu, maxMuGen, maxFS)
  for (i in m:maxMuGen) {
    p_m <- pmf_mu_gen[i + 1, 2]
    prob <- prob + dbinom(m, i, p_est) * p_m
  }
  return(prob)
}

#' Mutant Lineages Established Distribution
#'
#' This function calculates the mutant lineages established distribution.
#'
#' @param n Initial population size, a positive integer
#' @param R0 Reproduction number, a positive number
#' @param mu Mutation rate, a positive number
#' @param maxMuGen Maximum number of mutant lineages being calculated
#' @param maxFS Maximum number of final sizes being calculated
#' @return A data frame of mutant lineages established and corresponding probabilities
#' @examples
#' pmf.P_mu_est(2, 1.2, 0.4, 60, 60)
#' pmf.P_mu_est(3, 3.5, 1.6, 120, 130)
#' @export
pmf.P_mu_est <- function(n, R0, mu, maxMuGen, maxFS) {
  p_est <- 1 - 1 / R0
  df_mu_gen <- pmf.P_mu_n(n, R0, mu, maxMuGen, maxFS)
  prob <- rep(0, maxMuGen + 1)
  for (j in 0:maxMuGen) {
    p_m <- df_mu_gen[j + 1, 2]
    for (i in 0:j) {
      prob[i + 1] <- prob[i + 1] + dbinom(i, j, p_est) * p_m
    }
  }
  names(prob) <- c(0:maxMuGen)
  df <- data.frame(MuLinEst = names(prob), prob)
  df$MuLinEst <- as.numeric(df$MuLinEst)
  return(df)
}

#' Analytical Probability of Extinction
#'
#' This function calculates the probability of the whole lineage goes extinction
#' under a geometric distribution assumption.
#'
#' @param n Initial population size, a positive integer
#' @param R0 Reproduction number, a positive number
#' @return A number representing the probability
#' @examples
#' Px_analytical(2, 1.2)
#' Px_analytical(3, 3.5)
#' @export
Px_analytical <- function(n, R0) {
  px <- 1 / R0^n
  return(px)
}

#' Probability of Extinction
#'
#' This function calculates the probability of the whole lineage goes extinction
#' through multiplying the probability that both the wild-type and the mutant
#' lineages goes extinction.
#'
#' @param n Initial population size, a positive integer
#' @param R0 Reproduction number, a positive number
#' @param mu Mutation rate, a positive number
#' @param maxMuGen Maximum number of mutant lineages being calculated
#' @param maxFS Maximum number of final sizes being calculated
#' @return A number representing the probability
#' @examples
#' Px(2, 1.2, 0.4, 60, 60)
#' Px(3, 3.5, 1.6, 120, 130)
#' @export
Px <- function(n, R0, mu, maxMuGen, maxFS) {
  p_wt_ext <- 1 - P_wt_est(n, R0, mu)
  p_mu_est_0 <- P_mu_est(n, 0, R0, mu, maxMuGen, maxFS)
  prob <- p_mu_est_0 * p_wt_ext
  return(prob)
}

#' Probability of Getting 0 Clonal Mutation
#'
#' This function calculates the probability of getting 0 clonal mutation.
#'
#' @param n Initial population size, a positive integer
#' @param R0 Reproduction number, a positive number
#' @param mu Mutation rate, a positive number
#' @param maxMuGen Maximum number of mutant lineages being calculated
#' @param maxFS Maximum number of final sizes being calculated
#' @return A number representing the probability
#' @examples
#' P0(2, 1.2, 0.4, 60, 60)
#' P0(3, 3.5, 1.6, 120, 130)
#' @export
P0 <- function(n, R0, mu, maxMuGen, maxFS) {
  S_inf <- P_wt_est(n, R0, mu)
  S_2plus <- (1 - P_mu_est(n, 0, R0, mu, maxMuGen, maxFS) - P_mu_est(n, 1, R0, mu, maxMuGen, maxFS)) * (1 - S_inf)
  prob <- S_inf + S_2plus
  return(prob)
}

#' Probability of Getting More Than 1 Clonal Mutation
#'
#' This function calculates the probability of getting more than 1 clonal
#' mutations.
#'
#' @param n Initial population size, a positive integer
#' @param R0 Reproduction number, a positive number
#' @param mu Mutation rate, a positive number
#' @param maxMuGen Maximum number of mutant lineages being calculated
#' @param maxFS Maximum number of final sizes being calculated
#' @return A number representing the probability
#' @examples
#' P1plus(2, 1.2, 0.4, 60, 60)
#' P1plus(3, 3.5, 1.6, 120, 130)
#' @export
P1plus <- function(n, R0, mu, maxMuGen, maxFS) {
  p_wt_ext <- 1 - P_wt_est(n, R0, mu)
  prob <- P_mu_est(n, 1, R0, mu, maxMuGen, maxFS) * p_wt_ext
  return(prob)
}

#' Round Probability of Getting Clonal Mutation
#'
#' This function calculates the round probability of getting a certain number
#' of clonal mutations. Notice that this is not the exact probability.
#'
#' @param n Initial population size, a positive integer
#' @param m Number of clonal mutation, a positive integer
#' @param R0 Reproduction number, a positive number
#' @param mu Mutation rate, a positive number
#' @param maxMuGen Maximum number of mutant lineages being calculated
#' @param maxFS Maximum number of final sizes being calculated
#' @return A number representing the probability
#' @examples
#' pmf.P_clonal_round(2, 3, 1.2, 0.4, 60, 60)
#' pmf.P_clonal_round(3, 6, 3.5, 1.6, 120, 130)
#' @export
P_clonal_round <- function(n, m, R0, mu, maxMuGen, maxFS) {
  prob <- P1plus(n, R0, mu, maxMuGen, maxFS) * dpois(m, mu) / (1 - dpois(0, mu))
  return(prob)
}

#' Round Probability Distribution of Getting Clonal Mutation
#'
#' This function calculates the round probability distribution of getting a
#' certain number of clonal mutations. Notice that this is not the exact
#' probability.
#'
#' @param n Initial population size, a positive integer
#' @param R0 Reproduction number, a positive number
#' @param mu Mutation rate, a positive number
#' @param maxMuGen Maximum number of mutant lineages being calculated
#' @param maxFS Maximum number of final sizes being calculated
#' @param maxClonal Maximum number of clonal mutations being calculated
#' @return A data frame of clonal mutations and corresponding probabilities
#' @examples
#' pmf.P_clonal_round(2, 1.2, 0.4, 60, 60, 5)
#' pmf.P_clonal_round(3, 3.5, 1.6, 120, 130, 5)
#' @export
pmf.P_clonal_round <- function(n, R0, mu, maxMuGen, maxFS, maxClonal) {
  p1plus <- P1plus(n, R0, mu, maxMuGen, maxFS)
  prob <- rep(0, maxClonal + 2)
  prob[1] <- 1 / R0^n
  prob[2] <- P0(n, R0, mu, maxMuGen, maxFS)
  if (maxClonal != 0) {
    for (i in 1:maxClonal) {
      prob[i + 2] <- p1plus * dpois(i, mu) / (1 - dpois(0, mu))
    }
  }
  names(prob) <- c("x", 0:maxClonal)
  df <- data.frame(ClonalMu = names(prob), prob)
  return(df)
}

#' Probability of Getting Clonal Mutation
#'
#' This function calculates the round probability of getting a certain number
#' of clonal mutations.
#'
#' @param n Initial population size, a positive integer
#' @param m Number of clonal mutation, a non-negative integer
#' @param R0 Reproduction number, a positive number
#' @param mu Mutation rate, a positive number
#' @param maxMuGen Maximum number of mutant lineages being calculated
#' @param maxFS Maximum number of final sizes being calculated
#' @return A number representing the probability
#' @examples
#' P_clonal(2, 4, 1.2, 0.4, 60, 60)
#' P_clonal(3, 0, 3.5, 1.6, 120, 130)
#' @export
P_clonal <- function(n, m, R0, mu, maxMuGen, maxFS) {
  df_clonal_round <- pmf.P_clonal_round(n, R0, mu, maxMuGen, maxFS, m)[-1, ]
  px <- 1 / R0
  prob <- 0
  if (m == 0) {
    prob <- df_clonal_round[1, 2]
  } else {
    for (i in 1:m) {
      prob <- prob + df_clonal_round[i + 1, 2] * P_clonal(1, m - i, R0, mu, maxMuGen, maxFS) / (1 - px)
    }
  }
  return(prob)
}

#' Probability Distribution of Getting Clonal Mutation
#'
#' This function calculates the probability distribution of getting a certain
#' number of clonal mutations. The resulting pmf will contain probability of
#' extinction Px.
#'
#' @param n Initial population size, a positive integer
#' @param R0 Reproduction number, a positive number
#' @param mu Mutation rate, a positive number
#' @param maxMuGen Maximum number of mutant lineages being calculated
#' @param maxFS Maximum number of final sizes being calculated
#' @param maxClonal Maximum number of clonal mutations being calculated
#' @return A data frame of clonal mutations and corresponding probabilities
#' @examples
#' pmf.P_clonal_withPx(2, 1.2, 0.4, 60, 60, 5)
#' pmf.P_clonal_withPx(3, 3.5, 1.6, 120, 130, 8)
#' @export
pmf.P_clonal_withPx <- function(n, R0, mu, maxMuGen, maxFS, maxClonal) {
  df_clonal_round <- pmf.P_clonal_round(n,R0,mu,maxMuGen,maxFS,maxClonal)[-1, ]
  px <- 1 / R0
  prob <- rep(0, maxClonal + 2)
  prob[1] <- 1 / R0^n
  prob[2] <- df_clonal_round[1, 2]
  if (maxClonal != 0) {
    df_clonal_max_one_less <- pmf.P_clonal_withPx(1, R0, mu, maxMuGen, maxFS, maxClonal - 1)
    for (m in 1:maxClonal) {
      for (i in 1:m) {
        prob[m+2] <- prob[m+2]+df_clonal_round[i+1,2]*df_clonal_max_one_less[m-i+2,2]/(1-px)
      }
    }
  }
  names(prob) <- c("x", 0:maxClonal)
  df <- data.frame(ClonalMu = names(prob), prob)
  return(df)
}

#' Probability Distribution of Getting Clonal Mutation Normalized
#'
#' This function calculates the probability distribution of getting a certain
#' number of clonal mutations. The resulting probabilities are normalized by
#' establishment, 1.e. 1-Px, and Px will not be contained in the final data
#' frame.
#'
#' @param n Initial population size, a positive integer
#' @param R0 Reproduction number, a positive number
#' @param mu Mutation rate, a positive number
#' @param maxMuGen Maximum number of mutant lineages being calculated
#' @param maxFS Maximum number of final sizes being calculated
#' @param maxClonal Maximum number of clonal mutations being calculated
#' @return A data frame of clonal mutations and corresponding probabilities
#' @examples
#' pmf.P_clonal(2, 1.2, 0.4, 60, 60, 5)
#' pmf.P_clonal(3, 3.5, 1.6, 120, 130, 8)
#' @export
pmf.P_clonal <- function(n, R0, mu, maxMuGen, maxFS, maxClonal) {
  df_clonal_round <- pmf.P_clonal_round(n, R0, mu, maxMuGen, maxFS, maxClonal)[-1, ]
  px <- 1 / R0
  prob <- rep(0, maxClonal + 2)
  prob[1] <- 1 / R0^n
  prob[2] <- df_clonal_round[1, 2]
  if (maxClonal != 0) {
    df_clonal_max_one_less <- pmf.P_clonal_withPx(1, R0, mu, maxMuGen, maxFS, maxClonal - 1)
    for (m in 1:maxClonal) {
      for (i in 1:m) {
        prob[m + 2] <- prob[m + 2] + df_clonal_round[i + 1, 2] * df_clonal_max_one_less[m - i + 2, 2] / (1 - px)
      }
    }
  }
  names(prob) <- c("x", 0:maxClonal)
  df <- data.frame(ClonalMu = names(prob), prob)
  df$prob <- df$prob / (1 - df$prob[1])
  df <- df[-1, ]
  return(df)
}

#' Clonal Mutations with Various Parameters
#'
#' This function calculates for the clonal mutation distribution
#' across pairs of parameters N and mu.
#'
#' @param n_values A sequence of initial population sizes
#' @param R0 Reproduction number, a positive number
#' @param mu_values A sequence of mutation rate
#' @param maxMuGen Maximum number of mutant lineages being calculated
#' @param maxFS Maximum number of final sizes being calculated
#' @param clonal Maximum number of clonal mutations being calculated
#' @return A list of clonal mutation pmf with different parameters
#' @examples
#' list_clonal(1:10, 1.6, seq(0.1, 2, by = 0.2), 50, 50, 5)
#' list_clonal(seq(2, 10, by = 2), 0.8, seq(1.2, 2.5, by = 0.1), 60, 120, 7)
#' @export
list_clonal <- function(n_values, R0, mu_values, maxMuGen, maxFS, clonal) {
  listClonal <- 1:length(mu_values) %>%
    mclapply(function(x) {
      1:length(n_values) %>%
        lapply(function(y) {
          pmf.P_clonal(n_values[y], R0, mu_values[x], maxMuGen, maxFS, clonal)
        })
    }, mc.cores = detectCores())
  names(listClonal) <- mu_values
  for (i in 1:length(mu_values)) {
    names(listClonal[[i]]) <- n_values
  }
  return(listClonal)
}

#' Clonal Mutation Probability with Various Parameters
#'
#' This function calculates the probability of getting a certain
#' number of clonal mutations across pairs of parameters N and mu.
#'
#' @param clonal The number of clonal mutations being calculated, a positive integer
#' @param listClonal Returned list from function list_clonal with same parameters
#' @param n_values A sequence of initial population sizes
#' @param R0 Reproduction number, a positive number
#' @param mu_values A sequence of mutation rate
#' @param maxMuGen Maximum number of mutant lineages being calculated
#' @param maxFS Maximum number of final sizes being calculated
#' @return A data frame of N, mu values and corresponding probability of getting the input clonal mutations
#' @examples
#' pmf.clonal_params(0, listClonal, 1:10, 1.6, seq(0.1, 2, by = 0.2), 50, 50)
#' pmf.clonal_params(2, listClonal, seq(2, 10, by = 2), 0.8, seq(1.2, 2.5, by = 0.1), 60, 120)
#' @export
pmf.clonal_params <- function(clonal, listClonal, n_values, R0, mu_values, maxMuGen, maxFS) {
  df <- expand.grid(n = paste0("", n_values), mu = paste0("", mu_values))
  df$prob <- 0

  for (N in n_values) {
    for (MU in mu_values) {
      muloc <- which(names(listClonal)== MU)
      nloc <- which(names(listClonal[[muloc]])== N)
      prob <- listClonal[[muloc]][[nloc]]$prob[clonal+1]
      df$prob[which(df$n == N & df$mu == MU)] <- prob
    }
  }

  df$n <- as.numeric(df$n)
  df$mu <- as.numeric(paste(df$mu))
  return(df)
}

#' Initial Population Size Distribution
#'
#' This function calculates the poisson Distribution of initial population size
#' N given a mean of lambdaand also the adjusted distribution normalized by at
#' least one of the initial viral particles established.
#'
#' @param lambda Mean of Poisson distribution, a positive number
#' @param R0 Reproduction number, a positive number
#' @param maxIni Maximum initial population size being calculated
#' @return A data frame of N, distribution type and corresponding probability of getting that number of initial population size
#' @examples
#' pmf.ini_size(2.1, 1.6, 10)
#' pmf.ini_size(0.8, 2.4, 6)
#' @export
pmf.ini_size <- function(lambda, R0, maxIni) {
  df <- data.frame(N = 0:maxIni, pois = dpois(0:maxIni, lambda))
  df$adj_prob <- 0
  for (n in 0:maxIni) {
    Pest <- 1 - 1 / R0^n
    df$adj_prob[n + 1] <- df$pois[n + 1] * Pest
  }
  df$adj_prob <- df$adj_prob / sum(df$adj_prob)
  df <- melt(df, id.vars = "N")
  names(df) <- c("N", "type", "prob")
  df$type <- factor(df$type)
  levels(df$type) <- c("poisson", "adjusted")
  return(df)
}

#' Lambda Log Likelihood
#'
#' This function calculates the log likelihood of getting the observed data
#' points in the provided data set given a certain lambda, R0, and mu. We
#' usually assume that R0 is a known value and aims to estimate lambda and mu.
#'
#' @param df Data set with first column named "ClonalMu" recording number of clonal mutations, and the second column named "freq" recording number of times the clonal mutation is observed
#' @param listClonal Returned list from function list_clonal with same parameters
#' @param lambda Mean of Poisson distribution, a positive number
#' @param R0 Reproduction number, a positive number
#' @param mu Mutation rate, a positive number
#' @param maxMuGen Maximum number of mutant lineages being calculated
#' @param maxFS Maximum number of final sizes being calculated
#' @param maxIni Maximum initial population size being calculated
#' @return A number representing the probability
#' @examples
#' LL_lambda(sample_data, listClonal, 2.1, 1.6, 0.4, 50, 50, 10)
#' LL_lambda(flu_data, listClonal, 0.2, 10, 0.7, 60, 70, 8)
#' @export
LL_lambda <- function(df, listClonal, lambda, R0, mu, maxMuGen, maxFS, maxIni) {
  df$prob <- df$freq / sum(df$freq)
  dfPn <- pmf.ini_size(lambda, R0, maxIni)
  dfPn <- dfPn[dfPn$type == "adjusted",]
  prob <- 1

  list <- mclapply(df$ClonalMu,
                   pmf.clonal_params,
                   listClonal = listClonal,
                   n_values = 1:maxIni,
                   R0 = R0,
                   mu_values = mu,
                   maxMuGen = maxMuGen,
                   maxFS = maxFS,
                   mc.cores = detectCores())
  names(list) <- df$ClonalMu

  for (i in df$ClonalMu) {
    loc <- which(names(list) == i)
    dfClonalN <- list[[loc]]
    PclonalLam <- 0
    for (n in 1:maxIni) {
      Pn <- dfPn$prob[dfPn$N == n]
      PclonalN <- dfClonalN$prob[dfClonalN$n == n]
      PclonalLam <- PclonalLam + Pn*PclonalN
    }
    k <- df$freq[df$ClonalMu == i]
    prob <- prob * PclonalLam^k
  }
  prob <- log(prob)
  return(prob)
}


#' Log Likelihood Heat Map
#'
#' This function calculates the log likelihood of getting the observed data
#' points in the provided data set across a set of parameters lambda, mean N,
#' transmission bottleneck Nb, and mu. Since the actual distribution of initial
#' population size N does not follow a perfect Poisson distribution but an
#' altered version of Poisson distribution, we are also interested in estimating
#' the most likely mean N, which might be slightly different from lambda. We can
#' also estimate the mean of transmission bottleneck size Nb by applying a
#' binomial distribution to the initial population size.
#'
#' @param df Data set with first column named "ClonalMu" recording number of clonal mutations, and the second column named "freq" recording number of times the clonal mutation is observed
#' @param listClonal Returned list from function list_clonal with same parameters
#' @param lambda_values A sequence of mean of Poisson distribution
#' @param R0 Reproduction number, a positive number
#' @param mu_values A sequence of mutation rate
#' @param maxMuGen Maximum number of mutant lineages being calculated
#' @param maxFS Maximum number of final sizes being calculated
#' @param maxIni Maximum initial population size being calculated
#' @return A data frame with lambda, mean N, mean Nb, mu, and corresponding probabilities of getting the set of parameters given data points in the input data set
#' @examples
#' LL_meanNb.df(sample_data, listClonal, seq(0.1, 7, by = 0.2), 1.6, seq(0.1, 2, by = 0.2), 50, 50, 10)
#' LL_meanNb.df(flu_data, listClonal, seq(0.1, 1.4, by = 0.1), 10, seq(1.2, 2.5, by = 0.1), 60, 70, 8)
#' @export
LL_meanNb.df <- function(df, listClonal, lambda_values, R0, mu_values, maxMuGen, maxFS, maxIni, maxNb) {
  prob <- 1:length(mu_values) %>%
    mclapply(function(x) {
      1:length(lambda_values) %>%
        lapply(function(y) {
          LL_lambda(df, listClonal, lambda_values[y],
                    R0, mu_values[x], maxMuGen, maxFS, maxIni)
        })
    }, mc.cores = detectCores())
  names(prob) <- mu_values
  for (i in 1:length(mu_values)) {
    names(prob[[i]]) <- lambda_values
  }

  LLdf <- expand.grid(lambda = paste0("", lambda_values), mu = paste0("", mu_values))
  LLdf$prob <- unlist(prob)
  LLdf$lambda <- as.numeric(paste(LLdf$lambda))
  LLdf$mu <- as.numeric(paste(LLdf$mu))
  LLdf$meanN <- 0
  LLdf$meanNb <- 0

  for (j in lambda_values) {
    # calculate mean N
    dfPn <- pmf.ini_size(j, R0, maxIni)
    dfPn <- dfPn[dfPn$type == "adjusted",]
    dfPn$prob_N <- dfPn$N * dfPn$prob
    mean_n <- sum(dfPn$prob_N)
    LLdf[which(near(LLdf$lambda, j)),]$meanN <- mean_n

    # calculate mean Nb
    mean_nb <- 0
    for (k in 1:maxNb) {
      PNb <- 0
      for (n in 1:maxIni) {
        PNb <- PNb + dfPn$prob[n+1] * dbinom(k, n, 1-1/R0)
      }
      mean_nb <- mean_nb + k * PNb
    }
    LLdf[which(near(LLdf$lambda, j)),]$meanNb <- mean_nb
  }

  return(LLdf)
}


#' Expected Clonal Mutations
#'
#' This function calculates the expected number of clonal mutations.
#'
#' @param n Initial population size, a positive integer
#' @param R0 Reproduction number, a positive number
#' @param mu Mutation rate, a positive number
#' @param maxMuGen Maximum number of mutant lineages being calculated
#' @param maxFS Maximum number of final sizes being calculated
#' @param maxClonal Maximum number of clonal mutations being calculated
#' @return A number representing the expected clonal mutations
#' @examples
#' meanClonal(2, 1.2, 0.4, 60, 60, 5)
#' meanClonal(3, 3.5, 1.6, 120, 130, 5)
#' @export
meanClonal <- function(n, R0, mu, maxMuGen, maxFS, maxClonal) {
  df <- pmf.P_clonal(n, R0, mu, maxMuGen, maxFS, maxClonal)
  sum <- 0
  for (i in 1:(maxClonal + 1)) {
    sum <- sum + i * df[i + 1, 2]
  }
  return(sum)
}

#' Expected Clonal Mutations Across Parameters
#'
#' This function calculates the expected number of clonal mutations for a
#' sequence of R0 and mu values.
#'
#' @param n Initial population size, a positive integer
#' @param R0_values A sequence of reproduction number
#' @param mu_values A sequence of mutation rate
#' @param maxMuGen Maximum number of mutant lineages being calculated
#' @param maxFS Maximum number of final sizes being calculated
#' @param maxClonal Maximum number of clonal mutations being calculated
#' @return A number representing the expected clonal mutations
#' @examples
#' meanClonal.df(2, seq(1.1, 7.1, by = 0.1), seq(0.001, 0.101, by = 0.002), 60, 60, 5)
#' meanClonal.df(3, seq(1.05, 3, by = 0.05), seq(0.3, 0.5, by = 0.1), 120, 130, 5)
#' @export
meanClonal.df <- function(n, R0_values, mu_values, maxMuGen, maxFS, maxClonal) {
  mean <- 1:length(R0_values) %>%
    mclapply(function(x){
      1:length(mu_values) %>%
        lapply(function(y){
          meanClonal(n,R0_values[x], mu_values[y], maxMuGen, maxFS, maxClonal)
        })
    }, mc.cores = detectCores())

  mcdf <- expand.grid(mu = paste0("", mu_values), R0 = paste0("", R0_values))
  vs_slope = c()

  mcdf$mc <- unlist(mean)
  mcdf <- mcdf[, c("R0", "mu", "mc")]
  mcdf$logmc <- log(mcdf$mc)
  return(mcdf)
}
