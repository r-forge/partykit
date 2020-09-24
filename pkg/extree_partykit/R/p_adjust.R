p_adjust_bonferroni <- function(p, n = sum(!is.na(p)), log = FALSE) {
  ## number of tests to adjust for
  if(!is.integer(n)) n <- as.integer(round(n))

  ## no adjustment needed for a single test
  if(n == 1L) return(p)

  if(!log) p * n else p + log(n)
}

p_adjust_sidak <- function(p, n = sum(!is.na(p)), log = FALSE) {
  ## number of tests to adjust for
  if(!is.integer(n)) n <- as.integer(round(n))

  ## no adjustment needed for a single test
  if(n == 1L) return(p)

  ## Sidak adjustment
  padj <- if(!log) {
    1 - (1 - p)^n
  } else {
    log(1 - (1 - exp(p))^n)
  }
  
  ## Bonferroni adjustment for very small p-values to avoid numerical instabilities,
  ## handling discontinuity close to threshold
  thresh <- if(!log) 1e-6 else log(1e-6)
  if(any(small <- p < thresh, na.rm = TRUE)) {
    padj[which(small)] <- if(!log) { ## Bonferroni
      p[which(small)] * n
    } else {
      p[which(small)] + log(n)
    }
    plin <- if(!log) { ## transition
      (1 - (1 - thresh)^n) + (p - thresh)
    } else {
      log((1 - (1 - exp(thresh))^n) + (exp(p) - exp(thresh)))
    }
    close <- small & (padj > plin)
    padj[which(close)] <- plin[which(close)]
  }
  return(padj)
}
