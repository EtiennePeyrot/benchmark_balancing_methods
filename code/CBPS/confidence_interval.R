# CBPS just-identified propensity block.
#
# target = "ATE": psi_i(beta) = x_i * (A_i - e_i) / (e_i * (1 - e_i))
# target = "ATT": psi_i(beta) = c_att * x_i * (A_i - e_i) / (1 - e_i)
#
# By default, ATT uses the same overall scaling constant as the CBPS package,
# c_att = n / sum(A). This scaling does not change the root, but it does mirror
# the package convention more closely.
propensity_blocks_logit_CBPS_JI =
  function(A, U, eta, target = c("ATE", "ATT"), att_package_scaling = TRUE) {
    target = match.arg(target)
    
    n = length(A)
    V = cbind(1, U)
    eta_lin = as.vector(V %*% eta)
    e = plogis(eta_lin)
    r = exp(eta_lin) # = e / (1 - e)
    
    if (target == "ATE") {
      h = (2 * A - 1) * (1 + exp((1 - 2 * A) * eta_lin))
      psi_eta_i = V * h
      J_weight = exp( -(2 * A - 1) * eta_lin)
      J_eta_eta = -crossprod(V, V * J_weight) / n
    } else {
      c_att = if (att_package_scaling) n / sum(A) else 1
      h = c_att * ifelse(A == 1, 1, -r)
      psi_eta_i = V * h
      J_weight = ifelse(A == 1, 0, r)
      J_eta_eta = -c_att * crossprod(V, V * J_weight) / n
    }
    
    list(
      e = e,
      eta_lin = eta_lin,
      psi_eta_i = psi_eta_i,
      J_eta_eta = J_eta_eta,
      Xeta = V
    )
  }

outcome_blocks_binomial = function(Y, A, X) {
  n = length(Y)
  Z = cbind(1, as.matrix(X))
  p_z = ncol(Z)
  
  # Separate logistic outcome fits
  fit0 = glm(Y ~ Z - 1, subset = !A, family = binomial())
  fit1 = glm(Y ~ Z - 1, subset =  A, family = binomial())
  b0   = fit0$coefficients
  b1   = fit1$coefficients
  
  m0   = plogis(as.vector(Z %*% b0))
  m1   = plogis(as.vector(Z %*% b1))
  m0p  = m0 * (1 - m0)
  m1p  = m1 * (1 - m1)
  
  # Residual used in DR augmentation
  resid0 = Y - m0
  resid1 = Y - m1
  
  # Diagonal blocks
  J_b0_b0 = - crossprod(Z, Z * ((1 - A) * m0p)) / n
  J_b1_b1 = - crossprod(Z, Z * (A * m1p)) / n
  
  # Reusable pieces for M and J
  m0pZ = Z * m0p
  m1pZ = Z * m1p
  psi_b0_i = Z * (1 - A) * resid0
  psi_b1_i = Z * A * resid1
  
  list(
    Z = Z, p_z = p_z,
    m0 = m0, m1 = m1, m0p = m0p, m0pZ = m0pZ,
    b0 = b0, b1 = b1, m1p = m1p, m1pZ = m1pZ,
    resid0 = resid0, resid1 = resid1,
    psi_b0_i = psi_b0_i, psi_b1_i = psi_b1_i,
    J_b0_b0 = J_b0_b0, J_b1_b1 = J_b1_b1,
    pi_hat = mean(A)
  )
}

wls_ci = function(Y, A, ps_block, target=c("ATE","ATT"), alpha=0.05) {
  target = match.arg(target)
  n = length(Y)
  B = cbind(1, A)
  zero = function(r=1,c=1) matrix(0,nrow = r, ncol = c)
  
  ## propensity block
  e         = ps_block$e
  eta_lin   = ps_block$eta_lin
  psi_eta_i = ps_block$psi_eta_i
  J_eta_eta = ps_block$J_eta_eta
  Xeta      = ps_block$Xeta
  p_eta     = ncol(Xeta)
  
  ## weights and dw/deta
  r = exp(eta_lin) # e / (1 - e)
  if (target == "ATE") {
    w = 1 / ifelse(A == 1, e, 1 - e)
    dw_detai = (1 - 2 * A) * exp((1 - 2 * A) * eta_lin)
  } else {
    w = ifelse(A == 1, 1, r)
    dw_detai = ifelse(A == 1, 0, r)
  }
  
  ## WLS estimate
  BtWB = crossprod(B, B * w)
  BtWy = crossprod(B, Y * w)
  tau = solve(BtWB, BtWy)
  resid = as.numeric(Y - B %*% tau)
  
  # Build bread matrix J
  J_tau_tau = -BtWB / n
  J_tau_eta = crossprod(B, (resid * dw_detai) * Xeta) / n
  
  J = rbind(cbind(J_tau_tau    , J_tau_eta),
            cbind(zero(p_eta,2), J_eta_eta) )
  
  
  # build meat matrix M
  psi_tau_i = B * (w * resid)
  psi_i = cbind(psi_tau_i, psi_eta_i)
  M = crossprod(psi_i) / n
  
  
  # Sandwich variance
  Jinv = try(solve(J), silent = T)
  se_tau = if(inherits(Jinv, "try-error")) NA else {
    V = (1 / n) * Jinv %*% M %*% t(Jinv)
    se_tau = sqrt(V[2, 2])
  }
  tau_hat = tau[2]
  z = qnorm(1 - alpha / 2)
  ci = c(lower = tau_hat - z * se_tau,
         upper = tau_hat + z * se_tau)
  
  list(tau = tau_hat, ci = ci, se = se_tau)
}



dr_ci = function(A, outcome_block, ps_block, target = c("ATE","ATT"), alpha = 0.05) {
  target = match.arg(target)
  n = length(A)
  zero = function(r=1,c=1) matrix(0,nrow = r, ncol = c)
  
  # unpack
  Z         = outcome_block$Z
  p_z       = outcome_block$p_z
  m0        = outcome_block$m0
  m1        = outcome_block$m1
  m0pZ      = outcome_block$m0pZ
  m1pZ      = outcome_block$m1pZ
  resid0    = outcome_block$resid0
  resid1    = outcome_block$resid1
  J_b0_b0   = outcome_block$J_b0_b0
  J_b1_b1   = outcome_block$J_b1_b1
  psi_b0_i  = outcome_block$psi_b0_i
  psi_b1_i  = outcome_block$psi_b1_i
  pi_hat    = outcome_block$pi_hat
  
  e         = ps_block$e
  eta_lin   = ps_block$eta_lin
  psi_eta_i = ps_block$psi_eta_i
  J_eta_eta = ps_block$J_eta_eta
  Xeta      = ps_block$Xeta
  p_eta     = ncol(Xeta)
  
  # Build bread matrix J
  r = exp(eta_lin) # e / (1 - e)
  
  if (target == "ATE") {
    J_tau_tau = matrix(-1, 1, 1)
    J_tau_b0 = colMeans( ifelse(A == 1, -1, r) * m0pZ )
    J_tau_b1 = colMeans( ifelse(A == 1, -exp(-eta_lin), 1) * m1pZ )
    coeff_eta = - ifelse(A == 1, resid1, resid0) * exp((1 - 2 * A) * eta_lin)
  } else { # ATT
    J_tau_tau = matrix(-pi_hat, 1, 1)
    J_tau_b0 = -colMeans( ifelse(A == 1, 1, -r) * m0pZ )
    J_tau_b1 = rep(0, p_z)
    coeff_eta = ifelse(A == 1, 0, - resid0 * r)
  }
  J_tau_b0  = matrix(J_tau_b0, nrow = 1)
  J_tau_b1  = matrix(J_tau_b1, nrow = 1)
  J_tau_eta = matrix(colMeans(Xeta * coeff_eta), nrow = 1)
  
  
  J = rbind(
    cbind(J_tau_tau    , J_tau_b1       , J_tau_b0       , J_tau_eta),
    cbind(zero(p_z,1)  , J_b1_b1        , zero(p_z,p_z)  , zero(p_z,p_eta)),
    cbind(zero(p_z,1)  , zero(p_z,p_z)  , J_b0_b0        , zero(p_z,p_eta)),
    cbind(zero(p_eta,1), zero(p_eta,p_z), zero(p_eta,p_z), J_eta_eta)
  )
  
  # Build meat matrix M
  if (target == "ATE") {
    psi_tau_i = (m1 - m0) + ifelse(A == 1, resid1, resid0) * (2 * A - 1) * (1 + exp((1 - 2 * A) * eta_lin))
    tau_hat = mean(psi_tau_i)
    psi_tau = psi_tau_i - tau_hat
  } else { # ATT
    psi_tau_i = ifelse(A == 1, 1, -r) * resid0
    tau_hat = sum(psi_tau_i) / sum(A)
    psi_tau = psi_tau_i - A * tau_hat
  }
  
  M = crossprod(cbind(psi_tau, psi_b1_i, psi_b0_i, psi_eta_i)) / n
  
  # Sandwich variance for tau
  Jinv = try(solve(J), silent=T)
  se_tau = if(inherits(Jinv, "try-error")) NA else {
    Vcov = (1 / n) * Jinv %*% M %*% t(Jinv)
    se_tau = sqrt(Vcov[1, 1])
    }
  z = qnorm(1 - alpha / 2)
  ci = c(lower = tau_hat - z * se_tau,
         upper = tau_hat + z * se_tau)
  list(tau = tau_hat, ci = ci, se = se_tau)
}



wls_ci_optimist = function(Y, A, W, alpha=0.05) {
  n = length(Y)
  B = cbind(1, A)
  
  
  ## WLS estimate
  BtWB = crossprod(B, B * W)
  BtWy = crossprod(B, Y * W)
  tau = solve(BtWB, BtWy)
  resid = as.numeric(Y - B %*% tau)
  
  
  # Build bread matrix J
  J = -BtWB / n
  
  
  # build meat matrix M
  psi_i = B * (W * resid)
  M = crossprod(psi_i) / n
  
  
  # Sandwich variance
  Jinv = solve(J)
  V = (1 / n) * Jinv %*% M %*% t(Jinv)
  tau_hat = tau[2]
  se_tau = sqrt(V[2, 2])
  z = qnorm(1 - alpha / 2)
  ci = c(lower = tau_hat - z * se_tau,
         upper = tau_hat + z * se_tau)
  
  list(tau = tau_hat, ci = ci, se = se_tau)
}



dr_ci_optimist = function(A, outcome_block, W, target = c("ATE","ATT"), alpha = 0.05) {
  target = match.arg(target)
  n = length(A)
  zero = function(r=1,c=1) matrix(0,nrow = r, ncol = c)
  if (target == "ATT") W = W * sum(A) # scale weight for ATT to match psi_tau formula
  
  # unpack
  Z         = outcome_block$Z
  p_z       = outcome_block$p_z
  m0        = outcome_block$m0
  m1        = outcome_block$m1
  m0pZ      = outcome_block$m0pZ
  m1pZ      = outcome_block$m1pZ
  resid0    = outcome_block$resid0
  resid1    = outcome_block$resid1
  J_b0_b0   = outcome_block$J_b0_b0
  J_b1_b1   = outcome_block$J_b1_b1
  psi_b0_i  = outcome_block$psi_b0_i
  psi_b1_i  = outcome_block$psi_b1_i
  pi_hat    = outcome_block$pi_hat
  
  # Build bread matrix J
  if (target == "ATE") {
    J_tau_tau = matrix(-1, 1, 1)
    J_tau_b0  = matrix(colMeans((-1 + (1 - A) * W) * m0pZ), nrow = 1)
    J_tau_b1  = matrix(colMeans((1 - A * W) * m1pZ), nrow = 1)
  } else { # ATT  (one-sided, Hájek-normalized)
    J_tau_tau = matrix(-pi_hat, 1, 1)
    J_tau_b0  = matrix(colMeans(((1 - A) * W - A) * m0pZ), nrow = 1)
    J_tau_b1  = matrix(0, nrow = 1, ncol = p_z)
  }
  
  J = rbind(
    cbind(J_tau_tau    , J_tau_b1       , J_tau_b0),
    cbind(zero(p_z,1)  , J_b1_b1        , zero(p_z,p_z)),
    cbind(zero(p_z,1)  , zero(p_z,p_z)  , J_b0_b0)
  )
  
  # Build meat matrix M
  if (target == "ATE") {
    psi_tau_i = (m1 - m0) + (A * W * resid1 - (1 - A) * W * resid0)
    tau_hat   = mean(psi_tau_i)
    psi_tau   = psi_tau_i - tau_hat
  } else { # ATT
    psi_raw = (A - (1 - A) * W) * resid0
    tau_hat = mean(psi_raw) / pi_hat
    psi_tau = psi_raw - A * tau_hat
  }
  
  M = crossprod(cbind(psi_tau, psi_b1_i, psi_b0_i)) / n
  
  # Sandwich variance
  Jinv = solve(J)
  Vcov = (1 / n) * Jinv %*% M %*% t(Jinv)
  se_tau = sqrt(Vcov[1, 1])
  
  z  = qnorm(1 - alpha / 2)
  ci = c(lower = tau_hat - z * se_tau,
         upper = tau_hat + z * se_tau)
  
  list(tau = tau_hat, ci = ci, se = se_tau)
}
