data {
  int<lower=0> N;
  int<lower=0> P;
  matrix[P, P] K;
  vector[N]    y;
  real         lambda;
}
transformed data{
  // this block results verified with KRR_logit() analytical solution
  matrix[N, N] m;
  matrix[N, N] ident_N;
  ident_N = diag_matrix(rep_vector(1, rows(K)));
  m = K + lambda * ident_N;
  //print(m)
}
parameters {
  vector[N] alpha_hat;
}
transformed parameters{
  // this block results verified with KRR_logit() analytical solution
  vector[N] q;
  vector[N] e_;
  vector[N] diagW;
  vector[N] pi_;
  vector[N] spec;
  vector[N] Kalpha;
  vector[N] alpha;
  for (n in 1:N) {
    alpha[n] = 1.0 / N;
  }
  Kalpha = K * alpha;
  spec = 1 + exp(-Kalpha);
  for (n in 1:N){
    pi_[n] = 1.0 / spec[n];
  }
  for (n in 1:N){
    diagW[n] = pi_[n] * (1.0 - pi_[n]);
  }
  for (n in 1:N){
    e_[n] = (y[n] - pi_[n]) / diagW[n];
  }
  for (n in 1:N){
    q[n] = Kalpha[n] + e_[n];
  }
}
model {
  // this block results approximate KRR_logit() analytical solution
  alpha_hat ~ normal(0,10); // for starters
  q ~ normal(m * alpha_hat, 0.01);
}
generated quantities{
  vector[N] yhat1;
  vector[N] yhat2;
  yhat1 =  (1.0 + exp(-(K * alpha_hat)));
  for (n in 1:N){
    yhat2[n] = 1.0 / yhat1[n];
  }
}


// KRR_logit <- function(K,y,lambda){
//   #### Logistic KRR
//    N = nrow(K)
//   alpha = rep(1/N, N) # initial values of alpha, transformed Parameters block
//   Kalpha = as.vector(K %*% alpha) # as 1D matrix of vector? stan wants a vector it seems
//   spec = 1 + exp(-Kalpha) # transformed Parameters block
//   pi = 1 / spec # transformed Parameters block
//   diagW = pi * (1 - pi) # transformed Parameters block
//   e = (y - pi) / diagW # transformed Parameters block // errors started here
//   q = Kalpha + e  # transformed Parameters block // errors continued here
//   ident.N <- diag(rep(1,N)) # constructed in model block
//   theSol = solve(K + lambda * ident.N, q) # the objective
//   log_pred <- 1 / (1 + exp(-as.vector(K %*% theSol))) # generated quantities
//   return(list(pred = log_pred, alphas = theSol))
// }

