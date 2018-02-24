//# cleaned up code based on Bob Carpenter suggestions
//# http://discourse.mc-stan.org/t/group-size-parameters-in-ridge-regression-with-mean-embedding-kernel/3355/2
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
  vector[N] alpha;
  ident_N = diag_matrix(rep_vector(1, rows(K)));
  m = K + lambda * ident_N;
  for (n in 1:N) {
    alpha[n] = 1.0 / N;
  }
  print(m);
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
  Kalpha = K * alpha;
  print("trans params");
  // spec = 1 + exp(-Kalpha);
  spec = inv(inv_logit(Kalpha));
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
  print(Kalpha);
}
model {
  print("model");
  // this block results approximate KRR_logit() analytical solution
  alpha_hat ~ normal(0,100); // 10
  q ~ normal(m * alpha_hat, 0.1); // 0.01
}
generated quantities{
  vector[N] yhat1 = inv_logit(K * alpha_hat);
  // vector[N] yhat1 = inv(yhat2);
}

