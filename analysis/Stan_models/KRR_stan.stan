data {
  int<lower=0> N;
  matrix[N, N] K;
  vector[N] y;
  real lambda;
}
transformed data{
  matrix[N, N] m;
  matrix[N, N] ident_N;
  ident_N = diag_matrix(rep_vector(1, rows(K)));
  m = K + lambda * ident_N;
  //print(m);
}
parameters {
  vector[N] alpha;
}
model {
  alpha ~ normal(0,1); // for starters
  y ~ normal(m * alpha, 0.01);
}
generated quantities{
  vector[N] yhat;
  yhat = K * alpha;
}


