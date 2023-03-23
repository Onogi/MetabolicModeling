//Joint inference of flux based analysis and genomic prediction model
//Developed from MetabolicModeling10pi. Updates are
//Residual variances of Y and each Mv independently follow half_Cauchy distributions,

data {
  int J;//Number of reactions (J = Jc + Jn)
  int Jc;//Number of constrained reactions
  int Jn;//Number of reactions not constrained
  int K;//Number of metabolites
  int N;//Number of genotypes
  int Nn;//Number of non-missing values in Y
  matrix[N * K, J] M;//Stoichiometric matrix
  matrix[N, N] G;//Genomic relationship matrix
  vector[N] Y;//Biomass
  real W;//Weight
  vector[J] Muv;//Prior mean of V
  vector[J] Alpha;//Offset of flux
  matrix[N, K] Ma;//Response variable other than Y
  int Pos_Jc[Jc];//Positions of constrained reactions in 1:J
  int Pos_Jn[Jn];//Positions of non-constrained reactions in 1:J
  int Pos_NonMiss[Nn];//Positions of non-missing values in Y (1:N)
}

transformed data {
  matrix[N, N] L_G;
  real Marray[N, K, J];//Arranged M
  int Jk[N, K];//Number of non-zero elements in M
  int Pos[N, K, J];//Positions of non-zero elements in M
  real Cauchy_Mu = 0.0;
  real Cauchy_Sigma1_Y = 0.1;
  real Cauchy_Sigma2 = 10000.0;
  real Cauchy_Sigma3 = 10.0;
  real Cauchy_Sigma4 = 10.0;
  real Cauchy_Sigma1_M_shape = 0.1;
  real Cauchy_Sigma1_M_rate = 1.0;  
  
  for(i in 1:N)
    for(k in 1:K){
      Jk[i, k] = 0;
      for(j in 1:J){
        Marray[i, k, j] = M[(i - 1) * K + k, j];
        if(Marray[i, k, j] != 0.0){
          Jk[i, k] += 1;
          Pos[i, k, Jk[i, k]] = j;
        }
      }
    }

  L_G = cholesky_decompose(G);
}

parameters {
  vector<lower=0>[J] Sigma_U;//Genetic variance of flux
  vector<lower=0>[J] Sigma_R;//Residual variance of flux
  real<lower=0> Sigma_E;//Residual variance of biomass
  real<lower=0> Sigma_M;//Residual variance of Mv
  matrix[J, N] U_raw;//Raw samples of U
  matrix<lower=0>[Jc, N] Vc;//Flux (constrained)
  matrix[Jn, N] Vn;//Flux (not constrained)
  real A;//Intercept of biomass
  real<lower=0> B;//Slope of biomass
  real<lower=0> Cauchy_Sigma1_M;//Scaling parameter for the Cauchy prior for M
}

transformed parameters {
  matrix[J, N] U;//Genetic component of flux
  matrix[J, N] V;//Flux
  matrix[N, K] Mv;

  for(j in 1:J)
    U[j, ] = to_row_vector(rep_vector(Muv[j], N) + L_G * to_vector(U_raw[j, ]) * Sigma_U[j]); 

  for(j in 1:Jc)
    V[Pos_Jc[j], ] = Vc[j, ] + Alpha[Pos_Jc[j]];

  for(j in 1:Jn)
    V[Pos_Jn[j], ] = Vn[j, ] + Alpha[Pos_Jn[j]];

  for(i in 1:N)
    for(k in 1:K){
      Mv[i, k] = 0.0;
      for(j in 1:Jk[i, k])
        Mv[i, k] += Marray[i, k, Pos[i, k, j]] * V[Pos[i, k, j], i];
    }
}

model {
  Cauchy_Sigma1_M ~ gamma(Cauchy_Sigma1_M_shape, Cauchy_Sigma1_M_rate);
  Sigma_E ~ cauchy(Cauchy_Mu,Cauchy_Sigma1_Y);
  Sigma_M ~ cauchy(Cauchy_Mu,Cauchy_Sigma1_M);
  Sigma_U ~ cauchy(Cauchy_Mu,Cauchy_Sigma2);
  Sigma_R ~ cauchy(Cauchy_Mu,Cauchy_Sigma3);
  A ~ cauchy(Cauchy_Mu,Cauchy_Sigma4);
  B ~ cauchy(Cauchy_Mu,Cauchy_Sigma4);

  to_vector(U_raw) ~ normal(0, 1);
  for(j in 1:Jc)
    Vc[j, ] ~ normal(U[Pos_Jc[j], ], Sigma_R[Pos_Jc[j]]);
  for(j in 1:Jn)
    Vn[j, ] ~ normal(U[Pos_Jn[j], ], Sigma_R[Pos_Jn[j]]);

  target += normal_lpdf(Y[Pos_NonMiss] | A + B * to_vector(V[J, Pos_NonMiss]), Sigma_E);
  for(k in 1:K)
    target += W * normal_lpdf(to_vector(Ma[, k]) | to_vector(Mv[, k]), Sigma_M);
}

generated quantities{
  vector[J] Sigma2_U;
  vector[J] Sigma2_R;
  vector[J] h2;
  real Sigma2_E;
  real Sigma2_M;

  Sigma2_U = Sigma_U .* Sigma_U;
  Sigma2_R = Sigma_R .* Sigma_R;
  for(j in 1:J) h2[j] = Sigma2_U[j] / (Sigma2_U[j] + Sigma2_R[j]);

  Sigma2_E = Sigma_E * Sigma_E;
  Sigma2_M = Sigma_M * Sigma_M;
}
