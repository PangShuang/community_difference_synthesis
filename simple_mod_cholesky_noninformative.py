################################################################################
##  simple_mod_cholesky_noninformative.py: Bayesian analysis with non-informative priors.
##
##  Author: Nate Lemoine
##  Date created: September 5, 2018
################################################################################

import pickle
from pystan import StanModel


simple_model = """
data{
	int<lower=1> N;
	matrix[N,3] Time;
	matrix[N,2] Y;
	int<lower=1> experiment[N];
	int<lower=1> J;
	int<lower=1> K;
	matrix[J,K] trt_x;
	int<lower=1> community[J];
	int<lower=1> Z ;
	matrix[Z,3] anpp_rich;
}
parameters{
	matrix[K,3] U[2];

	matrix[3,J] effs_B[2];
	cholesky_factor_corr[3] L_Omega_B[2];
	vector<lower=0, upper=pi()/2>[3] tau_b_unif[2];

	matrix[3,Z] effs_G[2];
	cholesky_factor_corr[3] L_Omega_G[2];
	vector<lower=0, upper=pi()/2.>[3] tau_g_unif[2];

	matrix[3,3] D[2];

	cholesky_factor_corr[2] L_Omega_Y;
	vector<lower=0, upper=pi()/2>[2] tau_y_unif;
}
transformed parameters{
	matrix[J,3] B[2];
	matrix[J,3] Bhat[2];
	vector<lower=0>[3] tau_b[2];
	matrix[Z,3] G[2];
	matrix[Z,3] Ghat[2];
	vector<lower=0>[3] tau_g[2];


	for(m in 1:2){
		for(t in 1:3){
			tau_b[m,t] = 2.5 * tan(tau_b_unif[m,t]);
			tau_g[m,t] = 2.5 * tan(tau_g_unif[m,t]);
		}
		

		Ghat[m] = anpp_rich*D[m];
		G[m] = Ghat[m] + (diag_pre_multiply(tau_g[m], L_Omega_G[m])*effs_G[m])';

		Bhat[m] = G[m, community] + trt_x*U[m];
		B[m] = Bhat[m] + (diag_pre_multiply(tau_b[m], L_Omega_B[m])*effs_B[m])';
	}

	
}
model{
	matrix[N,2] Yhat;
	vector[2] tau_y;
	matrix[2,2] Sigma_Y;

	for(n in 1:N){
		for(m in 1:2){
			Yhat[n,m] = Time[n] * B[m,experiment[n]]';
		}
	}


	for(m in 1:2){
		tau_y[m] = 2.5 * tan(tau_y_unif[m]);
	}
	Sigma_Y = diag_pre_multiply(tau_y, L_Omega_Y);
	for(n in 1:N){
		Y[n] ~ multi_normal_cholesky(Yhat[n], Sigma_Y);
	}

	for(m in 1:2){
		L_Omega_B[m] ~ lkj_corr_cholesky(25.0);
		to_vector(effs_B[m]) ~ normal(0,1);
		L_Omega_G[m] ~ lkj_corr_cholesky(25.0);
		to_vector(effs_G[m]) ~ normal(0,1);
	}

	L_Omega_Y ~ lkj_corr_cholesky(25.0);
}
"""
simple_comp = StanModel(model_code=simple_model)

with open('simple_cholesky_noninformative.pkl', 'wb') as f:
	pickle.dump(simple_comp, f)
