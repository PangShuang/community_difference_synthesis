################################################################################
##  magnitude_regressions_noninf.py: Bayesian analysis with non-informative priors examining whether the magnitude of N, irrigation, or drought treatments affect the richness or compositonal differences.
##
##  Author: Nate Lemoine
##  Date created: January 16, 2018
################################################################################


import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
import pandas as pd
import seaborn as sns
sns.set(context='talk', style='ticks', font_scale=1.2, rc={'xtick.direction':'in', 'ytick.direction':'in', 'xtick.top': True, 'ytick.right': True})
import pystan as pyst
from patsy import dmatrix
from tabulate import tabulate


#%% LOAD DATA
irr_data = pd.read_csv('/nfs/corre-data/ForAnalysis_allAnalysisH2Omag_irr.csv')
dr_data = pd.read_csv('/nfs/corre-data/ForAnalysis_allAnalysisH2Omag_drought.csv')
N_data = pd.read_csv('/nfs/corre-data/ForAnalysis_allAnalysisNmag.csv')


#%% SET UP THE MODEL
stan_mod="""
data{
	int<lower=1> N;
	int<lower=1> K;
	vector[N] Y;
	matrix[N,K] X;
}
parameters{
	real<lower=0> sd_y;
	vector[K] B;
}
model{
	Y ~ normal(X*B, sd_y);
}
"""
reg_mod =pyst.StanModel(model_code=stan_mod)


#%% ANALYSIS FUNCTIONS
def std_func(x):
	mu = np.mean(x); std = np.std(x, ddof=1)
	return (x - mu)/std

def kde_plot(data):
	pd.DataFrame(data).plot(kind='kde')
	plt.show()

def analysis_func(dataset, response, predictor1, predictor2, filename):
	y = std_func(dataset[response])
	formula = '~std_func({0})*std_func({1})'.format(predictor1, predictor2)
	X = dmatrix(formula, data=dataset)
	mod_data = {'N': y.shape[0], 'K': X.shape[1], 'Y': y, 'X': np.asarray(X)}
	fit = reg_mod.sampling(data=mod_data, iter=20000, chains=4)
	b = fit.extract('B')['B']
	columns = X.design_info.column_names
	df = pd.DataFrame(b, columns=columns)
	df.to_csv('/nfs/corre-data/posteriors_{0}.csv'.format(filename), index=False)
	return df


#%% POST-PROCESSING FUNCTIONS
def results_table(data, ci):
	names = data.columns
	means = list(data.mean())
	low_ci = np.round((1. - ci)/2.0,3)
	high_ci = (1. - low_ci)
	low = list(data.quantile(low_ci))
	high = list(data.quantile(high_ci))
	pr = np.mean(data > 0)
	pr = list(np.maximum(pr, 1-pr))
	index = ['Mean', 'CI {0}'.format(low_ci), 'CI {0}'.format(high_ci), 'Pr(Effect)']
	tab = tabulate([means, low, high, pr], headers=names, showindex=index)
	print(tab)


#%% RUN THE ANALYSIS
n_mc_posts = analysis_func(N_data, 'mean_change', 'n', 'MAP', 'N_MeanChange')
n_rich_posts = analysis_func(N_data, 'S_PC', 'n', 'MAP', 'N_Richness')
irr_mc_posts = analysis_func(irr_data, 'mean_change', 'precip', 'MAP', 'Irr_MeanChange')
irr_rich_posts = analysis_func(irr_data, 'S_PC', 'precip', 'MAP', 'Irr_Richness')
dr_mc_posts = analysis_func(dr_data, 'mean_change', 'precip', 'MAP', 'Drought_MeanChange')
dr_rich_posts = analysis_func(dr_data, 'S_PC', 'precip', 'MAP', 'Drought_Richness')


#%% PROCESS RESULTS
results_table(n_mc_posts, 0.95)
results_table(n_rich_posts, 0.95)
results_table(irr_mc_posts, 0.95)
results_table(irr_rich_posts, 0.95)
results_table(dr_mc_posts, 0.95)
results_table(dr_rich_posts, 0.95)

# PLOT FUNCTIONS
def int_plot(data, posts, response, predx, pred_cat, cat_label, x_label, y_label, filepath):
	response2 = std_func(data[response])
	pred_cat2 = std_func(data[pred_cat])
	predx2 = std_func(data[predx])
	quantiles = [pred_cat2.min(), np.percentile(pred_cat2, 50), pred_cat2.max()]
	pal = sns.color_palette()
	xt = np.linspace(predx2.min(), predx2.max(), 100)
	labels = [x + ' {0}'.format(cat_label) for x in ['Low', 'Medium', 'High']]
	plt.plot(np.array(predx2), np.array(response2), 'ko', ms=8)
	for i in range(3):
		yhat = posts.apply(lambda x: pd.Series(x[0] + x[1]*quantiles[i] + x[2]*xt + x[3]*quantiles[i]*xt), axis=1) 
		plt.plot(xt, yhat.mean(axis=0), c=pal[i], label=labels[i])
		plt.fill_between(xt, np.percentile(yhat,2.5,axis=0), np.percentile(yhat,97.5,axis=0), 
			facecolor=pal[i], edgecolor='None', alpha=0.5)

	plt.legend(bbox_to_anchor=(1.01, 1), loc=2)
	plt.xlabel('Standardized {0}'.format(x_label))
	plt.ylabel('Standardized {0}'.format(y_label))
	plt.savefig('/nfs/corre-data/{0}.pdf'.format(filepath))
	plt.show()

def partial_plot(data, posts, response, pred, pred_label, ylabel, filename):
	formula = '~std_func({0})*std_func(MAP)'.format(pred)
	X = np.asarray(dmatrix(formula, data=data))
	response2 = std_func(data[response])
	yhat = X.dot(posts.T)
	err = response2[:,np.newaxis] - yhat
	# map
	yhat_err1 = X[:,[0,2]].dot(posts.iloc[:,[0,2]].T)
	err_map = err + yhat_err1
	map_std = std_func(data['MAP'])
	x_map = np.linspace(map_std.min(), map_std.max(), 100)
	yhat_map = posts.apply(lambda x: pd.Series(x[0]+x[2]*x_map),axis=1)
	plt.plot(map_std, err_map.mean(axis=1), 'ko', ms=8)
	plt.fill_between(x_map, np.percentile(yhat_map,2.5,axis=0), np.percentile(yhat_map, 97.5,axis=0), 
		facecolor='grey', edgecolor='None')
	plt.plot(x_map, yhat_map.mean(), 'r-', lw=3)
	plt.xlabel('Standardized MAP')
	plt.ylabel('{0}\nPartial Residuals'.format(ylabel))
	plt.savefig('/nfs/corre-data/{0}-partialMAP.pdf'.format(filename))
	plt.show()
	# other predictor
	yhat_err2 = X[:,[0,1]].dot(posts.iloc[:,[0,1]].T)
	err_o = err + yhat_err2
	o_std = std_func(data[pred])
	x_o = np.linspace(o_std.min(), o_std.max(), 100)
	yhat_o = posts.apply(lambda x: pd.Series(x[0] + x[1]*x_o), axis=1)
	plt.plot(o_std, err_o.mean(axis=1), 'ko', ms=8)
	plt.fill_between(x_o, np.percentile(yhat_o,2.5,axis=0), np.percentile(yhat_o,97.5,axis=0), 
		facecolor='grey', edgecolor='None')
	plt.plot(x_o, yhat_o.mean(), 'r-', lw=3)
	plt.xlabel('Standardized {0}'.format(pred_label))
	plt.ylabel('{0}\nPartial Residuals'.format(ylabel))
	plt.savefig('/nfs/corre-data/{0}-partial{1}.pdf'.format(filename,pred_label))
	plt.show()


# PLOTS
int_plot(N_data, n_mc_posts, 'mean_change', 'MAP', 'n', 'N', 'MAP', 'Mean Change', 'MeanChange_N_MAP')
partial_plot(N_data, n_rich_posts, 'S_PC', 'n', 'N', 'Richness', 'Richness_N_MAP')
partial_plot(irr_data, irr_rich_posts, 'S_PC', 'precip', 'Water Amount', 'Richness', 'Richness_Irr_MAP')