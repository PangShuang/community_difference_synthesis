################################################################################
##  simple_mod_work_effs_coding_20yr_absval.py: Bayesian analysis examining differences in the absolute value of richness and composition for datasets cutoff at 20 years.
##
##  Author: Nate Lemoine
##  Date created: January 16, 2018
################################################################################



import pandas as pd
import numpy as np
import pystan as pyst
import patsy
import pickle

corre_data = pd.read_csv('/nfs/corre-data/ForAnalysis_allAnalysis20yr_abs.csv')
corre_data = corre_data[~corre_data['anpp'].isnull()]
corre_data.loc[:,'community_code'] = corre_data['site_code'] + '::' + corre_data['project_name'] + '::' + corre_data['community_type']
corre_data.loc[:,'within_study_ID'] = corre_data['community_code'] + '::' + corre_data['treatment']

scale_fun = lambda x: x - x.min()
corre_data.loc[:,'time'] = corre_data.groupby('within_study_ID')['calendar_year'].transform(scale_fun)

std_fun = lambda x: (x - x.mean()) / x.std()
corre_data.loc[:,'mean_change_std'] = std_fun(corre_data['mean_change'])
corre_data.loc[:,'s_change_std'] = std_fun(corre_data['S_PC'])

corre_data.loc[:,'within_int'] = pd.factorize(corre_data['within_study_ID'])[0] + 1
corre_data.loc[:,'comm_INT'] = pd.factorize(corre_data['community_code'])[0] + 1

time_x = patsy.dmatrix('~time + np.power(time,2)', corre_data)
Y = np.array(corre_data[['mean_change_std', 's_change_std']])

study_data = corre_data.groupby('within_study_ID').apply(lambda x: x[x['time']==x['time'].max()])

study_data.loc[:,'ANPP_STD'] = std_fun(study_data['anpp'])
study_data.loc[:,'RICH_STD'] = std_fun(study_data['rrich']) 
trt_x = patsy.dmatrix("C(trt_type, Sum)", study_data)[:,1:]
trt_x_anpp = patsy.dmatrix("C(trt_type, Sum):ANPP_STD", study_data)[:,2:]
trt_x_rich = patsy.dmatrix("C(trt_type, Sum):RICH_STD", study_data)[:,2:]

trt_x2 = np.hstack((trt_x, trt_x_anpp, trt_x_rich))

comm_data = study_data.drop_duplicates('community_code')

anpp_rich_X = patsy.dmatrix('~ANPP_STD + RICH_STD', comm_data)

compMod = pickle.load(open('simple_cholesky.pkl', 'rb'))

simple_data_1 = {'N': Y.shape[0],
				'Time': np.array(time_x),
				'Y': Y,
				'experiment': np.array(corre_data['within_int']),
				'J': trt_x2.shape[0],
				'K': trt_x2.shape[1],
				'Z': study_data['comm_INT'].max(),
				'trt_x': np.array(trt_x2),
				'community': np.array(study_data['comm_INT']),
				'anpp_rich': np.array(anpp_rich_X),
				}
	
meanFit = compMod.sampling(data=simple_data_1, iter=10000, chains=4, 
  sample_file='/nfs/corre-data/simple_mod_effs_coding_interactions_20yrAbs')