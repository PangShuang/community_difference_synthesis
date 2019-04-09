################################################################################
##  timestd_PC.py: Bayesian analysis examining differences in the absolute value of richness and composition for datasets cutoff at 20 years.
##
##  Author: Nate Lemoine
##  Date created: January 16, 2018
################################################################################


#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import patsy
import pickle
import os
os.chdir('/nfs/corre-data/')

# type in the model name between quotes". Choose one of:
# cholesky_mod_noninformative
# cholesky_mod_N01
# cholesky_mod_normal-gamma
model_name = "cholesky_mod_noninformative"
# input dataframe name here WITHOUT the .csv extension
data_name = "timestdbytrt_20yr_PC"


#%% DATA PREP

# read in data and remove null ANPP values
corre_data = pd.read_csv("{0}.csv".format(data_name))
corre_data = corre_data[~corre_data['anpp'].isnull()]
# each site x project x community type gets its own unique ID
corre_data.loc[:,'community_ID'] = corre_data['site_code'] + '::' + corre_data['project_name'] + '::' + corre_data['community_type']
# each treatment within the community gets its own ID
corre_data.loc[:,'treatment_ID'] = corre_data['community_ID'] + '::' + corre_data['treatment']
# convert calendar year to treatment year
scale_fun = lambda x: x - x.min()
corre_data.loc[:,'time'] = corre_data.groupby('treatment_ID')['calendar_year'].transform(scale_fun)
# standardize the responses
std_fun = lambda x: (x - x.mean()) / x.std()
corre_data.loc[:,'mean_change_std'] = std_fun(corre_data['composition_diff'])
corre_data.loc[:,'s_change_std'] = std_fun(corre_data['S_PC'])
# turn the communities and treatments into integer IDs
corre_data.loc[:,'treat_INT'] = pd.factorize(corre_data['treatment_ID'])[0] + 1
corre_data.loc[:,'comm_INT'] = pd.factorize(corre_data['community_ID'])[0] + 1
# standardize time within each group
corre_data.loc[:,'time_std'] = corre_data.groupby(['comm_INT', 'treat_INT'])['time'].transform(lambda x: (x - x.mean()) / x.std())
# make the matrix for time
time_X = patsy.dmatrix('~time_std + np.power(time_std,2)', corre_data)
Y = np.array(corre_data[['mean_change_std', 's_change_std']])
# get the unique treatments
treatment_data = corre_data.groupby('treatment_ID').apply(lambda x: x[x['time']==x['time'].max()])
# get ANPP and richness for each treatment
treatment_data.loc[:,'ANPP_STD'] = std_fun(treatment_data['anpp'])
treatment_data.loc[:,'RICH_STD'] = std_fun(treatment_data['rrich']) 
# predictors at the treatment level
trt_X = patsy.dmatrix("C(trt_type, Sum)", treatment_data)[:,1:]
# get the unique communities
comm_data = treatment_data.drop_duplicates('community_ID')
anpp_rich_X = patsy.dmatrix('~ANPP_STD + RICH_STD', comm_data)

#%% model
compMod = pickle.load(open("{0}.pkl".format(model_name), 'rb'))


#%%
multi_data = {'N': Y.shape[0],
               'J': corre_data['treat_INT'].max(),
               'K': trt_X.shape[1],
               'Z': treatment_data['comm_INT'].max(),
               'Y': Y,
               'Time': np.array(time_X),
               'W': np.array(trt_X),
               'U': np.array(anpp_rich_X),
               'treatment': corre_data['treat_INT'],
               'community': treatment_data['comm_INT']
                }
    
multi_fit = compMod.sampling(data=multi_data, iter=10000, chains=4,
                             sample_file="/nfs/corre-data/{0}-{1}".format(model_name, data_name))


#%%
D = multi_fit.extract('D')['D']
E = multi_fit.extract('E')['E']


responses = ['MC', 'SPC']
effects = ['int', 'slope', 'quad']
for i in range(2):
    for j in range(3):
        D_df = pd.DataFrame(D[:,i,:,j])
        D_df.plot(kind='kde')
        plt.savefig("{0}-{1}_kdeplots_{2}_ANPP_RICH-{3}.pdf".format(model_name, data_name, responses[i], effects[j]))
        plt.close()
        E_df = pd.DataFrame(E[:,i,:,j])
        E_df.plot(kind='kde')
        plt.savefig("{0}-{1}_kdeplots_{2}_TRT-{3}.pdf".format(model_name, data_name, responses[i], effects[j]))
        plt.close()

#%%
B = multi_fit.extract('B')['B']
B_MC = B[:,0,:,:]
B_SPC = B[:,1,:,:]



for i in range(1,B_MC.shape[1]):
    tempdf = corre_data.loc[corre_data['treat_INT']==i]
    tempBMC = B_MC[:,i-1,:]
    tempBMC_mean = tempBMC.mean(axis=0)
    tempx = np.linspace(tempdf['time_std'].min(), tempdf['time_std'].max(), 100)
    tempy = tempBMC_mean[0] + tempBMC_mean[1]*tempx + tempBMC_mean[2]*tempx**2
    plt.plot(tempdf['time_std'], tempdf['mean_change_std'], 'ko')
    plt.plot(tempx, tempy, 'k-')
    plt.title(tempdf['treatment_ID'].unique())
    plt.savefig('PC_mc_figs_timestd_noninf/{}.png'.format(i))
    plt.close()
    tempBSPC = B_SPC[:,i-1,:]
    tempBSPC_mean = tempBSPC.mean(axis=0)
    tempx = np.linspace(tempdf['time_std'].min(), tempdf['time_std'].max(), 100)
    tempy = tempBSPC_mean[0] + tempBSPC_mean[1]*tempx + tempBSPC_mean[2]*tempx**2
    plt.plot(tempdf['time_std'], tempdf['s_change_std'], 'ko')
    plt.plot(tempx, tempy, 'k-')
    plt.title(tempdf['treatment_ID'].unique())
    plt.savefig('PC_s_figs_timestd_noninf/{}.png'.format(i))
    plt.close()
