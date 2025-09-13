
import pymc as pm
import pandas as pd
import numpy as np
import pickle
import scipy
from scipy import stats



class EeToSigma:
    GAS_CONST = 1.9872036
    def __init__(self, Rs, Ss, Ts, exp_ee, exp_sigma):
        self.Rs = Rs
        self.Ss = Ss
        N = Rs.shape[1]
        self.Ts = Ts
        self.T_mat = np.repeat(Ts, N).reshape(-1, N)
        self.EEs = exp_ee
        self.exp_sigma = exp_sigma
        self.model = pm.Model()

        with self.model:
            # prior
            calc_sigma = pm.Uniform("calc_sigma", 0,10)
            # model
            rs = pm.Normal('rs', mu = self.Rs, sigma = calc_sigma)
            ss = pm.Normal('ss', mu = self.Ss, sigma = calc_sigma)
            rexp = pm.math.sum(pm.math.exp( -1000  * rs / EeToSigma.GAS_CONST / self.T_mat  ), axis = 1)
            sexp = pm.math.sum(pm.math.exp( -1000  * ss / EeToSigma.GAS_CONST / self.T_mat  ), axis = 1)
            ee = 2 * 100 * rexp / (rexp + sexp)  - 100
            # observables
            ee_obs = pm.Normal("ee_obs", mu = ee, sigma = self.exp_sigma, observed = self.EEs)

    def fit(self, draws = 25000, var_names = ["calc_sigma"], **kwargs):
        with self.model:
            self.idata = pm.sample(draws = draws,
                                   var_names = var_names,chains=4,cores=4, **kwargs)
            self.sigma_distr = np.array(self.idata['posterior']['calc_sigma']).flatten()

# read data created in "data_preparation" as "wheeler.pcl" 
with open("wheeler.pcl", 'rb') as f:
    wheeler = pickle.load(f)

T = np.array(wheeler['Ts'], dtype=float)
EEs = np.array(wheeler['EEs'], dtype=float)
Exp_error = 2
catalysts = ['1a', '1b', '2', '3a', '3b', '4a', '4b', '4c', '4d', '4e', '5a', '5b', '5c', '6']
method = 'PCM-B97D_E'
Rs_meas = wheeler['Rs'][method]
Ss_meas = wheeler['Ss'][method]
results = []
#get distributions for different catalysts
for i, catalyst in enumerate(catalysts):
    Rs_cat = Rs_meas[i]
    Ss_cat = Ss_meas[i]
        
    model = EeToSigma(Rs_cat.reshape(1, -1), Ss_cat.reshape(1, -1), [T[i]], [EEs[i]], Exp_error)
    model.fit(draws=25000,  njobs=4, chains=4, target_accept=0.99)
    df = pd.DataFrame({'Distr': model.sigma_distr, 'Method': method, 'Cat':catalyst})
    results.append(df)
dfs = pd.concat(results)
dfs.to_csv("PCM-B97D_E_Cats.csv")






