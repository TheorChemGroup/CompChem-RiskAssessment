
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

    def fit(self, draws = 31250, var_names = ["calc_sigma"], **kwargs):
        with self.model:
            self.idata = pm.sample(draws = draws,
                                   var_names = var_names,chains=32,cores=32, **kwargs)
            self.sigma_distr = np.array(self.idata['posterior']['calc_sigma']).flatten()


# read data created in "data_preparation" as "wheeler.pcl" or as "wheeler_random.pcl" for data without random lowest TS
with open("wheeler.pcl", 'rb') as f:
    wheeler = pickle.load(f)

T = np.array(wheeler['Ts'], dtype=float)
EEs = np.array(wheeler['EEs'], dtype=float)
Exp_error = 2
selected_methods = ['PCM-B97D_E', 'PCM-B97D_H', 'PCM-B97D_G', 'PCM-B97D_qG', 'SMD-B97D','SMD-M06-2X', 'PCM-M06-2X', 'PCM-wB97XD', 'B97D']

#get distributions for different methods
results =[]
for method in selected_methods:
    Rs_meas = wheeler['Rs'][method]
    Ss_meas = wheeler['Ss'][method]
    model = EeToSigma(Rs_meas, Ss_meas, T, EEs, Exp_error)
    model.fit(draws=31250, target_accept=0.99)
    df = pd.DataFrame({'Distr': model.sigma_distr, 'Method': method})
    results.append(df)
    df.to_csv(f'Distr16.{method}.csv', index=False)



