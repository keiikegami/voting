import numpy as np
from scipy.special import factorial

class ExpectedChi:
    
    def __init__(self, lam = 1.0, T = 5, sim = 10, random_state = 1):
        self.lam = lam # parameter of exponential distribution which rho's follow
        self.T = T # time
        self.sim = sim # simulation number
        self.random_state = random_state
    
    def setparam(self):
        # params[0, :]　：　rho_eta
        # params[(1,2,3,4), :] : rho_chi for each candidate
        # params[(5,6,7,8), :] : mu_chi for each candidate
        # params[(9,10,11,12), :] : chi for each candidate
        r = np.random.RandomState(self.random_state)
        rho = r.exponential(self.lam, (5, self.sim))
        mu_chi = r.normal(0.0,1.0,(8, self.sim))
        params = np.vstack((rho, mu_chi))
        return params
        
    def signals(self):
        # make the 4 dimentional array : all_signals
        # each element in all_signals is 3 dimentional array whose shape is (T! * 4 * T )
        params = self.setparam()
        all_signals = np.empty([self.sim, min(100, factorial(self.T).astype(int)), 4, self.T])
        g = np.random.RandomState(self.random_state)
        for i in range(self.sim):
            for k in range(min(100, factorial(self.T).astype(int))):
                for j in range(4):
                    all_signals[i, k, j, :] = g.normal(params[j + 9, i], 1/params[0, i],  (1, self.T))
        return all_signals
    
    def cum_signals(self):
        # calculating cumulative signals for expected chi_k given omega
        sig = self.signals()
        for i in range(self.sim):
            for j in range(self.T):
                if j != 0:
                    sig[i, :, :, :][:, :, j] = sig[i, :, :, :][:, :, j] + sig[i, :, :, :][:, :, j -1]
        return sig
            
    def expected_chi(self):
        # calculating expected chi_k given omega
        cum = self.cum_signals()
        params = self.setparam()
        transform = np.empty((self.sim, min(100, factorial(self.T).astype(int)), 4, self.T))
        
        # first multiply cum by rho_eta 
        for i in range(self.sim):
            transform[i, :, :, :] = cum[i, :, :, :] * params[0, i]
            
        # second add rho_chi + mu_chi to transform
        for i in range(self.sim):
            for j in range(4):
                transform[i, :, j, :] = transform[i, :, j, :] + (params[j + 1, i] * params[j + 5, i])
        
        # third divide by rho_chi + t * rho_eta
        for i in range(self.sim):
            for j in range(4):
                for t in range(self.T):
                    transform[i, :, j, t] = transform[i, :, j, t] / (params[j + 1, i] + t * params[0, i])
        
        return transform
    
    def means(self, each_sim = True):
        # calculate mean
        if each_sim == True:
            # calculate the means of expected chi for each parameter set and t
            # return matrix
            return self.expected_chi().mean(axis = (1,2))
        
        else:
            # calcualte the over all mean
            return self.expected_chi().mean(axis = (0,1,2))
        
    def variances(self, each_sim = True):
        # calculate variance
        if each_sim == True:
            # calculate the variances of expected chi for each parameter set and t
            # return list
            return self.expected_chi().var(axis = (1,2))
        
        else:
            # calcualte the over all variances
            return self.expected_chi().var(axis = (0,1,2))