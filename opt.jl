# delete OPTIONS1
# replace x0 by inivalue
result = optimize(new_loglike, parameter, BFGS())
theta = Optim.minimizer(result)
likelihood = Optim.minimum(result)