# delete OPTIONS1
# replace x0 by inivalue

include("preparation.jl")
include("bayes.jl")
include("loglike.jl")

result = optimize(loglike, parameter, BFGS())

theta = Optim.minimizer(result)
likelihood = Optim.minimum(result)
