function parallel(nproc, seed)
   
    addprocs(nproc)
    @everywhere include("preparation.jl")
    @everywhere include("bayes.jl")
    @everywhere include("loglike.jl")
    
    srand(seed)
    initialvalues = randn(273, nproc)
    
    writetable(*("initialvalue", "_$seed", ".txt"), DataFrame(initialvalues))
    
    estparam = Array(Float64, 273, nproc)
    likelihood = Array(Float64, 1, nproc)
    result = Array(Any, nproc, 1)
    
    for i in 1:nproc
        result[i, 1] = remotecall(optimize, i+1, loglike, initialvalues[:, i], BFGS())
    end
    
    for j in 1:nproc
        wait(result[j, 1])
        z = fetch(result[j, 1])
        estparam[:, j] = Optim.minimizer(z)
        likelihood[1, j] = Optim.minimum(z)
    end
    
    writetable(*("estparam", "_$seed", ".txt"), DataFrame(estparam))
    writetable(*("likelihood", "_$seed", ".txt"), DataFrame(likelihood))

end

parallel(48, 1)