using Optim
using StatsFuns
using DataFrames

# inverse beta distribution function を行列に対応するように拡張
import StatsFuns.betainvcdf
betainvcdf(alpha::Number, beta::Number, x::Array) = reshape([betainvcdf(alpha, beta, i) for i in reshape(x, 1, size(x, 1)*size(x, 2)) ], size(x, 1), size(x, 2))

# maxをArray{String, 1}に対応するように拡張
# しかしArray型で入っているのでAnyに対応させる
import Base.max
max(number::Real, comparison::Any) = [max(number, parse(Int, i)) for i in comparison] 

# normpdfを配列に拡張
import StatsFuns.normpdf
normpdf(array::Array{Float64, 2}) = reshape([normpdf(i) for i in reshape(array, 1, size(array, 1)*size(array, 2))], size(array, 1), size(array, 2))