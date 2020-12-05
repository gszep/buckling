using ReverseDiff
include("./patches/StaticArrays.jl")

function loss( x::StaticVector{2,T}; kwargs... ) where T<:Number
    return log(norm( x - boundary(atan(x[2],x[1]);kwargs...) ) )
end

function loss( data::Vector{<:StaticVector}; kwargs... )
    return sum( x->loss(x;kwargs...), data ) / length(data)
end

function ∇loss( data::Vector{<:StaticVector}; weights=[0.0], origin=[0.0,0.0], orientation=[75.0,0.0], kwargs... )

    ∇weights = ReverseDiff.gradient(     θ->loss(data;weights=θ,origin=origin,orientation=orientation,kwargs... ), weights )
    ∇origin = ReverseDiff.gradient(      θ->loss(data;weights=weights,origin=θ,orientation=orientation,kwargs...), origin  )

    return Dict( :weights=>∇weights, :origin=>∇origin, :orientation=>[0.0,0.0] )
end

import Flux: update!
function update!(opt, x::Dict, x̄::Dict)
    for key ∈ keys(x)
        update!(opt,x[key],x̄[key])
    end
end