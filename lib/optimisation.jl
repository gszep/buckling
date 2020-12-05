using ForwardDiff: gradient!

function loss( x::Tuple{<:AbstractVector,<:AbstractVector}; kwargs... )
    data,model = x
    boundary!(model,data;kwargs...)
    return log(norm( data - model ))
end

function loss( data::Vector{<:Tuple{<:AbstractVector,<:AbstractVector}}; kwargs... )
    return sum( x->loss(x;kwargs...), data ) / length(data)
end

function gradients!( gradients::Dict, data::Vector{<:Tuple{<:AbstractVector,<:AbstractVector}};
         weights=[0.0], origin=[0.0,0.0], orientation=[75.0,0.0], kwargs... )

    gradient!( gradients[:weights], θ->loss(data;weights=θ,origin=origin,orientation=orientation,kwargs...), weights )
    gradient!( gradients[:origin],  θ->loss(data;weights=weights,origin=θ,orientation=orientation,kwargs...), origin  )

    return gradients
end

import Flux
function Flux.update!(opt, x::Dict, x̄::Dict)
    for key ∈ keys(x)
        Flux.update!(opt,x[key],x̄[key])
    end
end