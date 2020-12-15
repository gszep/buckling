function boundary!( model::AbstractVector, α::T; weights=[0.0], coordinates=[75.0,75.0,0.0,0.0,0.0], variance=0.01) where T<:Number
    X₀,θ = coordinates[3:4], coordinates[end]
    u,v = [cos(θ),sin(θ)], [-sin(θ),cos(θ)]
    a,b = coordinates[1:2]

    N = length(weights) # surface modes
    Rα = sum( n-> weights[n]*exp(-(α-(n/N-1/2)π)^2/variance), 1:N )

    model .= X₀ + (a+Rα)*cos(α-θ)*u + (b+Rα)*sin(α-θ)*v
    return model
end

function boundary!( model::AbstractVector, data::AbstractVector; coordinates=[75.0,75.0,0.0,0.0,0.0], kwargs...)
    x,y = coordinates[3:4]
    α = atan(data[2]-y,data[1]-x)
    return boundary!(model,α; kwargs...)
end

function boundary(parameters::Dict{Symbol,Vector{Vector{T}}}, weights::Union{Dict{Symbol,Vector{Vector{T}}},Nothing}=nothing; N=500) where T<:Number

    angles = range(-π,π,length=N)
    nImages = length(first(values(parameters)))
    boundaries = Dict([ channel=>[ zeros(N,2) for t ∈ 1:nImages ] for channel ∈ keys(parameters) ])

    for channel ∈ keys(parameters)
        for t ∈ 1:nImages
            for k ∈ 1:length(angles)
                boundaries[channel][t][k,:] = copy(boundary!( boundaries[channel][t][k,:], angles[k];
                    coordinates=parameters[channel][t], weights= isnothing(weights) ? [0.0] : weights[channel][t] ))
            end
        end
    end

    return boundaries
end