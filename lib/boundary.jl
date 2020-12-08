function boundary!( model::AbstractVector, α::T; weights=[0.0], coordinates=[75.0,75.0,0.0,0.0,0.0], variance=0.01) where T<:Number
    X₀,θ = coordinates[3:4], coordinates[end]
    u,v = [cos(θ),sin(θ)], [-sin(θ),cos(θ)]
    a,b = coordinates[1:2]

    N = length(weights) # surface modes
    Rα = sum( n-> weights[n]*exp(-(α-(n/N-1/2)π)^2/variance), 1:N )

    model .= X₀ + (a+Rα)*cos(α-θ)*u + (b+Rα)*sin(α-θ)*v
    return model
end

function boundary!( model::AbstractVector, data::AbstractVector; kwargs...)
    α = atan(data[2],data[1])
    return boundary!(model,α; kwargs...)
end