function boundary!( model::AbstractVector, α::T; weights=[0.0], origin=[0.0,0.0], 
                    orientation=[75.0,0.0], variance=0.01) where T<:Number

    radius,rotation = orientation
    N = length(weights)
    r = radius + sum( n-> weights[n]*exp(-(α-(n/N-1/2)π)^2/variance), 1:N )

    model[1] = r*cos(α-rotation) + origin[1]
    model[2] = r*sin(α-rotation) + origin[2]

    return model
end

function boundary!( model::AbstractVector, data::AbstractVector; kwargs...)
    α = atan(data[2],data[1])
    return boundary!(model,α; kwargs...)
end