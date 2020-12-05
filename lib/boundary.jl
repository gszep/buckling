function boundary( α::T; weights=[0.0], origin=[0.0,0.0], orientation=[75.0,0.0], variance=0.01) where T<:Number
    radius,rotation = orientation

    μ = length(weights) > 1 ? range(-π/2,π/2,length=length(weights)) : 0.0
    r = radius + sum( n-> weights[n]*exp(-(α-μ[n])^2/variance), 1:length(weights) )

    return r*[cos(α-rotation),sin(α-rotation)] + origin
end