function loss( index::CartesianIndex{2}, θ::AbstractVector{T} ) where T<:Number
    x,y = width[index.I[1]], height[index.I[2]]
    return norm( [x,y] - boundary( atan(y,x), θ))
end

function loss( indexes::Vector{CartesianIndex{2}}, θ::AbstractVector{T} ) where T<:Number
    return sum( idx -> loss(idx,θ), indexes ) / length(indexes)
end

function hermiteGaussian(x::T, n::Integer, α::T) where T<:Number
    n == 0 ? (return one(T)) : nothing
    y = α*x

    Hₖ₋₁, Hₖ = one(T), 2y
    for k ∈ 1:n-1 # calculate polynomial using recurrance relation
        Hₖ₋₁, Hₖ = Hₖ, ( 2y*Hₖ - 2k*Hₖ₋₁ )
    end

    # weighted by gaussian
    return √α * Hₖ * exp(-y^2/2) / √abs(2^n * factorial(n) * √π)
end