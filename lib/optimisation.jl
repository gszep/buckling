function loss( index::CartesianIndex{2}, θ::AbstractVector{T} ) where T<:Number
    x = images.axes[axisdim(images,Axis{:x})].val[index.I[1]]
    y = images.axes[axisdim(images,Axis{:y})].val[index.I[2]]
    return norm( [x,y] - boundary( atan(y,x), θ))
end

function loss( indexes::Vector{CartesianIndex{2}}, θ::AbstractVector{T} ) where T<:Number
    return sum( idx -> loss(idx,θ), indexes ) / length(indexes)
end