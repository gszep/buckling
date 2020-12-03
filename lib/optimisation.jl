using ReverseDiff
include("./patches/StaticArrays.jl")

function loss( x::StaticVector{2,T}; kwargs... ) where T<:Number
    return log(norm( x - boundary(atan(x[2],x[1]);kwargs...) ) )
end

function loss( data::Vector{<:StaticVector}; kwargs... )
    return sum( x->loss(x;kwargs...), data ) / length(data)
end

function ∇loss( data::Vector{<:StaticVector}; weights=[0.0], origin=[0.0,0.0], radius=75.0, rotation=0.0, kwargs... )
    return ReverseDiff.gradient( θ->loss(data;weights=θ,origin=origin,radius=radius,rotation=rotation, kwargs...), weights )
end



# function loss( x::T, y::T; kwargs... ) where T<:Quantity
#     x = images.axes[axisdim(images,Axis{:x})].val[index.I[1]]
#     y = images.axes[axisdim(images,Axis{:y})].val[index.I[2]]
#     return norm( [x,y] - boundary(atan(y,x);kwargs...) )
# end

# function loss( indexes::Vector{CartesianIndex{2}}; kwargs... )
#     return sum( idx -> loss(idx;kwargs...), indexes ) / length(indexes)
# end