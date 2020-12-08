using ForwardDiff: gradient!

function loss( x::Tuple{<:AbstractVector,<:AbstractVector}; kwargs... )
    data,model = x
    boundary!(model,data;kwargs...)
    return norm( data - model )^2
end

function loss( data::Vector{<:Tuple{<:AbstractVector,<:AbstractVector}}; kwargs... )
    return sum( x->loss(x;kwargs...), data ) / length(data)
end

function ellipse!( coordinates::AbstractVector, x::AbstractVector,y::AbstractVector;
         constraints = Array(sparse([2,1,3],[2,3,1],[1.0,-2.0,-2.0],6,6)),
         radius_offset = 0.0 )

    features = @. [ x^2 x*y y^2 x y one(x) ]
    scatter = features'features

    eig = eigen!(scatter,constraints)
    optimum = findfirst( @. ~isinf(eig.values) )
    A,B,C,D,E,F = eig.vectors[:,optimum] # cartesian coordinates

    ######################## to canonical coordinates
    ΔR = √((A-C)^2+B^2)
    Δ = B^2-4A*C

    # eccentricity
    coordinates[1] = radius_offset-√( 2*( A*E^2 + C*D^2 - B*D*E + Δ*F )*( A + C + ΔR ) ) / Δ
    coordinates[2] = radius_offset-√( 2*( A*E^2 + C*D^2 - B*D*E + Δ*F )*( A + C - ΔR ) ) / Δ

    # origin and orientation
    coordinates[3] = (2C*D-B*E)/Δ
    coordinates[4] = (2A*E-B*D)/Δ
    coordinates[5] = atan((C-A-ΔR)/B)

    return coordinates
end