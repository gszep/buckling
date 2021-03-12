function ellipse!( parameters::AbstractVector, x::AbstractVector,y::AbstractVector;
         constraints = Array(sparse([2,1,3],[2,3,1],[1.0,-2.0,-2.0],6,6)),
         radius_offset = 0.0 )

    features = @. [ x^2 x*y y^2 x y one(x) ]
    scatter = features'features

    eig = eigen!(scatter,constraints)
    optimum = findfirst( @. ~isinf(eig.values) )
    A,B,C,D,E,F = eig.vectors[:,optimum] # cartesian coordinates

    ######################## to ellipse parameters
    ΔR = √((A-C)^2+B^2)
    Δ = B^2-4A*C

    # eccentricity
    parameters[1] = radius_offset-√( 2*( A*E^2 + C*D^2 - B*D*E + Δ*F )*( A + C + ΔR ) ) / Δ
    parameters[2] = radius_offset-√( 2*( A*E^2 + C*D^2 - B*D*E + Δ*F )*( A + C - ΔR ) ) / Δ

    # origin and orientation
    parameters[3] = (2C*D-B*E)/Δ
    parameters[4] = (2A*E-B*D)/Δ
    parameters[5] = atan((C-A-ΔR)/B)

    return parameters
end

function ellipse( A::AxisArray{<:Bool}; kwargs... )
    parameters = Dict([ channel=>[ zeros(5) for t ∈ 1:nimages(A) ] for channel ∈ first(A.axes) ])

    for channel ∈ first(A.axes)
        for t ∈ 1:nimages(A)

            dataₜ = targets(A[Axis{:channel}(channel),Axis{:t}(t)])
            ellipse!( parameters[channel][t], map(x->x[1],dataₜ), map(x->x[2],dataₜ);
                radius_offset = (channel == :U) | (channel == :iM) ? 15.0 : 0.0,
                kwargs... )
        end
    end

    return parameters
end

function surface_modes( A::AxisArray{<:Bool}, parameters::Dict{Symbol,Vector{Vector{T}}}; N = 50, variance=0.01, kwargs... ) where T<:Number
    weights = Dict([ channel=>[ zeros(N) for t ∈ 1:nimages(A) ] for channel ∈ first(A.axes) ])
    residuals = Dict([ channel=>[ Float64[] for t ∈ 1:nimages(A) ] for channel ∈ first(A.axes) ])
    angles = Dict([ channel=>[ Float64[] for t ∈ 1:nimages(A) ] for channel ∈ first(A.axes) ])

    for channel ∈ first(A.axes)
        for t ∈ 1:nimages(A)

            dataₜ = targets(A[Axis{:channel}(channel),Axis{:t}(t),Axis{:z}(1)])
            R,X = parameters[channel][t][1:2], parameters[channel][t][3:4]

            R = map( x -> norm(x-X)-sum(R)/2, dataₜ )
            α = map( x -> atan(x[2]-X[2],x[1]-X[1]), dataₜ )
            Φ = hcat(map( n-> exp.(-(α.-(n/N-1/2)π).^2/variance), 1:N )...)

            weights[channel][t] = (Φ'Φ+5I(N))\Φ'R
            residuals[channel][t] = R
            angles[channel][t] = α
        end
    end

    return weights,residuals,angles
end