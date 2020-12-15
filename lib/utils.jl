using HDF5,Images
using ImageAxes: Axis

using SimpleTraits
@traitimpl TimeAxis{Axis{:t}}

using Unitful
const s = Unitful.s
const mins = Unitful.minute
const μm = Unitful.μm

function File(path::String; t=:, x=:, y=:, z=:, Δt=15s, Δr=(0.37μm,0.37μm,1.5μm) )
    filename,ext = splitext(path)
    c = 1 # number of channels

    if occursin("h5",ext)
        A = h5read(path,"data",(:,x,y,z,t))
        @assert((size(A,1)==1)&(ndims(A)==5),"allowed formats: single channel TIF (x,y,t) or HDF5 (x,y,z,t)")
        z = 1:size(A,4)

    elseif occursin("tif",ext)
        A = channelview(load(path))
        if ndims(A)==3 A = A[[CartesianIndex()],:,:,:] elseif ndims(A)==4 A = maximum(A,dims=1) end
        @assert(ndims(A)==4,"allowed formats: single channel TIF (x,y,t) or HDF5 (x,y,z,t)")

        A = A[:,x,y,t]
        z = [CartesianIndex()] # single z-slice

    else 
        throw("allowed formats: single channel TIF (x,y,t) or HDF5 (x,y,z,t)")
    end

    Δx, Δy, Δz = Δr
    images = A[c,:,:,z,:]

    return AxisArray( colorview(Gray,normedview(images)),

        Axis{:x}( size(images,1)>1 ? range(-Δx*size(images,1)/2,Δx*size(images,1)/2,length=size(images,1)) : 0*Δx ),
        Axis{:y}( size(images,2)>1 ? range(-Δy*size(images,2)/2,Δy*size(images,2)/2,length=size(images,2)) : 0*Δy ),
        Axis{:z}( size(images,3)>1 ? range(-Δz*size(images,3)/2,Δz*size(images,3)/2,length=size(images,3)) : 0*Δz ),
        Axis{:t}( size(images,4)>1 ? range(0*Δt,Δt*size(images,4),length=size(images,4)) : 0*Δt )
    )
end

function Files(path::String; t=:, x=:, y=:, z=:, Δt=15s, Δr=(0.37μm,0.37μm,1.5μm) )

    @assert( ~occursin(".",path), "expecting folder of .tif files for each channel")
    files = filter( file->endswith(file,".tif"), joinpath.( path, readdir(path)) )
    image_axes = File( first(files), t=t, x=x, y=y, z=z, Δt=Δt, Δr=Δr ).axes

    return AxisArray( StackedView([ File(file, t=t, x=x, y=y, z=z, Δt=Δt, Δr=Δr ) for file ∈ files ]...),
        Axis{:channel}([ Symbol(basename(first(splitext(file)))) for file ∈ files ]), image_axes...)
end

function targets( A::AxisArray{<:Bool}; downSample::Integer=1 )
    targets = Vector{Tuple{Vector{Float64},Vector{Float64}}}()

    for index ∈ findall(A)[1:downSample:end]
        push!(targets, ( [A.axes[k][index[k]].val for k ∈ 1:ndims(A)], zeros(ndims(A)) ) )
    end

    @assert(length(targets)<10^5,"too many targets; increase downSample")
    return targets
end