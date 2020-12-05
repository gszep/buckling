begin #################################### required libs
    using LinearAlgebra,SparseArrays
    using Images,ImageBinarization
    using Flux

    using ImageMorphology: opening!,closing!
    using AbstractPlotting.MakieLayout
    using GLMakie,Makie

    include("lib/boundary.jl")
    include("lib/optimisation.jl")
    include("lib/utils.jl")
end

@time begin #################################### load and binarize image
    name = "E2/M.tif"
    images = File(joinpath("data",name), z=40:60, Δr=[0.37μm,0.37μm,1.5μm], Δt=15s)

    binaryImages = AxisArray( Array{Bool}(undef,size(images)), images.axes )
    binarize!(binaryImages,images,Otsu())
    region = Tuple(1:ndims(images))

    #################################### de-noise binarization
    opening!(binaryImages,region); closing!(binaryImages,region) # convert to sparse(binaryImages) ?
    printstyled(color=:green,"Preprocessing Done")
end

@time begin ################################################################## optimise boundary
    parameters = Dict( :weights => zeros(50), :origin => [0.0,5.0], :orientation=>[75.0,0.0] )
    gradients =  Dict( :weights => zeros(50), :origin => [0.0,0.0], :orientation=>[0.0,0.0] )

    angles = range(-π,π,length=500)
    boundaries = [ zeros(length(angles),2) for t ∈ 1:nimages(images) ]

    optimiser = Flux.Momentum(1.0)
    for t ∈ 1:nimages(images)
        dataₜ = targets(binaryImages[Axis{:z}(1),Axis{:t}(t)])

        for _ ∈ ( t>1 ? UnitRange{Int}(1,10) : UnitRange{Int}(1,100) )

            gradients!(gradients,dataₜ;parameters...)
            Flux.update!(optimiser,parameters,gradients)
        end

        for k ∈ 1:length(angles) # record solution
            boundaries[t][k,:] = copy(boundary!( boundaries[t][k,:], angles[k]; parameters... ))
        end
    end
    printstyled(color=:blue,"Optimisation Done")
end

begin ################################################################## display image timeseries
    scene, layout = layoutscene(resolution = (800, 800))
    tSlider = layout[1,1] = LSlider( scene, range=1:nimages(images), startvalue=26 )
    ax = layout[2,1] = LScene(scene, camera = cam3d!)

    ######################### show image volume
    imageₜ = lift(tSlider.value) do t images[Axis{:t}(t)].data end
    x,y,z,_ = ( range( extrema(ustrip(ax.val))..., length=2) for ax ∈ images.axes )
    volume!( ax, x,y, minimum(z) ≠ maximum(z) ? z : range(-1,1,length=2), imageₜ)

    ######################### boundary model
    boundaryₜ = lift(tSlider.value) do t boundaries[t] end
    lines!( ax, boundaryₜ, linewidth=3, color=:gold )

    ######################### render all
    RecordEvents( scene, "output" )
end
