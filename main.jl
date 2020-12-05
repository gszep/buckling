begin #################################### required libs
    using LinearAlgebra,SparseArrays,StaticArrays
    using Images,ImageBinarization

    using ImageMorphology: opening!,closing!
    using AbstractPlotting.MakieLayout
    using GLMakie,Makie

    include("lib/boundary.jl")
    include("lib/optimisation.jl")
    include("lib/utils.jl")
end

begin #################################### load and binarize image
    name = "E2/M.tif"
    images = File(joinpath("data",name), z=40:60, Δr=[0.37μm,0.37μm,1.5μm], Δt=15s)

    binaryImages = AxisArray( Array{Bool}(undef,size(images)), images.axes )
    binarize!(binaryImages,images,Otsu())
    region = Tuple(1:ndims(images))

    #################################### de-noise binarization
    opening!(binaryImages,region); closing!(binaryImages,region) # convert to sparse(binaryImages) ?
    printstyled(color=:green,"Preprocessing Done")
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
    angles = range(-π,π,length=500)
    boundary_positions = hcat( boundary.(angles; weights=weights, radius=75.0, origin=[0.0,5] )... )'
    lines!( ax, boundary_positions, linewidth=3, color=:gold )

    ######################### render all
    RecordEvents( scene, "output" )
end

data = targets(binaryImages[Axis{:z}(1),Axis{:t}(26)])
using Flux: Momentum,update!

optimiser = Momentum(1.0)
for i ∈ 1:100
    gradients = ∇loss(data;weights=weights,radius=75.0,origin=[0.0,5])
    update!(optimiser,weights,gradients)
end