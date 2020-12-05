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
    name = "E2/M.h5"
    images = File(joinpath("data",name), z=48:52, Δr=[0.37μm,0.37μm,1.5μm], Δt=15s)

    binaryImages = AxisArray( Array{Bool}(undef,size(images)), images.axes )
    binarize!(binaryImages,images,Otsu())
    region = Tuple(1:ndims(images)-1)

    #################################### de-noise binarization
    opening!(binaryImages,region); closing!(binaryImages,region) # convert to sparse(binaryImages) ?
    printstyled(color=:green,"Preprocessing Done")
end

begin ################################################################## display image timeseries
    scene, layout = layoutscene(resolution = (800, 800))
    tSlider = layout[1,1] = LSlider( scene, range=1:nimages(images), startvalue=26 )
    ax = layout[2,1] = LScene(scene, camera = cam3d!, raw = false, xlabel = "μm", ylabel = "μm")

    ######################### scatter binarized image
    dataₜ = lift(tSlider.value) do t map( x-> Point(x[1],x[2],x[3]), targets(binaryImages[Axis{:t}(t)]) ) end
    meshscatter!( ax, dataₜ, color = :darkblue, markersize = 1 )

    ######################### boundary model
    angles = range(-π,π,length=500)
    weights = randn(50)

    boundary_positions = hcat( boundary.(angles; weights=weights*μm, radius=75.0μm, origin=zeros(2)μm )... )'
    lines!( ax, ustrip(uconvert.(μm,boundary_positions)), linewidth=3, color=:gold )

    ######################### render all
    RecordEvents( scene, "output" )
end

data = targets(binaryImages[Axis{:z}(1),Axis{:t}(26)])

loss(data;weights=weights)
∇loss(data;weights=weights)