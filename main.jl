begin #################################### required libs
    using LinearAlgebra,SparseArrays,StaticArrays
    using Images,ImageBinarization

    using ImageMorphology: opening!
    using AbstractPlotting.MakieLayout
    using GLMakie,Makie

    include("lib/boundary.jl")
    include("lib/optimisation.jl")
    include("lib/utils.jl")
end

begin #################################### load and binarize image
    name = "E2/M.tif"
    images = File(joinpath("data",name), Δr=[0.37μm,0.37μm,1.5μm], Δt=15s)

    binaryImages = AxisArray( Array{Bool}(undef,size(images)), images.axes )
    binarize!(binaryImages,images,Otsu())

    #################################### de-noise binarization
    opening!(binaryImages,Tuple(1:ndims(images))) # convert to sparse(binaryImages) ?
    printstyled(color=:green,"Preprocessing Done")
end

begin ################################################################## display image timeseries

    ######################### define time slider
    tSlider = slider( 1:nimages(images), sliderlength=700, position=(50,0), start=26,
        valueprinter=tSlider->"$(round(mins,step(timeaxis(images))*tSlider,digits=2))" )
    imageₜ = lift(tSlider[end].value) do i .~binaryImages[Axis{:t}(i),Axis{:z}(1)] end

    ######################### image sequence
    scene, layout = layoutscene(resolution = (800, 800))
    ax = layout[1, 1] = LAxis(scene, xlabel = "μm", ylabel = "μm")
    
    image!( ax,
        map(x-> uconvert(μm,x).val, images.axes[axisdim(images,Axis{:x})]),
        map(y-> uconvert(μm,y).val, images.axes[axisdim(images,Axis{:y})]),
        imageₜ, color=:black )

    ######################### boundary model
    angles = range(-π,π,length=500)
    weights = randn(50)

    boundary_positions = hcat( boundary.(angles; weights=weights*μm, radius=75.0μm, origin=zeros(2)μm )... )'
    lines!( ax, ustrip(uconvert.(μm,boundary_positions)), linewidth=3, color=:gold )

    ######################### render all
    RecordEvents( hbox( scene, vbox(tSlider), parent = Scene(resolution=(800, 800)) ), "output" )
end

data = targets(binaryImages[Axis{:z}(1),Axis{:t}(26)])

loss(data;weights=weights)
∇loss(data;weights=weights)