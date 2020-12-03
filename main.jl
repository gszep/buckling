begin #################################### required libs
    using LinearAlgebra,SparseArrays
    using Images,ImageBinarization

    using ImageMorphology: opening!
    using AbstractPlotting.MakieLayout
    using GLMakie,Makie

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

function boundary( α::T, W::Vector{U}; origin::Vector{U}=[0.0μm,0.0μm], rotation::T=0.0, variance::T=0.01) where {T<:Number,U<:Quantity}
    μ = range(-π/2,π/2,length=length(W)-1)

    r = W[1] + sum( n-> W[n]*exp(-(α-μ[n-1])^2/variance), 2:length(W) )
    return r*[cos(α-rotation),sin(α-rotation)] + origin
end

begin ################################################################## display image timeseries
    parameters = [75μm; randn(100)μm]
    ######################### define time slider
    tSlider = slider( 1:nimages(images), sliderlength=700, position=(50,0), start=26,
        valueprinter=tSlider->"$(round(mins,step(timeaxis(images))*tSlider,digits=2))" )
    imageₜ = lift(tSlider[end].value) do i .~binaryImages[Axis{:t}(i),Axis{:z}(1)] end

    #width,height,angle = range(-1,1,length=width(images)),range(-1,1,length=height(images)),
    ######################### image sequence

    scene, layout = layoutscene(resolution = (800, 800))
    ax = layout[1, 1] = LAxis(scene, xlabel = "μm", ylabel = "μm")
    
    image!( ax,
        map(x-> uconvert(μm,x).val, images.axes[axisdim(images,Axis{:x})]),
        map(y-> uconvert(μm,y).val, images.axes[axisdim(images,Axis{:y})]),
        imageₜ, color=:black )

    ######################### boundary model
    angles = range(-π,π,length=500)
    boundary_positions = hcat(boundary.(angles,Ref(parameters))...)'
    lines!( ax, ustrip(uconvert.(μm,boundary_positions)), linewidth=3, color=:gold )

    ######################### render all
    RecordEvents( hbox( scene, vbox(tSlider), parent = Scene(resolution=(800, 800)) ), "output" )
end

targets = findall(binaryImages[:,:,1,26])[1:10:end]
loss(targets,parameters)




# img = load("data/E2/M_allZ.tif")

# img = load("data/E2/U.tif")
# img[:,:,50]

# img = load("data/E5/Mb.tif")
# img[:,:,end]

# img = load("data/E11/iM.tif")
# img[:,:,end]