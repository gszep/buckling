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
    name = "E11/iM.tif"
    images = File(joinpath("data",name), z=40:60, t=1:120, Δr=[0.37μm,0.37μm,1.5μm], Δt=15s)

    binaryImages = AxisArray( Array{Bool}(undef,size(images)), images.axes )
    binarize!(binaryImages,images,Otsu())
    region = Tuple(1:ndims(images))

    #opening!(binaryImages,region); closing!(binaryImages,region) # convert to sparse(binaryImages) ?
    ############################################ fit spherical coordinate system
    coordinates = [ zeros(5) for t ∈ 1:nimages(images) ]
    for t ∈ 1:nimages(images)
        dataₜ = targets(binaryImages[Axis{:t}(t)])
        ellipse!( coordinates[t], map(x->first(x)[1],dataₜ), map(x->first(x)[2],dataₜ), radius_offset = contains(name,"U") | contains(name,"iM") ? 15.0 : 0.0 )
    end

    #################################### de-noise binarization
    opening!(binaryImages,region); closing!(binaryImages,region) # convert to sparse(binaryImages) ?
    printstyled(color=:green,"Preprocessing Done")
end

@time begin ################################################################## optimise boundary
    parameters = Dict( :weights => zeros(50), :coordinates => zeros(5) )
    gradients =  Dict( :weights => zeros(50) )

    angles = range(-π,π,length=500)
    boundaries = [ zeros(length(angles),2) for t ∈ 1:nimages(images) ]

    optimiser = Flux.Momentum(1.0)
    for t ∈ 1:nimages(images)

        dataₜ = targets(binaryImages[Axis{:z}(1),Axis{:t}(t)])
        parameters[:coordinates] = coordinates[t]

        for _ ∈ UnitRange{Int}(1,50)

            gradient!( gradients[:weights], θ -> loss( dataₜ; weights=θ, coordinates=parameters[:coordinates] ), parameters[:weights] )
            Flux.update!(optimiser,parameters[:weights],gradients[:weights])
        end

        for k ∈ 1:length(angles) # record solution
            boundaries[t][k,:] = copy(boundary!( boundaries[t][k,:], angles[k]; parameters... ))
        end
    end
    printstyled(color=:blue,"Optimisation Done")
end

begin ################################################################## display image timeseries
    scene, layout = layoutscene(resolution = (1024,512))
    tSlider = layout[1,1] = LSlider( scene, range=1:nimages(images), startvalue=26 )
    ax = layout[2,1] = LScene(scene, camera = cam3d!)

    ######################### show image volume
    imageₜ = lift(tSlider.value) do t images[Axis{:t}(t)].data end
    x,y,z,_ = ( range( extrema(ustrip(ax.val))..., length=2) for ax ∈ images.axes )
    volume!( ax, x,y, minimum(z) ≠ maximum(z) ? z : range(-1,1,length=2), imageₜ)

    # dataₜ = lift(tSlider.value) do t map(x->Point(first(x)...),targets(binaryImages[Axis{:t}(t)])) end
    # meshscatter!( ax, dataₜ, color=:gold, markersize=1)

    ######################### boundary model
    boundaryₜ = lift(tSlider.value) do t boundaries[t] end
    lines!( ax, boundaryₜ, linewidth=3, color=:gold )

    ax = layout[2,2] = LAxis(scene, xlabel="Time [mins]", ylabel="Circumference [μm]")
    boundaryLength = [  sum(sqrt.(sum(abs2.(diff(boundary,dims=1)),dims=2))) for boundary ∈ boundaries ]
    scatter!(ax, map( t->ustrip(uconvert(mins,t)), timeaxis(images) ), boundaryLength,
        markersize=4, color=:gold)

    area = [ 3*a*b*π / 10^2 for (a,b,x,y,θ) ∈ coordinates ]
    scatter!(ax, map( t->ustrip(uconvert(mins,t)), timeaxis(images) ), area,
        markersize=4, color=:blue)

    tValue = lift(tSlider.value) do t ustrip(uconvert(mins,timeaxis(images)[t])) end
    vlines!(ax,tValue,linestyle=:dot)

    ######################### render all
    RecordEvents( scene, "output" )
end
