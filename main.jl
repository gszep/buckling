begin #################################### required libs
    using LinearAlgebra,SparseArrays
    using Images,ImageBinarization

    using ImageMorphology: opening!,closing!
    using AbstractPlotting.MakieLayout
    using GLMakie,Makie

    include("lib/boundary.jl")
    include("lib/optimisation.jl")
    include("lib/utils.jl")

    printstyled("Libraries Loaded\n")
end

@time begin #################################### load images 
    images = Files("data/E2", z=40:60, t=1:120, Δr=[0.37μm,0.37μm,1.5μm], Δt=15s)
    binaryImages = AxisArray( Array{Bool}(undef,size(images)), images.axes )

    ############################################ fit ellipsiodal coordinate system
    binarize!(binaryImages,images,Otsu())
    parameters = ellipse( binaryImages )
    weights = nothing

    #################################### de-noise binarization
    region = Tuple(2:ndims(images))
    opening!(binaryImages,region); closing!(binaryImages,region)
    printstyled(color=:green,"Preprocessing Done")
end

begin ############################################################### plotting
    scene, layout = layoutscene(resolution = (700,700))
    ax = layout[3,1] = LScene(scene, camera = cam3d!)

    slider = layout[1,1] = LSlider( scene, camera = cam3d!, range=1:nimages(images), startvalue=26 )
    menu = layout[2,1] = LMenu(scene, options = [ (String(channel),channel) for channel ∈ first(images.axes) ] )

    menu.selection.val = first(first(images.axes))
    t,channel = slider.value, menu.selection

    ######################### show image volume
    imageₜ = @lift( images[Axis{:channel}($channel),Axis{:t}($t)].data )
    x,y,z,_ = ( range( extrema(ustrip(ax.val))..., length=2) for ax ∈ images.axes if eltype(ax) <: Quantity  )
    volume!( ax, x,y, minimum(z) ≠ maximum(z) ? z : range(-1,1,length=2), imageₜ)

    ######################### boundary model
    boundaries = boundary(parameters,weights)
    lines!( ax, @lift( boundaries[$channel][$t] ), linewidth=3, color=:gold )

    ax = layout[3,2] = LAxis(scene, xlabel="Time [mins]", ylabel="Circumference [μm]")
    circumferences = []
    colors = [:gold,:blue]
    for (i,channel) ∈ enumerate(first(images.axes))

        boundaryLength = [ sum(sqrt.(sum(abs2.(diff(boundary,dims=1)),dims=2))) for boundary ∈ boundaries[channel] ]
        line = lines!(ax, map( t->ustrip(uconvert(mins,t)), timeaxis(images) ), boundaryLength, color=colors[i] )
        push!(circumferences,line)
    end

    tMins = @lift( ustrip(uconvert(mins,timeaxis(images)[$t])) )
    vlines!(ax,tMins,linestyle=:dot)

    layout[3,3] = LLegend(scene, circumferences, [String(channel) for channel ∈ first(images.axes)])
    RecordEvents( scene, "output" )
    scene
end

@time begin ################################################################## optimise boundary
    weights = surface_modes(binaryImages,parameters)
    printstyled(color=:green,"Surface Mode Regression Done")
end
