begin #################################### required libs
    using LinearAlgebra,SparseArrays
    using Images,ImageBinarization

    using ImageMorphology: opening!,closing!
    using AbstractPlotting, GLMakie

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

    figure = Figure(resolution = (700,700))
    ax = LScene( figure[2,1], camera = cam3d!)
    colors = [ parse(RGBA,"#0000FF"), parse(RGBA,"#008000")]

    ######################### show image volume
    t = figure[1,2] = Slider( figure, range=1:nimages(images), startvalue=26 )
    x,y,z,_ = ( range( extrema(ustrip(ax.val))..., length=2) for ax ∈ images.axes if eltype(ax) <: Quantity  )

    for (i,channel) ∈ enumerate(first(images.axes))
        volume!( ax, x,y, minimum(z) ≠ maximum(z) ? z : range(-1,1,length=2),
            @lift( images[Axis{:channel}(channel),Axis{:t}($(t.value))].data ),
            colormap = cgrad( [ RGBA(0,0,0,0), colors[i] ], [0,eps(),1] ) )
    end

    ######################### boundary model
    boundaries = boundary(parameters,weights)
    boundaries = Dict([ channel => map( x->x[150:end-150,:], boundaries[channel]) for channel ∈ keys(boundaries) ]) # drop out butthole
    for (i,channel) ∈ enumerate(first(images.axes))
        lines!( ax, @lift( boundaries[channel][$(t.value)] ), linewidth=3, color=colors[i] )
    end

    ######################### circumference
    ax = figure[2,2] = AbstractPlotting.Axis(figure, xlabel="Time [mins]", ylabel="Circumference [μm]")
    circumferences = []
    for (i,channel) ∈ enumerate(first(images.axes))

        boundaryLength = [ sum(sqrt.(sum(abs2.(diff(boundary,dims=1)),dims=2))) for boundary ∈ boundaries[channel] ]
        line = lines!(ax, map( t->ustrip(uconvert(mins,t)), timeaxis(images) ), boundaryLength, color=colors[i] )
        push!(circumferences,line)
    end

    tMins = @lift( ustrip(uconvert(mins,timeaxis(images)[$(t.value)])) )
    vlines!(ax,tMins,linestyle=:dot)

    ######################### bucking
    ax = figure[2,3] = AbstractPlotting.Axis(figure, xlabel="Angle [radians]", ylabel="Amplitude [μm]")
    for (i,channel) ∈ enumerate(first(images.axes))
        scatter!(ax, @lift( [angles[channel][$(t.value)] residuals[channel][$(t.value)]] ),
        color=colors[i], markersize=2, strokewidth=0 )
    end

    figure[2,4] = Legend(figure, circumferences, [String(channel) for channel ∈ first(images.axes)])
    figure
end

@time begin ################################################################## optimise boundary
    weights,residuals,angles = surface_modes(binaryImages,parameters)
    printstyled(color=:green,"Surface Mode Regression Done")
end
