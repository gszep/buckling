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

@time begin #################################### load images and fit spherical coordinates

    image_axes = File("data/E2/M.tif", z=40:60, t=1:120, Δr=[0.37μm,0.37μm,1.5μm], Δt=15s).axes
    images = AxisArray( StackedView(

        File("data/E2/M.tif", z=40:60, t=1:120, Δr=[0.37μm,0.37μm,1.5μm], Δt=15s),
        File("data/E2/U.tif", z=40:60, t=1:120, Δr=[0.37μm,0.37μm,1.5μm], Δt=15s)

    ), Axis{:channel}([:M, :U]), image_axes...)

    binaryImages = AxisArray( Array{Bool}(undef,size(images)), images.axes )
    binarize!(binaryImages,images,Otsu())

    ############################################ fit spherical coordinate system
    parameters = Dict( :weights => zeros(50), :coordinates => zeros(5) )
    coordinates = [ zeros(5) for t ∈ 1:nimages(images) ]
    angles = range(-π,π,length=500)

    boundaries = Dict(
        :U=>[ zeros(length(angles),2) for t ∈ 1:nimages(images) ],
        :M=>[ zeros(length(angles),2) for t ∈ 1:nimages(images) ]
    )

    for channel ∈ keys(boundaries)
        for t ∈ 1:nimages(images)

            dataₜ = targets(binaryImages[Axis{:channel}(channel),Axis{:t}(t)])
            ellipse!( coordinates[t], map(x->first(x)[1],dataₜ), map(x->first(x)[2],dataₜ), radius_offset = (channel == :U) | (channel == :iM) ? 15.0 : 0.0 )

            parameters[:coordinates] = coordinates[t]
            for k ∈ 1:length(angles) # record solution
                boundaries[channel][t][k,:] = copy(boundary!( boundaries[channel][t][k,:], angles[k]; parameters... ))
            end
        end
    end
end

begin
    ############################################################### plotting
    scene, layout = layoutscene(resolution = (800,800))
    slider = layout[1,1] = LSlider( scene, range=1:nimages(images), startvalue=26 )

    menu = layout[2,1] = LMenu(scene, options = zip(["U","M"],[:U,:M]))
    menu.selection.val = :M

    ax = layout[3,1] = LScene(scene, camera = cam3d!)
    t,channel = slider.value, menu.selection

    ######################### show image volume
    imageₜ = @lift( images[Axis{:channel}($channel),Axis{:t}($t)].data )
    x,y,z,_ = ( range( extrema(ustrip(ax.val))..., length=2) for ax ∈ image_axes )
    volume!( ax, x,y, minimum(z) ≠ maximum(z) ? z : range(-1,1,length=2), imageₜ)

    ######################### boundary model
    boundaryₜ = @lift( boundaries[$channel][$t] )
    lines!( ax, boundaryₜ, linewidth=3, color=:gold )

    ax = layout[3,2] = LAxis(scene, xlabel="Time [mins]", ylabel="Circumference [μm]")
    boundaryLength = [ sum(sqrt.(sum(abs2.(diff(boundary,dims=1)),dims=2))) for boundary ∈ boundaries[:U] ]
    u = lines!(ax, map( t->ustrip(uconvert(mins,t)), timeaxis(images) ), boundaryLength,
        markersize=4, color=:gold)

    boundaryLength = [ sum(sqrt.(sum(abs2.(diff(boundary,dims=1)),dims=2))) for boundary ∈ boundaries[:M] ]
    m = lines!(ax, map( t->ustrip(uconvert(mins,t)), timeaxis(images) ), boundaryLength,
        markersize=4, color=:blue)

    tMins = @lift( ustrip(uconvert(mins,timeaxis(images)[$t])) )
    vlines!(ax,tMins,linestyle=:dot)

    layout[3,3] = LLegend(scene, [u,m], ["U", "M"])
    printstyled(color=:green,"Preprocessing Done")
    RecordEvents( scene, "output" )
    scene
end

@time begin ################################################################## optimise boundary
    gradients =  Dict( :weights => zeros(50) )

    #################################### de-noise binarization
    region = Tuple(1:ndims(images))
    opening!(binaryImages,region); closing!(binaryImages,region)

    for t ∈ 1:nimages(images)

        dataₜ = targets(binaryImages[Axis{:z}(1),Axis{:t}(t)])
        parameters[:coordinates] = coordinates[t]

        optimiser = Flux.Momentum(1.0)
        for _ ∈ UnitRange{Int}(1,50)

            gradient!( gradients[:weights], θ -> loss( dataₜ; weights=θ, coordinates=parameters[:coordinates] ), parameters[:weights] )
            Flux.update!(optimiser,parameters[:weights],gradients[:weights])
        end

        for k ∈ 1:length(angles) # record solution
            boundaries[t][k,:] = copy(boundary!( boundaries[t][k,:], angles[k]; parameters... ))
        end
    end

    ############################################################### plotting
    boundaryₜ = lift(tSlider.value) do t boundaries[t] end
    lines!( ax, boundaryₜ, linewidth=3, color=:gold )

    printstyled(color=:blue,"Optimisation Done")
    RecordEvents( scene, "output" )
    scene
end

begin ################################################################## display timeseries calculations

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
    scene
end
