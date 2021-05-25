begin
    import Pkg
    Pkg.activate(".")
end

begin #################################### required libs
    using ImageMorphology: opening!,closing!
    using GLMakie,CairoMakie
    using DataFrames

    using LinearAlgebra,SparseArrays
    using Images,ImageBinarization

    include("lib/boundary.jl")
    include("lib/optimisation.jl")
    include("lib/utils.jl")

    printstyled("Libraries Loaded\n")
end

@time begin #################################### load images 
    egg = "E2"

    images = Files("data/$egg", z=40:60, t=1:120, Δr=[0.37μm,0.37μm,1.5μm], Δt=15s)
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

@time begin ################################################################## optimise boundary
    weights,residuals,angles = surface_modes(binaryImages,parameters)
    printstyled(color=:green,"Surface Mode Regression Done")
end

begin ############################################################### show movie
    GLMakie.activate!()
    movie = Figure(resolution = (500,700))

    ax = LScene( movie[2,1:2], scenekw = ( camera = cam3d!, show_axis=false, show_grid=false ))
    colors = [ parse(RGBA,"#0000FF"), parse(RGBA,"#008000")]

    ######################### show image volume
    timeSlider = labelslider!(movie, "Time", 1:nimages(images); format = i -> "$(round(i*step(timeaxis(images)).val/60,digits=2)) mins", startvalue=26)
    t, movie[3,1] = timeSlider.slider, timeSlider.layout

    x,y,z,_ = ( range( extrema(ustrip(ax.val))..., length=2) for ax ∈ images.axes if eltype(ax) <: Quantity  )
    for (i,channel) ∈ enumerate(first(images.axes))

        volume!( ax, x,y, minimum(z) ≠ maximum(z) ? z : range(-1,1,length=2),
            @lift( images[Axis{:channel}(channel),Axis{:t}($(t.value))].data ),
            colormap = cgrad( [ RGBA(0,0,0,0), colors[i] ], [0,eps(),1] ) )
    end

    ######################### boundary model
    boundaries = boundary(parameters,weights)
    boundaries = Dict([ channel => map( x->x[150:end-150,:], boundaries[channel]) for channel ∈ keys(boundaries) ]) # drop out butthole

    legend = []
    for (i,channel) ∈ enumerate(first(images.axes))
    
        line = lines!( ax, @lift( boundaries[channel][$(t.value)] ), linewidth=3, color=colors[i] )
        push!(legend,line)
    end

    movie[3,2] = Legend(movie, legend,
        [String(channel) for channel ∈ first(images.axes)], framevisible = false)

    button = Button(movie[1,1:2], label = "Export")
    on(button.clicks) do event
        save("output/$egg@t=$(round(t.value[]*step(timeaxis(images)).val/60,digits=2))mins.png",movie)

        begin ############################################################### save statistics
            CairoMakie.activate!()
            svg = Figure(resolution = (1000,500))
        
            ######################### boundary length
            ax = svg[1,1] = CairoMakie.Axis(svg, xlabel="Time [mins]", ylabel="Vegetal Boundary Length [μm]", xgridvisible=false, ygridvisible=false)
            legend,lengths = [], DataFrame()

            time = map( t->ustrip(uconvert(mins,t)), timeaxis(images) )
            lengths[!,"time [mins]"] = time
            for (i,channel) ∈ enumerate(first(images.axes))
        
                boundaryLength = [ sum(sqrt.(sum(abs2.(diff(boundary,dims=1)),dims=2))) for boundary ∈ boundaries[channel] ]
                lengths[!,"$channel [μm]"] = boundaryLength

                line = lines!(ax, time, boundaryLength, color=colors[i] )
                push!(legend,line)
            end
        
            tMins = @lift( ustrip(uconvert(mins,timeaxis(images)[$(t.value)])) )
            vlines!(ax,tMins,linestyle=:dot)
        
            ######################### bucking
            ax = svg[1,2] = CairoMakie.Axis(svg, xlabel="Vegetal Position [radians]", ylabel="Buckling Amplitude [μm]", xgridvisible=false, ygridvisible=false)
            for (i,channel) ∈ enumerate(first(images.axes))
                scatter!(ax, @lift( [angles[channel][$(t.value)] residuals[channel][$(t.value)]] ),
                color=colors[i], markersize=2, strokewidth=0 )
            end
        
            svg[1,3] = Legend(svg, legend, [String(channel) for channel ∈ first(images.axes)], framevisible = false)
            save("output/$egg@t=$(round(t.value[]*step(timeaxis(images)).val/60,digits=2))mins.svg",svg)
            GLMakie.activate!()

            #################################################################### saving csv files
            save("output/$egg.boundary-lengths.csv",lengths)
            for channel ∈ first(images.axes)
                ellipses = DataFrame("time [mins]"=>time)

                ellipses[!,"a [μm]"] = map( x->x[1] , parameters[channel])
                ellipses[!,"b [μm]"] = map( x->x[2] , parameters[channel])

                ellipses[!,"x0 [μm]"] = map( x->x[3] , parameters[channel])
                ellipses[!,"y0 [μm]"] = map( x->x[4] , parameters[channel])
                ellipses[!,"θ [radians]"] = map( x->x[5] , parameters[channel])

                save("output/$egg.ellipse-parameters.$channel.csv",ellipses)
            end
        end
        printstyled(color=:blue,"Export Done\n")
    end
    display(movie)
    wait()
end