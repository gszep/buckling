using Images,ImageBinarization
using GLMakie,Makie

M = binarize(load("data/E2/M.tif"), Otsu())
M = opening(opening(M))

begin
    t = slider(1:size(M,3),sliderlength=700,position=(50,0),start=26,valueprinter=x->"$(round(15*x/60,digits=1)) min")
    Mₜ = lift(t[end].value) do i .~M[:,:,i] end

    Figure = image(range(-1,1,length=size(M,1)),range(-1,1,length=size(M,2)),Mₜ,show_axis=false)

    θ = range(-π,π,length=100)
    lines!(3cos.(θ)/4,3sin.(θ)/4,color=:gold,linewidth=3)

    RecordEvents( hbox( Figure, vbox(t), 
        parent = Scene(resolution=(800, 800))),
    "output" )
end


# img = load("data/E2/M_allZ.tif")

# img = load("data/E2/U.tif")
# img[:,:,50]

# img = load("data/E5/Mb.tif")
# img[:,:,end]

# img = load("data/E11/iM.tif")
# img[:,:,end]