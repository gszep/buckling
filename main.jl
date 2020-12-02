using Images,ImageBinarization
using GLMakie,Makie

M = binarize(load("data/E2/M.tif"), Otsu())
M = opening(opening(M))

begin
    t = slider(1:size(M,3),sliderlength=400,position=(50,0),start=26)
    Mₜ = lift(t[end].value) do i M[:,:,i] end

    RecordEvents( hbox( image(Mₜ,axes=false), vbox(t), 
        parent = Scene(resolution=(500, 500))),
    "output" )
end


# img = load("data/E2/M_allZ.tif")

# img = load("data/E2/U.tif")
# img[:,:,50]

# img = load("data/E5/Mb.tif")
# img[:,:,end]

# img = load("data/E11/iM.tif")
# img[:,:,end]