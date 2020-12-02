using LinearAlgebra,SparseArrays
using Images,ImageBinarization

using ImageMorphology: opening!
using GLMakie,Makie

include("lib/optimisation.jl")
include("lib/utils.jl")
const Δt = 15 # seconds

#################################### load and binarize image
#images = File("data/E2/M.h5",z=49:51)
images = File("data/E2/M.tif")

binaryImages = Array{Bool}(undef, size(images))
binarize!(binaryImages,images,Otsu())

#################################### apply sparse opening
# sparse(binarize(File("ters.h5"), Otsu()))
opening!(binaryImages) # scales like 19*Nz seconds

function boundary( α::T, W::AbstractVector{T}) where T<:Number
    k = 1:length(W)

    r = 3/4+W'hermiteGaussian.(α,k,5)
    return r*[cos(α),sin(α)]
end

parameters = [0.1,.1,.1,.1,.1,0,0,0,0,0.1]
begin
    t = slider( 1:size(M,4), sliderlength=700, position=(50,0), start=26, valueprinter=t->"$(round(Δt*t/60,digits=1)) min" )
    imageₜ = lift(t[end].value) do i .~binaryImages[:,:,1,i] end

    width,height,angle = range(-1,1,length=size(M,1)),range(-1,1,length=size(M,2)),range(-π,π,length=500)
    Figure = image(width,height,imageₜ,show_axis=false,color=:black,transparency=true)

    lines!( hcat(boundary.(angle,Ref(parameters))...)', linewidth=3, color=:gold )
    RecordEvents( hbox( Figure, vbox(t), parent = Scene(resolution=(800, 800))), "output" )
end

dataIdx = findall(binaryImages[:,:,1,26])[1:10:end]
loss(dataIdx,parameters)






# img = load("data/E2/M_allZ.tif")

# img = load("data/E2/U.tif")
# img[:,:,50]

# img = load("data/E5/Mb.tif")
# img[:,:,end]

# img = load("data/E11/iM.tif")
# img[:,:,end]