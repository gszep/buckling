using HDF5
function File(path::String; t=:, x=:, y=:, z=:)
    filename,ext = splitext(path)

    if occursin("h5",ext)
        A = h5read(path,"data",(:,x,y,z,t))
        @assert((size(A,1)==1)&(ndims(A)==5),"allowed formats: single channel TIF (x,y,t) or HDF5 (x,y,z,t)")

        return A[1,:,:,:,:]

    elseif occursin("tif",ext)
        A = Gray.(load(path))
        @assert((ndims(A)==3),"allowed formats: single channel TIF (x,y,t) or HDF5 (x,y,z,t)")
        
        newaxis = [CartesianIndex()]
        return rawview(channelview(A))[:,:,newaxis,:]
    else 
        throw("allowed formats: single channel TIF (x,y,t) or HDF5 (x,y,z,t)")
    end
end