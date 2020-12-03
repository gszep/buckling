using ReverseDiff: record_minus,TrackedArray
using StaticArrays: StaticArray

import Base: -
-( x::StaticArray, y::TrackedArray{V,D} ) where {V,D} = record_minus(x,y,D)
-( x::TrackedArray{V,D}, y::StaticArray ) where {V,D} = record_minus(x,y,D)