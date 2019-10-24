"""
     floatvecs(seq_of_coords, T) where T is Float64 (default) or Float32
    
```    
vec_of_coords = [(1, 2), (2, 4)]
floatvecs(vec_of_coords) = [[1.0, 2.0], [2.0, 4.0]]

tup_of_coords = ((1, 2), (2, 4))
floatvecs(tup_of_coords, Float32) = ([1.0f0, 2.0f0], [2.0f0, 4.0f0])
```
"""
floatvecs(seq_of_coords, ::Type{T}=Float64) where {T} = map(x->[T.(x)...,], seq_of_coords)
