## hints
    
- If your points are disaggregated (e.g. all the `xs` in vec_of_xs, all the `ys` in vec_of_ys)
    - aggregate them this way `points = collect(zip(xs, ys, zs))`

- Often, the `xs` respect `x[i-1] < x[i] < x[i+1]` or `x[i-1] > x[i] > x[i+1]`     
    - when the path is a closed curve, one triplet may follow    
    - `x[i-1] < x[i] > x[i+1]` or `x[i-1] > x[i] < x[i+1]`

