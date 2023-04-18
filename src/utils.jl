export radix

@doc raw"""
    radix(n,r)

The `length(r)` least significant digits of the integer `n` in mixed radix `r = [r1,r2,...,rk]`. 
# Examples
```jldoctest
julia> radix(5,[2 2 2 2 2])
5-element Vector{Int64}:
 0
 0
 1
 0
 1
julia> radix(86400,[7 24 60 60])
4-element Vector{Int64}:
 1
 0
 0
 0
```
"""
radix(n,r) = [div(n,prod(r[k+1:end])) % r[k] for k=1:length(r)]

