export guess_int_null_vec

@doc """
    guess_int_null_vec( x::Vector{BigFloat})

Use LLL to find an integer relation among the elements of `x`.

# Examples

```jldoctest
julia> guess_int_null_vec(BigFloat.([1; sin(big(pi)/8)^2; sin(big(pi)/4)]))
3-element Vector{BigInt}:
 -1
  2
  1
```
"""
function guess_int_null_vec(x::Vector{BigFloat})
    prec = precision(x[1])
    t = ZZ(2)^prec
    n = length(x)
    v = ZZ.(round.(t .* x))
    L = matrix_space(ZZ, n, n + 1)
    T = L([ZZ(j == k) + (k == n + 1) * v[j] for j = 1:n, k = 1:n+1])

    return BigInt.(lll(T)[1, :])[1:n]
end
