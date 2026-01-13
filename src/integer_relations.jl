export guess_int_null_vec

@doc """
    guess_int_null_vec(x::Vector{BigFloat})
    guess_int_null_vec(x::Vector{BigFloat}, prec::Integer; int_tol::Integer=5, verbose::Bool=false)

Attempt to find a (nonzero) integer relation among the entries of `x` using LLL.

Returns a vector `a::Vector{BigInt}` such that `dot(a, x)` is numerically close to zero. If no
relation passes a verification check at the requested precision, returns `nothing`.

## Method

1. **Candidate via LLL (at `3*prec÷4` bits):** form `v = round(2^(3*prec÷4) .* x)`, build an
`n × (n+1)` lattice basis, and extract a candidate `a` from the first LLL-reduced row.
2. **Verification (at `prec` bits):** recompute `v = round(2^prec .* x)` and accept only if

`abs(dot(a, v)) ≤ int_tol * sum(abs, a)`.

`prec` is in **bits** and defaults to `precision(x[1])`. Smaller `int_tol` makes acceptance stricter.
If `verbose=true`, a warning is emitted when verification fails.

# Examples

```jldoctest
julia> guess_int_null_vec([one(BigFloat); sin(big(pi)/8)^2; sin(big(pi)/4)])
3-element Vector{BigInt}:
 -1
  2
  1
```
"""
function guess_int_null_vec(x::Vector{BigFloat})
    return guess_int_null_vec(x, precision(x[1]))
end

function guess_int_null_vec(x::Vector{BigFloat}, prec::Integer; int_tol::Integer=5, verbose::Bool=false)
    n = length(x)
    L = matrix_space(ZZ, n, n + 1)

    t = ZZ(2)^(3 * prec ÷ 4)
    v = ZZ.(round.(t .* x))
    M = L([ZZ(j == k) + (k == n + 1) * v[j] for j = 1:n, k = 1:n+1])
    a = BigInt.(lll!(M)[1, :])[1:n]

    # increase precision and verify
    s = sum(abs, a)
    t = ZZ(2)^(prec)
    v = ZZ.(round.(t .* x))

    if abs(dot(a, v)) > int_tol * s
        verbose && (@warn "No integer relation found")
        return nothing
    end

    return a
end
