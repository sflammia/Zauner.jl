export reduced_hj_orbit, minimal_hj_stabilizer

@doc """
    _hj_step(Q,n=0)

Takes a quadratic form Q and applies the Hirzebruch-Jung map to convert a Euclidean reduced form to an HJ reduced form at level `n`.
"""
_hj_step(Q::QuadBin,n=0) = QuadBin([-1 1; -n-1 n]'*qmat(Q)*[-1 1; -n-1 n])
# inverse of [n -1; n+1 -1].

@doc """
    reduced_hj_orbit(q::QuadBin{ZZRingElem})

Returns the reduced orbit under HJ reduction.

# Examples

```jldoctest
julia> reduced_hj_orbit(QuadBin(1,-4,1))
3-element Vector{QuadBin{ZZRingElem}}:
 Binary quadratic form over ZZ: x^2 - 4*x*y + y^2
 Binary quadratic form over ZZ: 3*x^2 - 6*x*y + 2*y^2
 Binary quadratic form over ZZ: 2*x^2 - 6*x*y + 3*y^2
 ```
"""
function reduced_hj_orbit(q::QuadBin{ZZRingElem})
    @assert discriminant(q) > 0
    rootD = sqrt(1.0*discriminant(q))
    Q = reduction(q)
    V = cycle(Q)
    n0 = map( x -> ceil(Int,(2.0*abs(x.a))/(-1.0x.b+rootD))-2, V )
    R = typeof(V[1])[]
    for k = 1:length(n0), n=0:n0[k]
        push!( R, _hj_step(V[k],n) )
    end
    R
end

@doc """
    minimal_hj_stabilizer(V::Vector{QuadBin{ZZRingElem}},d)

Returns a minimal form `Q` among those in `V`.
Specifically, if `V` is a reduced HJ orbit then the output `Q` minimizes the length of the HJ expansion of the stabilizer of `Q` modulo `d`.
To break ties on length, the algorithm further chooses a minimal form by minimizing according first to the sum of the absolute coefficients, then by the maximum absolute coefficient.
See the internal function `_quadcompare_sum_then_max` for details of the comparison function used for sorting.
"""
function minimal_hj_stabilizer(V::Vector{QuadBin{ZZRingElem}},d)
    Q = deepcopy(V)
    # Flip to an equivalent form with Q.a > 0.
    for q in Q
        x = sign(q.a)
        q.a, q.c = x*q.a, x*q.c
    end
    L = stabilizer.(Q)
    n = sl2zorder.(L,d)
    A = L.^n
    W = psl2word.(A)
    l = minimum(length.(W))
    k = findall(x->length(x) == l,W)
    # these are the HJ reduced forms of shortest word length
    Q = Q[k]
    sort!(Q; lt = (x,y) -> _quadcompare_sum_then_max(x,y))
    Q[1]
end


@doc """
    _quadcompare_sum_then_max(P::QuadBin,Q::QuadBin)

Compare two forms `P` and `Q` and return `true` or `false` according to the following comparison.
Let `p = abs.([P.a,P.b,P.c])` and `q = abs.([Q.a,Q.b,Q.c])`.
Return `true` if `sum(p) < sum(q)`, and `false` if `sum(p) > sum(q)`.
If `sum(p) == sum(q)`, then return the value `maximum(p) < maximum(q)`.
"""
function _quadcompare_sum_then_max(P::QuadBin,Q::QuadBin)
    p = abs.([P.a, P.b, P.c])
    q = abs.([Q.a, Q.b, Q.c])
    diff = sum(p .- q)
    if diff < 0 # sum(p) < sum(q)
        return true
    elseif diff > 0
        return false
    else
        return maximum(p) < maximum(q)
    end
end
