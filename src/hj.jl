export reduced_hj_orbit, minimal_hj_stabilizer, form_sign_fix

@doc """
    _hj_step(Q,n=0)

Takes a quadratic form Q and applies the Hirzebruch-Jung map to convert a Euclidean reduced form to an HJ reduced form at level `n`.
"""
_hj_step(Q::QuadBin, n=0) = QuadBin([-1 1; -n-1 n]' * qmat(Q) * [-1 1; -n-1 n])
# inverse of [n -1; n+1 -1].


@doc """
    reduced_hj_orbit(q::QuadBin{ZZRingElem})

Returns the reduced orbit under HJ reduction.

# Examples

```jldoctest
julia> reduced_hj_orbit(binary_quadratic_form(1,-4,1))
3-element Vector{QuadBin{ZZRingElem}}:
 Binary quadratic form over ZZ: x^2 - 4*x*y + y^2
 Binary quadratic form over ZZ: 3*x^2 - 6*x*y + 2*y^2
 Binary quadratic form over ZZ: 2*x^2 - 6*x*y + 3*y^2
```
"""
function reduced_hj_orbit(q::QuadBin{ZZRingElem})
    @assert discriminant(q) > 0
    rootD = sqrt(1.0 * discriminant(q))
    # Euclidean reduce the form q and ensure q.a > 0
    Q = form_sign_fix(reduction(q))
    # find the Euclidean reduced cycle with a > 0.
    V = cycle(Q)
    Vroots = map(x -> (rootD - 1.0 * x.b) / (2.0 * x.a), V)
    n0 = maximum(ceil.(Int, 1.0 ./ Vroots) .- 1)
    R = QuadBin{ZZRingElem}[]
    # Note: this could be faster by dropping elements of Vroots at the first n for which they fail.
    for n = 0:n0-1
        for k = 1:length(V)
            if Vroots[k] < 1 / (n + 1)
                push!(R, _hj_step(V[k], n))
            end
        end
    end
    R
end



@doc """
    form_sign_fix(q::QuadBin{ZZRingElem})

Take a form q = ⟨a,b,c⟩ and if a < 0 return the form r = ⟨-a,b,-c⟩.
Note that r is Euclidean reduced iff q is Euclidean reduced.
"""
function form_sign_fix(q::QuadBin{ZZRingElem})
    (q.a < 0 ? binary_quadratic_form(-q.a, q.b, -q.c) : q)
end



@doc """
    minimal_hj_stabilizer(V::Vector{QuadBin{ZZRingElem}},d)

Returns a minimal form `Q` among those in `V`.
Specifically, if `V` is a reduced HJ orbit then the output `Q` minimizes the length of the HJ expansion of the stabilizer of `Q` modulo `d`.
The function assumes that Q.a > 0 for each form in V.
To break ties on length, the algorithm further chooses a minimal form by minimizing according first to the sum of the absolute coefficients, then by the maximum absolute coefficient.
See the internal function `_quadcompare_sum_then_max` for details of the comparison function used for sorting.
"""
function minimal_hj_stabilizer(V::Vector{QuadBin{ZZRingElem}}, d)
    Q = deepcopy(V)
    L = stabilizer.(Q)
    n = sl2zorder.(L, d)
    A = L .^ n
    W = psl2word.(A)
    l = minimum(length.(W))
    k = findall(x -> length(x) == l, W)
    # these are the HJ reduced forms of shortest word length
    Q = Q[k]
    sort!(Q; lt=(x, y) -> _quadcompare_sum_then_max(x, y))
    Q[1]
end


@doc """
    _quadcompare_sum_then_max(P::QuadBin,Q::QuadBin)

Compare two forms `P` and `Q` and return `true` or `false` according to the following comparison.
Let `p = abs.([P.a,P.b,P.c])` and `q = abs.([Q.a,Q.b,Q.c])`.
Return `true` if `sum(p) < sum(q)`, and `false` if `sum(p) > sum(q)`.
If `sum(p) == sum(q)`, then return the value `maximum(p) < maximum(q)`.
"""
function _quadcompare_sum_then_max(P::QuadBin, Q::QuadBin)
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
