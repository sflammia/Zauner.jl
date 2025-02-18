import Base.==

export radix, hash, AdmissibleTuple, is_admissible, is_antiunitary, is_antiunitary_with_generator, has_fa_symmetry

# Define a hash for QuadBin type. Useful for a Dict{QuadBin}.
Base.hash(q::QuadBin{ZZRingElem}, h::UInt) = hash(q.a, hash(q.b, hash(q.c, h)))


@doc """
    radix(n,r)

The `length(r)` least significant digits of the integer `n` in mixed radix `r = [r1; r2; ...; rk]`.
# Examples
```jldoctest
julia> radix(5,[2; 2; 2; 2; 2])
5-element Vector{Int64}:
 0
 0
 1
 0
 1

julia> radix(86400,[7; 24; 60; 60])
4-element Vector{Int64}:
 1
 0
 0
 0
```
"""
radix(n, r) = [div(n, prod(r[k+1:end])) % r[k] for k = 1:length(r)]


# convert expansion n back from mixed radix r
# fromradix(n,r) = dot( cumprod(reverse(push!(r,1)))[end-1:-1:1], n)


@lazy struct AdmissibleTuple
    d::Integer                     # dimension
    r::Integer                     # rank
    n::Integer                     # (d^2-1)/(r*(d-r))
    D::Integer                     # fundamental discriminant of K
    f::Integer                     # sequence of conductors f_j
    q::Integer                     # conductor of Q
    j::Integer                     # grid vertical
    m::Integer                     # grid horozontal
    k::Integer                     # L^k = I mod d
    K::AbsSimpleNumField           # associated field K = ℚ(√n(n-4))
    a::AbsSimpleNumFieldElem       # associated field generator √D
    u::AbsSimpleNumFieldOrderElem  # totally positive fundamental unit > 1
    @lazy H::AbsSimpleNumField     # ring class field for q*Z(K)
    @lazy g::NumFieldHom           # Galois automorphism in H s.t. g(√D) = -√D
    h::Integer                     # class number, or degree of H/K
    G::Vector{Integer}             # orders of class group generators
    c::Vector{QuadBin}             # basis of generator Q's for the class group
    Q::QuadBin                     # integral binary quadratic form
    L::Matrix{ZZRingElem}          # S(Q) generator
    A::Matrix{ZZRingElem}          # S_d(Q) generator
    x::BigFloat                    # positive root of Q
    R::BigFloat                    # log of u
end


Base.show(io::IO, ::MIME"text/plain", F::AdmissibleTuple) = print(io, "AdmissibleTuple( d = $(F.d), ", (F.r == 1 ? "" : "r = $(F.r), "), "K = ℚ(√$(F.D)), q = $(F.q), Q = ⟨$(F.Q.a),$(F.Q.b),$(F.Q.c)⟩, h = $(F.h) )")
Base.show(io::IO, F::AdmissibleTuple) = print(io, "AdmissibleTuple( d = $(F.d), ", (F.r == 1 ? "" : "r = $(F.r), "), ", K = ℚ(√$(F.D)), q = $(F.q), Q = ⟨$(F.Q.a),$(F.Q.b),$(F.Q.c)⟩, h = $(F.h) )")


@doc """
    AdmissibleTuple( d::Integer [, Q::QuadBin])
    AdmissibleTuple( d::Integer, r::Integer [, Q::QuadBin])
    AdmissibleTuple( D::Integer, j::Integer, m::Integer [, Q::QuadBin])

    AdmissibleTuple( dQ::Tuple{Integer, QuadBin{ZZRingElem}})
    AdmissibleTuple( drQ::Tuple{Integer, Integer, QuadBin{ZZRingElem}})

Data type for the arithmetic data defining a set of ghost overlaps.

If only one integer is given, it is interpreted as a dimension and assumed that the rank is `r = 1`.
The syntax that takes three integers `(D,j,m)` requires that `D` is a *fundamental* discriminant for a real quadratic field.

In all three cases, if the optional argument `Q` is left unspecified, then it defaults to a principal form given by `Q = binary_quadratic_form( 1, 2-n, 1)`, where `n` is the integer `(d^2-1)/(r(d-r))`, since this has ``\\mathrm{disc}(Q) = n(n-4)``.

The final way to specify an input is as a tuple `(d,Q)` or `(d,r,Q)`.

The defined (and precomputed) fields in an `AdmissibleTuple` are given by:

```
    d ::Integer                     # dimension
    r ::Integer                     # rank
    n ::Integer                     # (d^2-1)/(r*(d-r))
    K ::AbsSimpleNumField           # associated field K = ℚ(√n(n-4))
    D ::Integer                     # fundamental discriminant of K
    f ::Integer                     # conductor, where n*(n-4) = D*f^2
    q ::Integer                     # conductor of Q
    j ::Integer                     # grid vertical position
    m ::Integer                     # grid horozontal position
    a ::AbsSimpleNumFieldElem       # associated field generator √D
    u ::AbsSimpleNumFieldOrderElem  # Zauner unit
    H ::AbsSimpleNumField           # (lazy) ring class field for q*Z(K)
    g ::NumFieldAut                 # (lazy) Galois automorphism in H s.t. g(√D) = -√D
    h ::Integer                     # class number, or degree of H/K
    G ::Vector{Integer}             # orders of class group generators
    c ::Vector{QuadBin}             # basis of generator Q's for the class group
    Q ::QuadBin                     # binary quadratic form
    L ::Matrix{ZZRingElem}          # S(Q) generator
    k ::Integer                     # L^k = I mod d
    A ::Matrix{ZZRingElem}          # S_d(Q) generator
    x ::BigFloat                    # positive root of Q
    R ::BigFloat                    # log of u
```

The fields marked '(lazy)', namely 'H' and 'g', are not initialized at first since their computation is substantially more expensive than the other fields.
Once they have been computed their values are memoized and do not have to be recomputed.

From the way that `H` is computed currently, `K` is not recognized as a subfield.
This could lead to some non-intuitive results.
For example, the quadratic generator `a` is an element of `K`, but *not* a recognized element of `H`.
If an explicit embedding is needed, it can be obtained using `is_subfield(K,H)`.

# Examples

```jldoctest
julia> AdmissibleTuple(5)
AdmissibleTuple( d = 5, K = ℚ(√12), q = 1, Q = ⟨1,-4,1⟩, h = 1 )
```
```jldoctest
julia> AdmissibleTuple(11,3)
AdmissibleTuple( d = 11, r = 3, K = ℚ(√5), q = 1, Q = ⟨1,-3,1⟩, h = 1 )
```
```jldoctest
julia> AdmissibleTuple(5,1,1)
AdmissibleTuple( d = 4, K = ℚ(√5), q = 1, Q = ⟨1,-3,1⟩, h = 1 )
```
```jldoctest
julia> AdmissibleTuple(11,binary_quadratic_form(3,-12,4))
AdmissibleTuple( d = 11, K = ℚ(√24), q = 2, Q = ⟨3,-12,4⟩, h = 2 )
```
"""
AdmissibleTuple(d::Integer) = AdmissibleTuple(d, 1)

function AdmissibleTuple(d::Integer, r::Integer)
    @req 0 < 2r < (d - 1) "r must satisfy 0 < 2r < d-1."
    n = Int((d^2 - 1) // (r * (d - r))) # throws an error if (d,r) is not admissible
    @req n > 4 "n must be > 4."
    D, f = Int.(coredisc(n * (n - 4)))
    K, a = quadratic_field(D)
    u, R = pellreg(D)
    t = sqrt(BigFloat(D)) * f / 2
    j = round(Int, asinh(t) / R)
    m = round(Int, asinh(t * r) / asinh(t))
    Q = binary_quadratic_form(1, 2 - n, 1) # disc(Q) = n*(n-4)
    q = Int(conductor(Q))
    # H = ghostclassfield(K,q)
    # g = signswitch(H,D)
    h, G, c, _, _ = quadclassunit(q^2 * D)
    h = Int(h)
    G = Int.(G)
    L = [ZZ(n)-2one(ZZ) -one(ZZ); one(ZZ) zero(ZZ)]
    k = Int(sl2zorder(L, d))
    A = L^k
    x = (-BigInt(Q.b) + sqrt(BigInt(discriminant(Q)))) / (2BigInt(Q.a))
    return AdmissibleTuple(d, r, n, D, f, q, j, m, k, K, a, u, uninit, uninit, h, G, c, Q, L, A, x, R)
end

AdmissibleTuple(d::Integer, Q::QuadBin) = AdmissibleTuple(d, 1, Q)

AdmissibleTuple(dQ::Tuple{Integer,QuadBin}) = AdmissibleTuple(dQ...)
AdmissibleTuple(drQ::Tuple{Integer,Integer,QuadBin}) = AdmissibleTuple(drQ...)

function AdmissibleTuple(d::Integer, r::Integer, Q::QuadBin)
    @req 0 < 2r < (d - 1) "r must satisfy 0 < 2r < d-1."
    n = Int((d^2 - 1) // (r * (d - r))) # throws an error if (d,r) is not admissible
    @req n > 4 "n must be > 4."
    D, f = Int.(coredisc(n * (n - 4)))
    DQ, q = coredisc(Q)
    q = Int(q)
    @req DQ == D "Fundamental discriminant of Q and K must match."
    @req f % q == 0 "Conductor q of Q must divide f_j."
    K, a = quadratic_field(D)
    u, R = pellreg(D)
    # H = ghostclassfield(K,q)
    # g = signswitch(H,D)
    h, G, c, _, _ = quadclassunit(q^2 * D)
    h = Int(h)
    G = Int.(G)
    t = sqrt(BigFloat(D)) * f / 2
    j = round(Int, asinh(t) / R)
    m = round(Int, asinh(t * r) / asinh(t))
    L = stabilizer(Q)
    k = Int(sl2zorder(L, d))
    A = L^k
    x = (-BigInt(Q.b) + sqrt(BigInt(discriminant(Q)))) / (2BigInt(Q.a))
    return AdmissibleTuple(d, r, n, D, f, q, j, m, k, K, a, u, uninit, uninit, h, G, c, Q, L, A, x, R)
end


function AdmissibleTuple(D::Integer, j::Integer, m::Integer)
    @req D > 1 && is_fundamental_discriminant(D) "D must be a positive fundamental discriminant."
    @req j > 0 && m > 0 "j and m must be > 0."
    K, a = quadratic_field(D)
    u, R = pellreg(D)
    f = Int(ZZ((u^j - u^(-j)) // a))
    r = Int(ZZ((u^(j * m) - u^(-j * m)) // (f * a)))
    d = Int(ZZ((u^(j * (m + 1)) - u^(-j * (m + 1))) // (f * a))) + r
    n = Int((d^2 - 1) // (r * (d - r))) # throws an error if (d,r) is not admissible
    Q = binary_quadratic_form(1, 2 - n, 1) # disc(Q) = n*(n-4)
    q = Int(conductor(Q))
    # H = ghostclassfield(K,q)
    # g = signswitch(H,D)
    h, G, c, _, _ = quadclassunit(q^2 * D)
    h = Int(h)
    G = Int.(G)
    L = [ZZ(n)-2one(ZZ) -one(ZZ); one(ZZ) zero(ZZ)]
    k = Int(sl2zorder(L, d))
    A = L^k
    x = (-BigInt(Q.b) + sqrt(BigInt(discriminant(Q)))) / (2BigInt(Q.a))
    return AdmissibleTuple(d, r, n, D, f, q, j, m, k, K, a, u, uninit, uninit, h, G, c, Q, L, A, x, R)
end


function AdmissibleTuple(D::Integer, j::Integer, m::Integer, Q::QuadBin)
    @req D > 1 && is_fundamental_discriminant(D) "D must be a positive fundamental discriminant."
    @req j > 0 && m > 0 "j and m must be > 0."
    K, a = quadratic_field(D)
    u, R = pellreg(D)
    f = Int((u^j - u^(-j)) // a)
    DQ, q = coredisc(Q)
    q = Int(q)
    @req DQ == D "Fundamental discriminant of Q and K must match."
    @req f % q == 0 "Conductor of Q must divide f_j."
    # H = ghostclassfield(K,q)
    # g = signswitch(H,D)
    h, G, c, _, _ = quadclassunit(q^2 * D)
    h = Int(h)
    G = Int.(G)
    r = Int((u^(j * m) - u^(-j * m)) // (f * a))
    d = Int((u^(j * (m + 1)) - u^(-j * (m + 1))) // (f * a)) - r
    n = Int((d^2 - 1) // (r * (d - r))) # throws an error if (d,r) is not admissible
    L = stabilizer(Q)
    k = Int(sl2zorder(L, d))
    A = L^k
    x = (-BigInt(Q.b) + sqrt(BigInt(discriminant(Q)))) / (2BigInt(Q.a))
    return AdmissibleTuple(d, r, n, D, f, q, j, m, k, K, a, u, uninit, uninit, h, G, c, Q, L, A, x, R)
end


@doc raw"""
    is_admissible( d::Integer [, Q::QuadBin])
    is_admissible( d::Integer, r::Integer [, Q::QuadBin])
    is_admissible( D::Integer, j::Integer, m::Integer [, Q::QuadBin])

Test whether the arguments form an admissible tuple.

# Examples

```jldoctest
julia> is_admissible(5,1)
true

julia> is_admissible(5,2)
false

julia> is_admissible(5,1,1)
true

julia> is_admissible(5,2,1,binary_quadratic_form(1,-7,1))
true

julia> is_admissible(11,1,binary_quadratic_form(3,-12,4))
true
```
"""
is_admissible(d::Integer) = is_admissible(d, 1)

function is_admissible(d::Integer, r::Integer)
    n = Int((d^2 - 1) // (r * (d - r)))
    0 < 2r < (d - 1) && n > 4
end

is_admissible(d::Integer, Q::QuadBin) = is_admissible(d, 1, Q)

function is_admissible(d::Integer, r::Integer, Q::QuadBin)
    n = Int((d^2 - 1) // (r * (d - r)))
    D, f = Int.(coredisc(n * (n - 4)))
    DQ, q = coredisc(discriminant(Q))

    is_admissible(d, r) && (D == DQ) && (f % q == 0)
end

function is_admissible(D::Integer, j::Integer, k::Integer)
    D > 1 && is_fundamental_discriminant(D) && j > 0 && k > 0
end

function is_admissible(D::Integer, j::Integer, k::Integer, Q::QuadBin)
    DQ, q = coredisc(discriminant(Q))
    K, a = quadratic_field(D)
    u = pell(D)
    f = Int(ZZ((u^j - u^(-j)) // a))

    is_admissible(D, j, k) && (D == DQ) && (f % q == 0)
end


# Equality testing for AdmissibleTuples
# need to add a tests for :H and :g with conditions for if they are defined.
function ==(F::AdmissibleTuple, G::AdmissibleTuple)
    test = F.Q.a == G.Q.a
    test &= F.Q.b == G.Q.b
    test &= F.Q.c == G.Q.c
    test &= coordinates(F.u) == coordinates(G.u)
    test &= F.x ≈ G.x
    test &= isapprox(F.R, G.R, atol=2^-120) # only good to about half BigFloat default prec.
    for p in (:d, :r, :n, :D, :f, :q, :j, :m, :k, :K, :a, :h, :L, :A)
        test &= getfield(F, p) == getfield(G, p)
    end
    return test
end




# Def 7.20, Lemma 7.21, and Def. 4.31 of main.tex
# Anti-unitary if and only if the following hold:
#     j_min must be odd;
#     d_min - 3 must be a square, k^2;
#     kf_Q must divide f_{j_min}.
@doc """
    is_antiunitary(F::AdmissibleTuple)

Test if the tuple `F` has antiunitary symmetry.
"""
function is_antiunitary(F::AdmissibleTuple)
    @assert F.r == 1 "Only rank = 1 is supported."

    fQ = conductor(F.Q)
    jmin = 0
    dmin = 0
    fmin = 0
    # calculate j_min, d_min, f_min
    for j = 1:F.j
        dmin = Integer(trace(F.u^j)) + 1
        _, fmin = coredisc((dmin + 1) * (dmin - 3))
        if fmin % fQ == 0
            jmin = j
            break
        end
    end

    if iseven(jmin)
        return false
    else
        test, k = is_square_with_sqrt(dmin - 3)
        return test && (fmin % (k * fQ) == 0)
    end

end

@doc """
    is_antiunitary_with_generator(F::AdmissibleTuple)

Test if the tuple `F` has antiunitary symmetry and return a tuple: if the first output is `true`, the second element of the tuple is a symmetry generator, and if `false`, the second output is the zero matrix.
"""
function is_antiunitary_with_generator(F::AdmissibleTuple)
    @assert F.r == 1 "Only rank = 1 is supported."

    fQ = conductor(F.Q)
    jmin = 0
    dmin = 0
    fmin = 0
    # calculate j_min, d_min, f_min
    for j = 1:F.j
        dmin = Integer(trace(F.u^j)) + 1
        _, fmin = coredisc((dmin + 1) * (dmin - 3))
        if fmin % fQ == 0
            jmin = j
            break
        end
    end

    if iseven(jmin)
        return false, ZZ.([0 0; 0 0])
    else
        test, k = is_square_with_sqrt(dmin - 3)
        if test && (fmin % (k * fQ) == 0)
            M = ZZ.(k // 2 * [1 0; 0 1] + [0 -1; 1 0] * qmat(F.Q) .* fmin // (k * fQ))
            return true, M
        else
            return false, ZZ.([0 0; 0 0])
        end
    end

end


@doc """
    has_fa_symmetry(F::AdmissibleTuple)

Test for ``F_a`` symmetry as opposed to Zauner (``F_z``) symmetry.
"""
function has_fa_symmetry(F::AdmissibleTuple)
    rem(F.d, 9) == 3 && rem(F.f ÷ F.q, 3) == 0
end
