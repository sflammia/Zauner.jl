export wh, coredisc, conductor, pell, pellreg, ghostclassfield, signswitch, quadclassunit, numsics, ghostbasis, ghostelements

@doc """
    wh( p::Vector{<:Integer}, d::Integer [, T::Type = BigFloat])
    wh( m::Integer, n::Integer, d::Integer [, T::Type = BigFloat])

    wh( p::Vector{<:Integer}, v::Vector)
    wh( m::Integer, n::Integer, v::Vector)
    wh( n::Integer, v::Vector)

Weyl-Heisenberg displacement operators.
We have `wh(p,q,d) == wh([p,q],d)` acts on the standard basis as ``|k\\rangle \\to v_d^{p q} ω_d^{q}|k+p\\rangle``, where arithmetic inside the ket is modulo ``d`` and where ``ω_d = v_d^2`` and ``v_d = -\\mathrm{e}^{i π/d}``.
These forms explicitly construct the matrix that acts this way where `0` is the first element of the basis.

The forms `wh(p,v)` or `wh(m,n,v)` give the action onto the vector `v` without explicitly forming the matrix.
This is much faster when working with high dimensions or high precision.
The form `wh(n,v)` with `d = length(v)` simply converts the integer `n` into the base-`d` expansion ``n = n_1 d + n_0`` and calls `wh([n1,n0],v)`.
Tuples can also be used in place of vectors in these function calls.

# Examples
```jldoctest
julia> wh(2,3,4)
4×4 Matrix{Complex{BigFloat}}:
 0.0-0.0im  -0.0+0.0im  0.0+1.0im  0.0-0.0im
 0.0-0.0im  -0.0+0.0im  0.0+0.0im  1.0-0.0im
 0.0-1.0im  -0.0+0.0im  0.0+0.0im  0.0-0.0im
 0.0-0.0im  -1.0+0.0im  0.0+0.0im  0.0-0.0im
```
```jldoctest
julia> v = [1; 0; 0; 0];

julia> wh(1,2,v) ≈ [0.0; 1.0im; 0.0; 0.0]
true
```
"""
wh( m::Integer, n::Integer, d::Integer, T::Type=BigFloat) =
    [ (-e(T(1)/(2*d)))^(m*n) * (rem(j-k-m,d) == 0) * e((k-1)*T(n)/d) for j=1:d, k=1:d]

wh( p::Vector{<:Integer}, d::Integer, T::Type=BigFloat) = wh(p[1],p[2],d,T)

wh( p::Tuple{Integer,Integer}, d::Integer, T::Type=BigFloat) = wh(p[1],p[2],d,T)

wh( m::Integer, n::Integer, v::Vector) =
    [ (-e( eltype(v)(1)/(2*length(v))) )^((2k-m)*n) * v[mod(k-m,length(v))+1] for k=0:length(v)-1]

wh( p::Vector{<:Integer}, v::Vector) = wh(p[1],p[2],v)

wh( p::Tuple{Integer,Integer}, v::Vector) = wh(p[1],p[2],v)

wh( n::Integer, v::Vector) = wh( divrem( n, length(v)), v)


@doc """
    coredisc(D)
    coredisc(Q::QuadBin)

Outputs a tuple `(Δ,f)` of the fundamental discriminant `Δ` and conductor `f` of `D`, or
of a binary quadratic form with discriminant `D`.

# Examples

Here are two examples with fundamental discriminant 5 and conductor 3.
```jldoctest
julia> coredisc(45) # 45 = 3^2 * 5
(5, 3)
```
```jldoctest
julia> coredisc( QuadBin(1,-7,1) ) # -4*det( [1 -7//2; -7//2 1]) == 3^2 * 5
(5, 3)
```

See also `conductor`.
"""
coredisc(D) = ( fundamental_discriminant(D), conductor(D) )
coredisc(Q::QuadBin) = coredisc(discriminant(Q))



@doc """
    conductor(Q::QuadBin)

The conductor of discriminant of the integral binary quadratic form Q.
"""
conductor(Q::QuadBin) = conductor(discriminant(Q))



@doc """
    pell(D)
    pell(Q::QuadBin)
    pell(Zω::AbsSimpleNumFieldOrder)

Finds the least unit > 1 having norm 1 for the discriminant `D`,
for the binary quadratic form with discriminant `D`,
or the quadratic order ℤ(ω) where ``ω = \\frac{1}{2}\\bigl( (D \\mod 4) +  \\sqrt{D}\\bigr)``.

# Example

The fundamental unit for ``ℚ(\\sqrt{5})`` is ``(\\sqrt{5}+1)/2``, but this has norm -1, so we return the square of this.

```jldoctest
julia> pell(5)
1//2*sqrt(5) + 3//2
```
"""
function pell(D)
    @req D>0 && is_discriminant(D) "D must be a positive discriminant."
    Δ, f = coredisc(D)
    K, a = quadratic_field(Δ)
    ω = (D%4 + f*a)//2
    # generate the quadratic order from the standard basis
    Zω = Order([one(ω), ω])
    pell(Zω)
end
function pell(Zω::AbsSimpleNumFieldOrder)
    @req degree(Zω)==2 "Order must have degree 2."

    Ug, fu = unit_group(Zω)
    u = fu(Ug[2])
    x, y = Int.(sign.(coordinates(one(Zω.nf)*u))) # coordinates wrt K basis, not Zω.
    n = norm(u)
    u = ( x*y > 0 ? x*u : x*n*u^(-1)) # choose totally positive unit
    if n == -1 u=u^2 end # if the fundamental unit has norm -1, we square it
    u
end
pell(Q::QuadBin) = pell(discriminant(Q))


@doc """
    pellreg(D)
    pellreg(Q::QuadBin)
    pellreg(Zω::AbsSimpleNumFieldOrder)

Finds the least unit > 1 having norm 1 for the discriminant `D`,
for the binary quadratic form with discriminant `D`,
or the quadratic order ℤ(ω) where ``ω = \\frac{1}{2}\\bigl( (D \\mod 4) +  \\sqrt{D}\\bigr)``.
The output is a tuple with the unit and a `BigFloat` of the log of that unit (which is just the regulator of ℤ(ω), or twice that if the fundamental unit has norm -1).

!!! warning
    A limitation of the current implementation is that the regulator is only correct to about 128 bits, even though it returns a number to the current precision of `BigFloat`.

See also `pell`.
"""
function pellreg(D)
    @req D>0 && is_discriminant(D) "D must be a positive discriminant."
    Δ, f = coredisc(D)
    K, a = quadratic_field(Δ)
    ω = (D%4 + f*a)//2
    # generate the quadratic order from the standard basis
    Zω = Order([one(ω), ω])
    pellreg(Zω)
end
function pellreg(Zω::AbsSimpleNumFieldOrder)
    @req degree(Zω)==2 "Order must have degree 2."

    Ug, fu = unit_group(Zω)
    u = fu(Ug[2])
    x, y = Int.(sign.(coordinates(one(Zω.nf)*u))) # coordinates wrt K basis, not Zω.
    n = norm(u)
    u = ( x*y > 0 ? x*u : x*n*u^(-1)) # choose totally positive unit
    if n == -1 u=u^2 end
    u, BigFloat(regulator(Zω)*(3-n)/2) # only accurate to about 127 bits, in fact.
end
pellreg(Q::QuadBin) = pellreg(discriminant(Q))


# Note: no unit tests yet for this function
@doc """
    ghostclassfield( K::AbsSimpleNumField, q)
    ghostclassfield( F::AdmissibleTuple)

Compute the ring class field for the order q*ℤ(K) where ℤ(K) is the maximal order in K.
The output is an LLL-reduced `AbsSimpleNumField`, and so is an absolute field rather than a relative extension.

Given an `AdmissibleTuple` `F` it initializes the field `F.H` to be this ring class field.
"""
function ghostclassfield( K::AbsSimpleNumField, q::Integer)
    # @assert degree(K) == 2
    rcf = Hecke.ring_class_field( Order(K, q*basis(maximal_order(K))))
    simplify(absolute_simple_field(number_field(rcf))[1])[1]
end

ghostclassfield( F::AdmissibleTuple ) = ( (@isinit F.H) ? F.H : (@init! F.H = ghostclassfield(F.K,F.q)); F.H)

# Note: no unit tests yet for this function
@doc """
    signswitch( H::AbsSimpleNumField, D::Integer)
    signswitch( F::AdmissibleTuple)

If `H` is the (absolute) ring class field for a ghost with fundamental discriminant `D` with some conductor, then this finds a sign-switching Galois automorphism `g` on `H`, that is ``g\\bigl(\\sqrt{D}\\bigr) = -\\sqrt{D}``.

If `F` is an `AdmissibleTuple` then is initializes the field `F.g`.
This requires that the field `F.H` has already been initialized with `ghostclassfield`.
"""
function signswitch( H::AbsSimpleNumField, D::Integer)
    _, x = H["x"]
    f = collect(keys((factor(x^2 - D).fac)))
    @assert length(f) == 2 # (x^2 - D) should factor completely in H
    r = Hecke.evaluate(f[1],0) # a root of x^2 – D
    @assert r^2 == D # sanity check
    auts = automorphism_list(H)
    return auts[findfirst([ s(r) == -r for s in auts])]
end

signswitch( F::AdmissibleTuple) = ( (@isinit F.g) ? F.g : (@init! F.g = signswitch(F.H,F.D)); F.g)


@doc raw"""
    quadclassunit(D)

Returns class group and unit group data for the quadratic order with disciminant `D`.
Let ``\omega = \bigl(\Delta\bmod 4 + \sqrt{\Delta}\bigr)/2``, so that a ``\mathbb{Z}``-basis is ``\mathbb{Z}+\omega\mathbb{Z}``.
Then the output is a tuple `(h,c,b,u)`, where
 - `h` is the class number,
 - `c` is an integer vector for the cycle stucture of the class group, ``\mathbb{Z}/{n_1} + \ldots + \mathbb{Z}/{n_r}``.
 - `b` is a vector of binary quadratic forms that generate the corresponding factor in the class group,
 - `u` is the totally positive fundamental unit with norm 1, written as `[x,y]` in the above basis.
"""
function quadclassunit(D)
    # fundamental discriminant and conductor
    Δ,f = coredisc(D)
    K, a = quadratic_field(Δ)
    ω = (D%4 + f*a)//2
    # generate the quadratic order from the standard basis
    Zω = Order([one(ω), ω])
    # fundamental totally positive unit of norm 1
    u = pell(Zω)

    # compute the class group and generator map
    cg, cm = picard_group(Zω)

    # the class number
    gbcn = order(cg)

    # the generators, cycle decomposition, and (reduced) quadratic form generators
    gcg = gens(cg)
    gbcyc = (gbcn == 1 ? [ZZ(1)] : diagonal(rels(cg)) )
    gbgens = (length(gcg) == 0 ? [quadbinid(D)] : [QuadBin(cm(g)) for g in gcg])
    gbgens = reduction.(gbgens)

    # compute the unit group and the expansion in the standard Z basis for Zω
    y = ZZ(f*trace(a*u)//D)
    x = ZZ((trace(u) - y*(D%4))//2)

    (gbcn,gbcyc,gbgens,u,[x,y])
end



@doc """
    numsics(d)

Number of WH 1-SICs in dimension d, modulo EC orbits.

# Example

Here are the number of inequivalent 1-SICs for dimensions 4-15:
```jldoctest
julia> foreach(d -> println( (d, numsics(d)) ), 4:15)
(4, 1)
(5, 1)
(6, 1)
(7, 2)
(8, 2)
(9, 2)
(10, 1)
(11, 3)
(12, 2)
(13, 2)
(14, 2)
(15, 4)
```
"""
function numsics(d)
    # fundamental discriminant and conductor
    D = (d+1)*(d-3)
    Δ, f0 = coredisc(D)
    F = sort(divisors(f0))
    s = 0
    for f in F
        D = f^2*Δ
        _, a = quadratic_field(D)
        ω = (D%4 + a)//2
        # generate the quadratic order from the standard basis
        Zω = Order( [one(ω), ω])

        # compute the class group
        cg, _ = picard_group(Zω)

        # the class number
        s += order(cg)
    end
    s
end



@doc """
    ghostbasis(d)

Compute a generating set (under Gauss composition) for the allowed quadratic forms in dimension d.

This returns a dictionary where each key is a divisor `f` of the conductor of `D = (d+1)*(d-3)` and each value is the tuple `(c,Q)` where `c[k]` is the order of `Q[k]` under Gauss composition.

# Examples

```jldoctest
julia> ghostbasis(15)
Dict{ZZRingElem, Tuple{Vector{ZZRingElem}, Vector{QuadBin}}} with 3 entries:
  4 => ([2], [Binary quadratic form over ZZ: 3*x^2 + 12*x*y - 4*y^2])
  2 => ([1], [Binary quadratic form over ZZ: x^2 + 6*x*y - 3*y^2])
  1 => ([1], [Binary quadratic form over ZZ: x^2 + 2*x*y - 2*y^2])
```
In this example, `Q = QuadBin(3,12,-4)` has order 2, so we have `Q^2 == QuadBin(1,12,-12)` since this is the principle form with the same discriminant.
"""
function ghostbasis(d)
    D = (d+1)*(d-3)
    Δ,f0 = coredisc(D)
    F = sort(divisors(f0))
    B = Dict{ZZRingElem,Tuple{Vector{ZZRingElem},Vector{QuadBin}}}()
    for f in F
        D = f^2*Δ
        K, a = quadratic_field(D)
        ω = (D%4 + a)//2
        # generate the quadratic order from the standard basis
        Zω = Order( [one(ω), ω])
        # compute the class group and generator map
        cg, cm = picard_group(Zω)

        # the generators, cycle decomposition, and (reduced) quadratic form generators
        gcg = gens(cg)
        gbcyc = (order(cg) == 1 ? [ZZ(1)] : diagonal(rels(cg)) )
        gbgens = (length(gcg) == 0 ? [quadbinid(D)] : [QuadBin(cm(g)) for g in gcg])

        gbgens = reduction.(gbgens)
        B[f] = (gbcyc, gbgens)
    end
    B
end


@doc """
    ghostelements(d)

Compute a representative quadratic form for each ghost class in dimension `d`.
There are a total of `numsics(d)` such forms.
These forms are merely Euclidean reduced, not Hirzebruch-Jung reduced.
"""
function ghostelements(d)
    D = (d+1)*(d-3)
    Δ,f0 = coredisc(D)
    F = sort(divisors(f0))
    B = QuadBin{ZZRingElem}[]
    for f in F
        D = f^2*Δ
        K, a = quadratic_field(D)
        ω = (D%4 + a)//2
        # generate the quadratic order from the standard basis
        Zω = Order( [one(ω), ω])
        # compute the class group and generator map
        cg, cm = picard_group(Zω)
        cgcn = order(cg)

        # the generators, cycle decomposition, and (reduced) quadratic form generators
        gcg = gens(cg)
        gbcyc = (cgcn == 1 ? [ZZ(1)] : diagonal(rels(cg)) )
        gbgens = (length(gcg) == 0 ? [quadbinid(D)] : [QuadBin(cm(g)) for g in gcg])

        gbgens = reduction.(gbgens)

        for k = 0:cgcn-1
            a = radix(k,gbcyc)
            push!( B, prod(gbgens .^ a) )
        end
    end
    B
end
