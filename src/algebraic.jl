export WH, coredisc, conductor, pell, pellreg, quadclassunit


@doc raw"""
    wh( p::Vector{Integer}, d::Integer [, T::Type = BigFloat])
    wh( m::Integer, n::Integer, d::Integer [, T::Type = BigFloat])

Weyl-Heisenberg displacement operators.

----------
    wh( p::Vector{<:Integer}, v::Vector)
    wh( m::Integer, n::Integer, v::Vector)

Action of `wh(p,length(v))` or `wh(m,n,length(v))` onto the vector `v` without explicitly forming the matrix. 
The output type is `typeof(v)`. 
"""
wh( m::Integer, n::Integer, d::Integer, T::Type=BigFloat) = 
    [ (-e(T(1)/(2*d)))^(m*n) * (rem(j-k-m,d) == 0) * e((k-1)*T(n)/d) for j=1:d, k=1:d]

wh( p::Vector{<:Integer}, d::Integer, T::Type=BigFloat) = wh(p[1],p[2],d,T)

wh( m::Integer, n::Integer, v::Vector) = 
    [ (-e( eltype(v)(1)/(2*length(v))) )^((2k-m)*n) * v[mod(k-m,length(v))+1] for k=0:length(v)-1]

wh( p::Vector{<:Integer}, v::Vector) = wh(p[1],p[2],v)



@doc """
    coredisc(D)
    coredisc(Q::QuadBin)

Outputs a tuple `(Δ,f)` of the fundamental discriminant and conductor of `D`, or 
of a binary quadratic form with discriminant D.
"""
coredisc(D) = ( fundamental_discriminant(D), conductor(D) )
coredisc(Q::QuadBin) = coredisc(discriminant(Q))



@doc raw"""
    conductor(Q::QuadBin)

The conductor of discriminant of the integral binary quadratic form Q.
"""
conductor(Q::QuadBin) = conductor(discriminant(Q))



@doc raw"""
    pell(D)
    pell(Q::QuadBin)
    pell(Zω::AbsSimpleNumFieldOrder)

Finds the least unit > 1 having norm 1 for the discriminant `D`, 
for the binary quadratic form with discriminant `D`, 
or the quadratic order Zω.
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
    if n == -1 u=u^2 end
    u
end
pell(Q::QuadBin) = pell(discriminant(Q))


@doc raw"""
    pellreg(D)
    pellreg(Q::QuadBin)
    pellreg(Zω::AbsSimpleNumFieldOrder)

Finds the least unit > 1 having norm 1 for the discriminant `D`, 
for the binary quadratic form with discriminant `D`, 
or the quadratic order Zω. 
The output is a tuple with the unit and the log of that unit (which is either the regulator, or twice the regulator) as a BigFloat.
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
@doc raw"""
    ghostclassfield( K::AbsSimpleNumField, q)

Compute the ring class field for the order q*Z(K) where Z(K) is the maximal order in K. 
The output is an `AbsSimpleNumField`, and so is an absolute field rather than a relative extension. 
"""
function ghostclassfield( K::AbsSimpleNumField, q::Integer)
    # @assert degree(K) == 2
    rcf = Hecke.ring_class_field( Order(K, q*basis(maximal_order(K))))
    simplify(absolute_simple_field(number_field(rcf))[1])[1]
end

ghostclassfield( F::AdmissibleTuple ) = ( (@isinit F.H) ? F.H : (@init! F.H = ghostclassfield(F.K,F.q)); F.H)

# Note: no unit tests yet for this function
@doc raw"""
    signswitch( H::AbsSimpleNumField, D::Integer)

If `H` is the (absolute) ring class field for a ghost with fundamental discriminant `D` with some conductor, then this finds a sign-switching Galois automorphism `g` on `H`, that is `g(√D) = -√D`. 
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

Returns class group and unit group data for the quadratic order with disciminant D. 
Let ``\omega = \bigl(\Delta\bmod 4 + \sqrt{\Delta}\bigr)/2``, so that a ``\mathbb{Z}`` basis is ``\mathbb{Z}+\mathbb{Z}[\omega]``.  
Then the output is a tuple `(c,e,b,u)`, where 
 - `c` is the class number, 
 - `e` is an integer vector for the cycle stucture of the class group, ``\mathbb{Z}^{e_1} + \ldots + \mathbb{Z}^{e_r}``.
 - `b` is a vector of binary quadratic forms that generate the corresponding factor in the class group,
 - `u` is a totally positive fundamental unit with norm 1, written as `[x,y]` in the above basis.
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



@doc raw"""
    sicnum(d)

Number of WH SICs in dimension d, modulo EC orbits.
"""
function sicnum(d)
    # fundamental discriminant and conductor
    D = (d+1)*(d-3)
    Δ,f0 = coredisc(D)
    F = sort(divisors(f0))
    s = 0
    for f in F
        D = f^2*Δ
        K, a = quadratic_field(D)
        ω = (D%4 + a)//2
        # generate the quadratic order from the standard basis
        Zω = Order( [one(ω), ω])
    
        # compute the class group and generator map
        cg, cm = picard_group(Zω)
    
        # the class number
        s += order(cg)
    end
    s
end



@doc raw"""
    ghostbasis(d)

Compute a basis for the allowed quadratic forms in dimension d. 

For each conductor, 
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


@doc raw"""
    ghostelements(d)

Compute a class representative for each ghost class in dimension d.
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



@doc raw"""
    hj(Q,n=0)

Takes a quadratic form Q and applies the Hirzebruch-Jung map to convert a Euclidean reduced form to an HJ reduced form at level `n`.
"""
hj(Q::QuadBin,n=0) = QuadBin([-1 1; -n-1 n]'*qmat(Q)*[-1 1; -n-1 n]) 
# inverse of [n -1; n+1 -1].


function reduced_HJ_orbit(q::QuadBin{ZZRingElem})
    rootD = sqrt(1.0*discriminant(q))
    Q = reduction(q)
    V = cycle(Q)
    n0 = map( x -> ceil(Int,(2.0*abs(x.a))/(-1.0x.b+rootD))-2, V )
    R = typeof(V[1])[]
    for k = 1:length(n0), n=0:n0[k]
        push!( R, hj(V[k],n) )
    end
    R    
end


function minimal_HJ_stabilizer(V::Vector{QuadBin{ZZRingElem}},d)
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
    # Q, L, n, A, W = Q[k], L[k], n[k], A[k], W[k]
    Q = Q[k]    
    sort!(Q; lt = (x,y) -> quadcompare_sum_then_max(x,y))
    Q[1]
end


# use this to sort a list of binary quadratic forms
function quadcompare_sum_then_max(P::QuadBin,Q::QuadBin)
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
