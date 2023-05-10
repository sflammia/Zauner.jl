export WH, coredisc, conductor, pell, towerh, quadclassunit


@doc raw"""
  WH(p,d)
  WH(m,n,d)

Weyl-Heisenberg displacement operators.
"""
WH(m,n,d,T::Type=BigFloat) = [ (-e(T(1)/(2*d)))^(m*n) * (rem(j-k-m,d) == 0) * e((k-1)*T(n)/d) for j=1:d, k=1:d]
WH(p,d,T::Type=BigFloat) = WH(p[1],p[2],d,T)


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
  pell(Zω::NfOrd)

Finds a totally positive fundamental unit of norm 1 for the discriminant `D`, 
for the binary quadratic form with some discriminant `D`, 
or the quadratic order Zω.
"""
function pell(D)
    @assert D>0 && is_discriminant(D) error("D must be a positive discriminant.")
    Δ, f = coredisc(D)
    K, a = quadratic_field(Δ)
    ω = (D%4 + f*a)//2
    # generate the quadratic order from the standard basis
    Zω = Order([one(ω), ω])
    pell(Zω)
end
function pell(Zω::NfOrd)
    @assert degree(Zω)==2 # only do this for quadratic orders
    
    Ug, fu = unit_group(Zω)
    u = fu(Ug[2])
    x, y = Int.(sign.(coordinates(one(Zω.nf)*u))) # coordinates wrt K basis, not Zω.
    u = ( x*y > 0 ? x*u : x*norm(u)*u^(-1))
    if norm(u) == -1 u=u^2 end
    u
end
pell(Q::QuadBin) = pell(discriminant(Q))



@doc raw"""
  towerh(d)
  towerh(d, u::NfOrdElem)


Tower height of a unit `u` for dimension `d`. 
If no unit is given, compute the tower height of dimension `d` for the fundamental discriminant by computing a fundamental unit.
"""
function towerh(u::NfOrdElem,d)
    em = complex_embedding(u.parent.nf, sqrt(Int(discriminant(u.parent.nf)/4)) )
    reg = Float64(log(abs(em))(u, 128))
    round(Int,acosh((d-1)/2)/reg)
end
function towerh(d)
    D = (d+1)*(d-3)
    Δ = fundamental_discriminant(D)
    u = pell(Δ)

    em = complex_embedding(u.parent.nf, sqrt(Float64(Δ)) )
    reg = Float64(log(abs(em))(u, 128))
    round(Int,acosh((d-1)/2)/reg)
end




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

For each conductor, the output adheres to Smith normal form, in the sense that it outputs 
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
    n0 = ceil(Int,1/maximum(map( x -> (-1.0x.b+rootD)/(2*abs(x.a)), V ))-1)
    R = hj.(V)
    for n=1:n0
        append!(R,hj.(V,n))
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
