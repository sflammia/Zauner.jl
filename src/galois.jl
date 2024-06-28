export galois_orbit, galois_normal_form, galois_elements, galois_order_orbit, centralizer_elements, stabilizer_elements

@doc """
    galois_orbit(F::AdmissibleTuple)
    galois_orbit(F::AdmissibleTuple, T::Matrix{ZZRingElem}, n::Integer)

Compute a maximal Galois orbit.
Output is a vector of `Vector{Integer}` elements lying on a maximal orbit for the tuple `F`.
The second form computes the Galois orbit with the centralizer subgroup generated by `T` of order `n` removed.
"""
function galois_orbit(F::AdmissibleTuple)
    d = Int(F.d)
    dd = 2^iseven(d)*d

    # compute the nontrivial centralizer orbits
    cents = centralizer_elements(F)
    inds = [ [x,y] for x = 0:dd-1, y = 0:dd-1 ][2:end] # nonzero orbits
    orbs = Vector{Vector{Integer}}[]
    while !isempty(inds)
        v = popfirst!(inds)
        vorb = unique(map( X -> mod.(X*v,dd) , cents))
        push!( orbs, vorb )
        setdiff!( inds, vorb)
    end

    # compute the elements of the stabilizer group
    # but don't keep the identity element
    stab = stabilizer_elements(F)[2:end]

    # quotient by the stabilizer
    for orb in orbs
        for p in orb
            # remove the equivalent points in the stabilizer orbit
            equiv = map( X -> mod.(X*p,dd) , stab)
            setdiff!( orb, equiv)
        end
    end

    # return a maximal orbit
    argmax( length, orbs)
end


# Compute the Galois orbit with the centralizer subgroup generated by T of order n removed
function galois_orbit(F::AdmissibleTuple, T::Matrix{ZZRingElem}, n::Integer)
    d = Int(F.d)
    dd = 2^iseven(d)*d

    # compute the nontrivial centralizer orbits
    cents = centralizer_elements(F)
    inds = [ [x,y] for x = 0:dd-1, y = 0:dd-1 ][2:end] # nonzero orbits
    orbs = Vector{Vector{Integer}}[]
    while !isempty(inds)
        v = popfirst!(inds)
        vorb = unique(map( X -> mod.(X*v,dd) , cents))
        push!( orbs, vorb )
        setdiff!( inds, vorb)
    end

    # compute the elements of the stabilizer group
    # but don't keep the identity element
    stab = stabilizer_elements[2:end]

    # these are the non-identity elements of E_j and Stab.
    Ejstab = [ T^k for k=1:n-1 ]
    for j=1:length(stab), k=0:n-1
        push!( Ejstab, mod.( stab[j]*T^k, dd) )
    end
    Ejstab = Matrix{Integer}.(Ejstab)

    # quotient by the stabilizer
    for orb in orbs
        for p in orb
            # remove the equivalent points in the stabilizer orbit
            equiv = map( X -> mod.(X*p,dd) , Ejstab)
            setdiff!( orb, equiv)
        end
    end

    # return a maximal orbit
    argmax( length, orbs)
end


@doc raw"""
    centralizer_elements(F::AdmissibleTuple)

Compute the elements of ``\mathrm{GL}(2,\mathbb{Z}/d')`` that are in the centralizer of the stability group for `F`.

# Examples

The centralizer for ``d=5`` has 24 elements.
```jldoctest
julia> F = AdmissibleTuple(5);

julia> cent = centralizer_elements(F);

julia> length(cent)
24
```
"""
function centralizer_elements(F::AdmissibleTuple)
    # No support presently for F_a orbits.
    has_fa_symmetry(F) && error("F_a orbits not supported yet")

    d = Int(F.d)
    dd = 2^iseven(d)*d #  d'
    Z = F.L
    eye = [ 1  0;  0  1]

    allM = map( x -> mod.(x[1]*eye+x[2]*Z,dd), CartesianIndices((0:dd-1,0:dd-1)))[:]
    elem = Matrix{Integer}[]
    while !isempty(allM)
        M = popfirst!(allM)
        # If p I + q Z is an element of GL(2,ℤ/dd)
        if gcd( dd, mod( M[1]*M[4]-M[2]*M[3], dd) ) == 1
            push!( elem, M)
        end
    end
    elem
end



@doc raw"""
    stabilizer_elements(F::AdmissibleTuple)

Compute the elements of ``\mathrm{GL}(2,\mathbb{Z}/d')`` that are in the complete stabilizer for `F`, including the extra antiunitary symmetry (if present).
The elements are ordered so that the identity element is first.

# Examples

```jldoctest
julia> F = AdmissibleTuple(7)
AdmissibleTuple( d = 7, K = ℚ(√8), q = 2, Q = ⟨1,-6,1⟩, h = 1 )

julia> stabilizer_elements(F)
3-element Vector{Matrix{ZZRingElem}}:
 [1 0; 0 1]
 [6 6; 1 0]
 [0 1; 6 6]
```
"""
function stabilizer_elements(F::AdmissibleTuple)
    d = Int(F.d)
    dd = 2^iseven(d)*d
    stab = typeof(F.L)[ ]
    test, L = is_antiunitary_with_generator(F)
    # if F is unitary, then F.L is the stabilizer generator.
    if !test
        L = F.L
    end
    stab = push!( stab, mod.(L,dd) )
    k=1
    while stab[k] != [1 0; 0 1]
        push!(stab, mod.(L*stab[k], dd))
        k += 1
    end
    return circshift(stab,1) # put the identity first.
end



@doc raw"""
    galois_normal_form(F::AdmissibleTuple)

Compute the normal form of the group ``\mathcal{M}/\mathcal{S}`` where ``\mathcal{M}`` is a maximal abelian subgroup of ``\mathrm{GL}(2,\mathbb{Z}/d')`` and ``\mathcal{S}`` is the stabilizer of `F`.
This is group is (conjecturally) isomorphic to the Galois group.
The output is a tuple `(g,n)` where `g` is a vector of generator matrices and `n` is a vector of the orders of those matrices in ``\mathcal{M}/\mathcal{S}``.
This corresponds to the canonical decomposition ``\mathbb{Z}/n_1 ⊕ \mathbb{Z}/n_2 ⊕ ... ⊕ \mathbb{Z}/n_k`` where the ``n_j`` are each prime powers.
The factors are canonical up to permutation (though the generators are not).

# Examples
```jldoctest
julia> F = AdmissibleTuple( 7, QuadBin(2,-4,1))
AdmissibleTuple( d = 7, K = ℚ(√8), q = 1, Q = ⟨2,-4,1⟩, h = 1 )

julia> g, n = galois_normal_form(F)
(Matrix{ZZRingElem}[[6 0; 0 6], [2 0; 0 2]], [2, 3])

julia> F = AdmissibleTuple(9)
AdmissibleTuple( d = 9, K = ℚ(√60), q = 1, Q = ⟨1,-8,1⟩, h = 2 )

julia> g, n = galois_normal_form(F)
(Matrix{ZZRingElem}[[8 5; 4 3], [8 0; 0 8], [7 0; 0 7]], [3, 2, 3])
```
"""
function galois_normal_form(F::AdmissibleTuple)
    d = Int(F.d)
    dd = 2^iseven(d)*d

    # G is the centralizer of H (the stabilizer)
    # compute the nontrivial centralizer orbits
    G = Matrix{ZZRingElem}.(centralizer_elements(F))
    sizeofG = length(G)

    # compute the elements of the stabilizer group H
    H = stabilizer_elements(F)
    H0 = copy(H) # save a copy of the stabilizer group for later
    # stabilizer elements are now in `H` with generator `L`
    # print(H0)

    gens = Matrix{ZZRingElem}[]     # generators
    ords = Int[]                    # orders
    setdiff!( G, H)

    while !isempty(G)
        p = popfirst!(G)
        push!( gens, p) # add this generator to the list
        k = 0
        old_length = 0
        while length(H) > old_length
            old_length = length(H)
            union!( H, map( h -> mod.(h*p, dd), H) )
            k += 1
        end
        push!( ords, k) # order of p
        setdiff!( G, H)  # trim elements from G
    end
    # we need the Smith normal form of this matrix
    A = diagonal_matrix(ZZ.(ords))

    # Fill in the relations to get the lower triangular part of A
    for k = 1:length(ords)-1
        p = mod.(gens[k+1]^ords[k+1],dd)
        j = -1
        r = zeros(Int,k)
        g = [0 0; 0 0]
        while !(g ∈ H0)
            j += 1
            r = radix( j , ords[1:k])
            g = mod.( foldl( (x,y) -> mod.(x*y, dd) , gens[1:k].^r[1:k]) * p, dd)
        end
        for j=1:k
            A[k+1,j] = r[j]
        end
    end
    # Now A offers complete information about the generators and relations

    # Smith normal form: U*A*V = D, where D is diagonal and D[j,j] | D[j+1,j+1].
    D, _, V = snf_with_transform(A)
    V = mod.(Int.(inv(V)), sizeofG) # make sure all entries are nonnegative
    q = [ foldl( (x,y) -> mod.(x*y, dd) , gens[:].^V[j,:] ) for j=1:length(gens)]

    # We only need diagonal entries > 1, so drop the rest
    D = diagonal(D)
    q = q[ D.>1 ]
    D = D[ D.>1 ]

    # Finally, split things further into elements with prime power order
    g = Matrix{ZZRingElem}[]
    n = Int[]
    for k=1:length(D)
        f = factor(D[k]).fac
        pp = Int.(keys(f).^values(f))
        P = prod(pp)
        for p in pp
            push!( g, mod.(q[k]^Int(P/p),dd) )
            push!( n, p)
        end
    end

    # output is a list of independent generators and their orders
    return g, n
end



@doc raw"""
    galois_elements(F::AdmissibleTuple)
    galois_elements(d::Integer, gens::Vector{Matrix{ZZRingElem}}, ords::Vector{Integer})

Compute one element from each coset of ``\mathcal{M}/\mathcal{S}`` where ``\mathcal{M}`` is a maximal abelian subgroup of ``\mathrm{GL}(2,\mathbb{Z}/d')`` and ``\mathcal{S}`` is the stabilizer of `F`.
The second form takes the dimension `d` and the generators and orders of generators for the galois group, as is computed by `galois_normal_form`.
The output is an array of `Matrix{ZZRingElem}` elements whose size is the the same as `ords`.
"""
function galois_elements(d::Integer, gens::Vector{Matrix{ZZRingElem}}, ords::Vector{Int64})
    dd = 2^iseven(d)*d
    n = prod(ords)
    T = reshape([ Matrix{ZZRingElem}(undef,2,2) for _ in 1:n ], ords...)
    for m in [ radix(k,ords) for k=0:n-1]
        T[(m.+1)...] .= mod.( prod(gens.^m), dd)
    end
    return T
end

galois_elements(F::AdmissibleTuple) = galois_elements(F.d,galois_normal_form(F)...)


@doc raw"""
    galois_order_orbit(F::AdmissibleTuple)

Compute the orders of the elements of the Galois group and a maximal orbit, output as a tuple.
"""
function galois_order_orbit(F::AdmissibleTuple)
    gens, ords = galois_normal_form(F)
    d = Int(F.d)
    dd = 2^iseven(d)*d
    T = galois_elements(d, gens, ords)
    p0 = galois_orbit(F)[1]
    porb = map( X -> Integer.(mod.(X*p0, dd)), T)
    return ords, porb
end
