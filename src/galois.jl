export galois_orbit, galois_normal_form

# Compute the maximal Galois orbit 
function galois_orbit(F::AdmissibleTuple)
    d = F.d
    if (d % 9) == 3 && rem(ZZ(F.f//F.q),3) == 0
        error("F_a orbits not supported yet")
        # this condition for F_a orbits is on p. 83 of main.tex, Thm 7.28.
    end
    # not supporting even dimensions yet.
    if iseven(d)
        error("Even dimensions not supported yet")
    end

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
    # if F is anti-unitary, this gives the generator.
    test, L = is_antiunitary_with_generator(F)
    # if F is unitary, then F.L is the stabilizer generator.
    if !test
        L = Matrix{Integer}(F.L)
    end
    stab = [ mod.(L,dd) ]
    k=1
    while stab[k] != [1 0; 0 1]
        push!(stab, mod.(L*stab[k], dd))
        k += 1
    end
    stab = stab[1:end-1] # don't keep the identity element

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


@doc raw"""
    centralizer_elements(F::AdmissibleTuple)

\\
Compute the elements of ``GL(2,ℤ/d')`` that are in the centralizer of the stability group for `F`. 


# Examples
First compute a ghost fiducial at 128-bit precision.
```jldoctest
julia> F = AdmissibleTuple(5)
AdmissibleTuple( d = 5, r = 1, K = ℚ(√12), Q = ⟨1,-4,1⟩ )

julia> centralizer_elements(F)
24-element Vector{Matrix{Integer}}:
 [1 0; 0 1]
 [2 0; 0 2]
 [3 0; 0 3]
 [4 0; 0 4]
 [4 4; 1 0]
 [0 4; 1 1]
 [1 4; 1 2]
 [2 4; 1 3]
 [3 4; 1 4]
 [3 3; 2 0]
 [4 3; 2 1]
 [0 3; 2 2]
 [1 3; 2 3]
 [2 3; 2 4]
 [2 2; 3 0]
 [3 2; 3 1]
 [4 2; 3 2]
 [0 2; 3 3]
 [1 2; 3 4]
 [1 1; 4 0]
 [2 1; 4 1]
 [3 1; 4 2]
 [4 1; 4 3]
 [0 1; 4 4]
```
"""
function centralizer_elements(F::AdmissibleTuple)
    # Need a caveat here about F_a orbits. Will add support later.
    d = F.d
    # this condition for F_a orbits is on p. 83 of main.tex, Thm 7.28.
    rem(d,9) == 3 && rem(ZZ(F.f//F.q),3) == 0 && error("F_a orbits not supported yet")
    # not supporting even dimensions yet.
    iseven(d) && error("Even dimensions not supported yet")

    dd = 2^iseven(d) * d #  d'
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
    galois_normal_form(F::AdmissibleTuple)

\\
Compute the normal form of the group `M/S` where `M` is a maximal abelian subgroup of `GL(2,ℤ/d')` and `S` is the stabilizer of `F`. 
This is group is (conjecturally) isomorphic to the Galois group. 
The output is a tuple `(g,n)` where `g` is a vector of generator matrices and `n` is a vector of the orders of those matrices in `M/S`. 
This corresponds to the canonical decomposition `ℤ/n_1 ⊕ ℤ/n_2 ⊕ ... ⊕ ℤ/n_k` where the `n_j` are each prime powers. 
The factors are canonical up to permutation (though the generators are not). 

# Examples
```jldoctest
julia> F = AdmissibleTuple( 7, QuadBin(2,-4,1))
AdmissibleTuple( d = 7, r = 1, K = ℚ(√8), q = 1, Q = ⟨2,-4,1⟩, h = 1 )

julia> g, n = galois_normal_form(F)
(Matrix{ZZRingElem}[[6 0; 0 6], [2 0; 0 2]], [2, 3])

julia> F = AdmissibleTuple(9)
AdmissibleTuple( d = 9, r = 1, K = ℚ(√60), q = 1, Q = ⟨1,-8,1⟩, h = 2 )

julia> g, n = galois_normal_form(F)
(Matrix{ZZRingElem}[[8 5; 4 3], [8 0; 0 8], [7 0; 0 7]], [3, 2, 3])
```
"""
function galois_normal_form(F::AdmissibleTuple)
    d = F.d
    if (d % 9) == 3 && rem(ZZ(F.f//F.q),3) == 0
        error("F_a orbits not supported yet")
        # this condition for F_a orbits is on p. 83 of main.tex, Thm 7.28.
    end
    # not supporting even dimensions yet.
    if iseven(d)
        error("Even dimensions not supported yet")
    end

    dd = 2^iseven(d)*d

    # G is the centralizer of H (the stabilizer)
    # compute the nontrivial centralizer orbits
    G = Matrix{ZZRingElem}.(centralizer_elements(F))
    sizeofG = length(G)
    
    # compute the elements of the stabilizer group H
    # if F is anti-unitary, this gives the generator.
    test, L = is_antiunitary_with_generator(F)
    
    # if F is unitary, then F.L is the stabilizer generator.
    if !test
        L = F.L
    end
    H = [ mod.(L,dd) ]
    k=1
    while H[k] != [1 0; 0 1]
        push!(H, mod.(L*H[k], dd))
        k += 1
    end
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
            push!( g, mod.(q[k]^Int(P/p),F.d) )
            push!( n, p)
        end
    end

    # output is a list of independent generators and their orders
    return g, n
end

