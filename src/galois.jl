export galois_orbit

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
    # NOTE: You only need the Zauner matrix to compute the centralizer. 
    # The stabilizer needs L of course. 
    # Need a caveat here about F_a orbits. Will add support later.
    d = F.d
    if rem(d,9) == 3 && rem(ZZ(F.f//F.q),3) == 0
        error("F_a orbits not supported yet")
        # this condition for F_a orbits is on p. 83 of main.tex, Thm 7.28.
    end
    # not supporting even dimensions yet.
    if rem(d,2) == 0
        error("Even dimensions not supported yet")
    end

    dd = 2^iseven(d) * d #  d'
    Z   = [0 -1;  1  0]*qmat(F.Q) # this is one choice of Zauner; could use [0 -1; 1 -1] instead
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




# This isn't used currently
function galois_reps(F::AdmissibleTuple)
    # Need a caveat here about F_a orbits. Will add support later.
    d = F.d
    if rem(d,9) == 3 && rem(ZZ(F.f//f.q),3) == 0
        error("F_a orbits not supported yet")
        # this condition for F_a orbits is on p. 83 of main.tex, Thm 7.28.
    end
    # not supporting even dimensions yet.
    if rem(d,2) == 0
        error("Even dimensions not supported yet")
    end
    
    # Need to test for anti-unitary symmetry
    # Def 7.20, Lemma 7.21, and Def. 4.31
    if is_antiunitary(F)
        error("Anti-unitary symmetry not supported yet")
    end

    # first compute the elements of the centralizer
    cents = centralizer_elements(F)

    # compute the elements of the stabilizer group
    L = Matrix{Integer}(F.L)
    dd = 2^iseven(d)*d
    stab = [ mod.(L,dd) ]
    k=1
    while stab[k] != [1 0; 0 1]
        push!(stab, mod.(L*stab[k], dd))
        k += 1
    end

    reps = Matrix{Integer}[]
    while !isempty(cents)
        M = popfirst!(cents)
        push!( reps, M)
        orb = map( x -> mod.(x*M,dd) , stab)
        setdiff!( cents, orb)
    end
    reps
end
