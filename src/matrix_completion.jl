export matrix_completion

@doc raw"""
    matrix_completion( nu::AbstractArray, F::AdmissibleTuple)

\\
Given the SIC phase overlaps `nu`, compute the associated SIC fiducial vector ψ. 
If the input does not correspond to a valid SIC with admissible data `F`, then the output is unpredictable.

(Need to try `circshift(nu,(a,b,c))` in general, where `(a,b,c)` is a tuple of shifts of the appropriate size for the Galois group of the SIC with data `F`.)

# Examples
First compute the ghost for `d = 5`.
```jldoctest
F = AdmissibleTuple(5)
1+1

# output

2

```
"""
function matrix_completion( nu::AbstractArray, F::AdmissibleTuple; 
        verbose::Bool = false, shift::Integer = 0, test::Bool = false, 
        max_iters = 1e5, eps_abs = 1e-7, eps_infeas = 1e-7, eps_rel = 1e-7)
    dd = 2^iseven(F.d)*F.d
    p0 = galois_orbit(F)[1]
    stab = stabilizer_elements(F)
    gens, ords = galois_normal_form(F)
    # ords == size(nu)
    gnu = circshift( nu, radix( shift, ords)) 

    X = HermitianSemidefinite(F.d)
    obj = real(tr(X))
    cons = Convex.EqConstraint[ ]
    # Possible idea to speed this up: only need about 4d - O(1) constraints 
    # for injectivity with POVMs (Heinosaari, Mazzarella & Wolf). 
    # Instead of a maximal orbit, we could sample a random sparse set 
    # of constraints, say 6d of them to be safe. 
    # Since the maximal orbit is generally quite large, sampling might be faster. 
    for m in [ radix(k+shift,ords) for k=0:prod(ords)-1 ]
        for s in stab
            nextp = Integer.(mod.( s*prod(gens.^m)*p0, dd))
            push!( cons, tr(wh(-nextp,F.d)*X) == gnu[(1 .+ m)...])
        end
    end
    opt = Convex.MOI.OptimizerWithAttributes( 
                SCS.Optimizer, 
                Convex.MOI.Silent() => !verbose,
                "max_iters" => max_iters,
                "eps_abs" => eps_abs,
                "eps_infeas" => eps_infeas,
                "eps_rel" => eps_rel,
          )
    prob = minimize(obj,cons)
    Convex.solve!(prob, opt) 
    XX = Convex.evaluate(X)/sqrt(F.d+1)
    verbose && println("Eigenvalues λ_1 and λ_2 are ",sort(real.(eigvals(XX)))[[end;end-1]])
    ψ = XX[:,1]/norm(XX[:,1])

    # if `test`, run a brute-force test of SIC overlap property
    if test
        sictest = maximum(abs.([ abs(tr(ψ'*wh(radix(p,[F.d;F.d]),ψ))) for p = 1:F.d^2-1 ] .- 1/sqrt(F.d+1)))
        verbose && println("Deviation of overlaps = $sictest")
    end
    return ψ
end