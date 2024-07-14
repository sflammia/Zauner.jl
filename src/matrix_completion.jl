export matrix_completion

@doc raw"""
    matrix_completion( nu::AbstractArray, F::AdmissibleTuple
    [; verbose::Bool = false, shift::Integer = 0, test::Bool = false,
    max_iters = 1e5, eps_abs = 1e-7, eps_infeas = 1e-7, eps_rel = 1e-7])

Given the SIC phase overlaps `nu` on a maximal galois orbit, compute the associated SIC fiducial vector `ψ`.
If the input does not correspond to a valid SIC with admissible data `F`, then the output is unpredictable.
"""
function matrix_completion(nu::AbstractArray, F::AdmissibleTuple;
    verbose::Bool=false, shift::Integer=0, test::Bool=false,
    max_iters=1e5, eps_abs=1e-7, eps_infeas=1e-7, eps_rel=1e-7)
    dd = 2^iseven(F.d) * F.d
    p0 = galois_orbit(F)[1]
    stab = stabilizer_elements(F)
    gens, ords = galois_normal_form(F)
    # ords == size(nu)
    gnu = circshift(nu, radix(shift, ords))

    X = HermitianSemidefinite(F.d)
    obj = real(tr(X))
    cons = typeof(tr(X) == 0)[]
    # Possible idea to speed this up: only need about 4d - O(1) constraints
    # for injectivity with POVMs (Heinosaari, Mazzarella & Wolf).
    # Instead of a maximal orbit, we could sample a random sparse set
    # of constraints, say 6d of them to be safe.
    # Since the maximal orbit is generally quite large, sampling might be faster.
    # Could also potentially speed this up and remove dependencies by using fast non-convex methods (SVT, Burer-Montiero, FPCS, etc.) from the theory of low-rank matrix sensing.
    for m in [radix(k + shift, ords) for k = 0:prod(ords)-1]
        for s in stab
            nextp = Integer.(mod.(s * prod(gens .^ m) * p0, dd))
            push!(cons, tr(wh(-nextp, F.d) * X) == gnu[(1 .+ m)...])
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
    prob = minimize(obj, cons)
    Convex.solve!(prob, opt)
    XX = Convex.evaluate(X) / sqrt(F.d + 1)
    verbose && println("Norm difference ||ψ^2 - ψ||_2 = ", norm(XX * XX - XX))
    ψ = XX[:, 1] / norm(XX[:, 1])

    # if `test`, run a brute-force test of SIC overlap property
    if test
        sictest = maximum(abs.([abs(tr(ψ' * wh(p, ψ))) for p = 1:F.d^2-1] .- 1 / sqrt(F.d + 1)))
        verbose && println("Deviation of overlaps = $sictest")
    end
    return ψ
end
