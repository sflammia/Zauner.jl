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
    d = Int(F.d)
    dd = 2^iseven(d) * d
    p0 = galois_orbit(F)[1]
    stab = stabilizer_elements(F)
    gens, ords = galois_normal_form(F)
    # ords == size(nu)
    gnu = circshift(nu, radix(shift, ords))

    X = HermitianSemidefinite(d)
    obj = real(tr(X))
    cons = typeof(tr(X) == 0)[]

    for m in [radix(k + shift, ords) for k = 0:prod(ords)-1]
        for s in stab
            nextp = Integer.(mod.(s * prod(gens .^ m) * p0, dd))
            push!(cons, tr(wh(-nextp, d) * X) == gnu[(1 .+ m)...])
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
    XX = Convex.evaluate(X) / sqrt(d + 1)
    verbose && println("Norm difference ||ψ^2 - ψ||_2 = ", norm(XX * XX - XX))
    ψ = XX[:, 1] / norm(XX[:, 1])

    # if `test`, run a brute-force test of SIC overlap property
    if test
        sictest = maximum(abs.([abs(tr(ψ' * wh(p, ψ))) for p = 1:d^2-1] .- 1 / sqrt(d + 1)))
        verbose && println("Deviation of overlaps = $sictest")
    end
    return ψ
end
