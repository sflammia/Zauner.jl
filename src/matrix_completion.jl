export initial_p_orbit, initial_B_estimate, rank_one_least_squares, shift_search


@doc raw"""
    initial_p_orbit(F::AdmissibleTuple; verbose::Bool = false)

Output a maximal Galois orbit of p-tuples in the shape of the Galois group, with the stabilizer group added in front.
"""
function initial_p_orbit(F::AdmissibleTuple;
    verbose::Bool=false)
    d = Int(F.d)
    dd = 2^iseven(d) * d
    p0 = galois_orbit(F)[1]
    verbose && println("Initial p-tuple:", p0)
    stab = stabilizer_elements(F)
    gens, ords = galois_normal_form(F)
    S = length(stab)
    p_size = (S, ords...)

    p = Array{Tuple{Integer,Integer}}(undef, p_size...)
    idx = [radix(k, ords) for k = 0:prod(ords)-1]
    for m in idx
        for s = 1:S
            p[s, (m .+ 1)...] = Tuple(Int.(mod.(stab[s] * prod(gens .^ m) * p0, dd)))
        end
    end

    return p
end

@doc raw"""
    initial_B_estimate(u::AbstractArray,
    p::AbstractArray{NTuple{2,I},N},
    d::Integer;
    verbose::Bool=false,
    shift::Integer=0) where {I<:Integer,N}

Output an initial estimate hermitian, unit-trace matrix.
"""
function initial_B_estimate(
    u::AbstractArray,
    p::AbstractArray{NTuple{2,I},N},
    d::Integer;
    verbose::Bool=false,
    shift::Integer=0
) where {I<:Integer,N}

    T = eltype(u)
    X = wh((0, 0), d, T) / d
    nu = circshift(u, radix(shift, size(u)))
    for a in axes(p, 1)
        @inbounds for idx in CartesianIndices(nu)
            pix = p[a, idx]
            nuval = nu[idx]
            X .+= nuval / (d * (d + 1)) .* wh(pix, d, T)
        end
    end
    return (X + X') / 2
end


"""
    rank_one_least_squares(p, B; maxit=1000, η=1e-1, tol=1e-12)

Solve the least squares problem via gradient descent on the complex unit sphere.
Here `p` is a list of observed Weyl-Heisenberg displacements which is assumed to be closed under reflection by -1 (to ensure hermitian output).
The input `B` is a hermitian matrix such that `b[p] = tr(B * wh(-p, d))`.
"""
function rank_one_least_squares(
    p::AbstractArray{NTuple{2,I},N},
    B::AbstractMatrix;
    maxit=1000, η=1e-1, tol=1e-12,
    verbose=false
) where {I<:Integer,N}

    T = eltype(B)
    d = size(B, 1)
    y = LinearAlgebra.eigvecs(B)[:, end]

    # reusable work space
    a = similar(y)
    b = similar(y)
    s = similar(y)

    function wh_grad!(y)
        s = zeros(T, d)
        @inbounds for q in p
            wh!(b, -q[1], -q[2], y)
            c = (y' * b) / d
            s .+= c * wh!(a, q[1], q[2], y) + conj(c) * b
        end
        return s
    end

    By = similar(y)

    for it in 1:maxit
        # gradient wrt y; simplified since B is hermitian
        LinearAlgebra.mul!(By, B, y)
        y_new = y .- η .* (wh_grad!(y) - 2 * By)
        y_new /= LinearAlgebra.norm(y_new)
        if it % 10 == 0
            del = LinearAlgebra.norm(y_new - y)
            verbose && println("At iteration $it, ‖y_new-y‖ = ", del)
            if del < tol
                y = y_new
                break
            end
        end
        y = y_new
    end

    y = cis(-angle(y[1])) .* y
    return y
end


function shift_search(
    u::AbstractArray,
    p::AbstractArray{NTuple{2,I},N},
    d::Integer;
    olp_goal=1e-10,
    tol=1e-12,
    maxit=1000,
    verbose::Bool=false
) where {I<:Integer,N}

    shift = 0
    test_val = 1.0
    v = zeros(eltype(u), d)
    while shift < prod(size(u))
        B = initial_B_estimate(u, p, d; shift=shift)
        v = rank_one_least_squares(p, B; tol=tol, maxit=maxit, verbose=verbose)
        test_val = sic_overlap_test(v)
        # println("shift overlap error = ", round(test_val, sigdigits=5))
        if test_val < olp_goal
            return v, shift
        end
        shift += 1
    end
    error("No solution found")
end
