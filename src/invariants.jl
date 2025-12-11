export necromancy

function ghost_orbit(F::AdmissibleTuple, prec::Integer;
    base::Integer=10,
    verbose::Bool=false)

    verbose && println(F)

    # the normal form orders of the Galois group and a maximal p-orbit
    ords, porb = galois_order_orbit(F)
    verbose && println("Galois group with orders ", ords)

    # get a low-precision ghost
    INIT_PREC = 64
    verbose && println("Computing the initial ghost to low $(INIT_PREC) bits.")
    setprecision(BigFloat, INIT_PREC; base=2)
    verbose ? (@time ψ = ghost(F)) : ψ = ghost(F)

    verbose && println("Precision of initial ghost is ", precision(real(ψ[1])), " bits")

    ψ = precision_bump(ψ, prec; base=base, verbose=verbose)
    ϕ = circshift(reverse(ψ), 1)
    ϕ .*= (BigFloat(F.d + 1)) / ϕ'ψ # include normalization factors
    verbose && println("Ghost precision is now ", precision(real(ψ[1])), " bits")

    # compute the ghost overlaps
    verbose && println("Computing the high-precision ghost overlaps.")
    (verbose ? (@time K = [real(ϕ'wh(p, ψ)) for p in porb]) : K = [real(ϕ'wh(p, ψ)) for p in porb])

    return K
end


@doc raw"""
    _ghost_invariants(K::AbstractArray{BigFloat})

Internal function to compute numerical approximations to the ghost invariants starting from the array `K` of ghost overlaps.
"""
function _ghost_invariants(K::AbstractArray{BigFloat})
    ords = size(K)
    r, n = length(ords), prod(ords)
    prec = precision(K[1])
    THRESHOLD = BigFloat(2)^(10 - prec / 2)

    V = Vandermonde{BigFloat}
    a = Matrix{BigFloat}[]
    b = Vector{BigFloat}[]
    s = Vector{BigFloat}[]
    for j = 1:r
        notj = NTuple{r - 1,Int}(setdiff(1:r, j))
        # first compute "l_j".
        for l = 1:n÷ords[j]
            c = dropdims(sum(K .^ l; dims=notj); dims=notj)
            # If these are distinct and nonzero then we've found l_j.
            # We test this with at least 10 bits above 1/2 of our base precision.
            if _distinct_nonzero(c, THRESHOLD)
                V = Vandermonde(c)
                push!(b, V \ circshift(c, -1))
                # note: s is called "e" in the draft, but that conflicts with e(z).
                # These are the power sums in the c_{j,t,l_j} over each t.
                push!(s, [sum(c .^ k) for k = 1:ords[j]])
                break
            end
        end
        # now compute a
        push!(a, zeros(BigFloat, (ords[j], n ÷ ords[j])))
        for l = 1:n÷ords[j]
            x = dropdims(sum(K .^ l; dims=notj); dims=notj)
            a[j][:, l] = V \ x
        end
    end

    return a, b, s
end


@doc """
    _distinct_nonzero(v::Vector{<:T},prec::T) where T<:Real

Internal function to test if the elements in a real vector `v` are distinct and nonzero to precision `prec`.
"""
function _distinct_nonzero(v::Vector{<:T}, prec::T) where {T<:Real}
    all(abs.([diff(sort(v)); v]) .> prec)

function class_field_bases(F::AdmissibleTuple, prec::Int)
    # Make sure that the class field and sign-switching automorphism are initialized
    ghostclassfield(F)
    signswitch(F)
    hb = lll_basis(maximal_order(F.H)) # LLL-reduced basis for H
    gb = F.g.(hb) # the Galois-conjugate basis

    eH = real_embeddings(F.H)[1] # fix a real embedding

    # Evaluate at the embedding eH.
    _fH(x) = BigFloat.(real.(evaluation_function(eH, prec).(x)))
    primalbasis = _fH.(hb)
    dualbasis = _fH.(gb)

    # return numerical bases to precision prec.
    return primalbasis, dualbasis
end



function pseudonecromancy(F::AdmissibleTuple, prec::Int; verbose=false, base=10)
    u = ghost_orbit(F, prec; verbose=verbose, base=base)
    x, p, y = ghost_basis(u; verbose=verbose)
    A = Zauner._multiplication_matrix(x, p, y)
    szu = size(u)

    # sign-switch A and construct z
    Hprimal, Hdual = class_field_bases(F, prec)

    z = similar.(p)
    for i = 1:length(z)
        B = Matrix{Float64}(undef, size(A[i])...)
        for t in CartesianIndices(A[i])
            # println("t = $t")
            R = guess_int_null_vec([Hprimal; A[i][t]], 3 * prec ÷ 4)
            # println(R)
            R2 = guess_int_null_vec([Hprimal; A[i][t]], prec)
            # println(R2)
            if R == R2 || R == -R2
                B[t] = Float64(-dot(Hdual, R[1:end-1]) / R[end])
            else
                error("Unreliable sign-switch of A matrix.")
            end
        end
        vals, vecs = eigen(B)
        z[i] = (vals[1] / vecs[1, 1]) .* vecs[:, 1]
    end

    Y = (length(y) == 1 ? y[1] : kron(reverse(y)...))
    s = zeros(Float64, szu)

    for t = 0:prod(szu)-1
        shift = radix(t, szu)
        c = dot(circshift(u, shift)[:], Y)

        if log10.(abs.(c)) .< -prec
            continue
        end
        R = guess_int_null_vec([Hprimal; c], 3 * prec ÷ 4)
        R2 = guess_int_null_vec([Hprimal; c], prec)

        if R == R2 || R == -R2
            s[(1 .+ shift)...] = Float64((Hdual' * R[1:end-1]) / (-R[end]))
        else
            error("Unreliable sign-switch of c vector.")
        end
    end

    Z = ComplexF64.(length(z) == 1 ? z[1] : kron(reverse(z)...))
    uu = zeros(ComplexF64, szu)
    for t = 0:prod(szu)-1
        shift = radix(t, szu)
        uu[(1 .+ shift)...] = dot(circshift(s, shift)[:], Z)
    end

    @assert all(abs2.(uu) .≈ F.d + 1) "Reconstructed SIC units are not proportional to complex phase units, √(d+1)×exp(iθ)."

    return uu
end

function pseudonecromancy(
    F::AdmissibleTuple;
    init_prec::Integer=2^8,
    max_prec::Integer=2^14,
    base=10,
    verbose=false,
)

    basename = _base_name(base)
    if init_prec < log(base, 10) * 2^6
        verbose && @info "Increasing initial precision from $init_prec $basename to $(2^6) digits."
    end
    prec = max(init_prec, round(Integer, log(base, 10) * 2^6, RoundUp))

    verbose && @info "Starting pseudonecromancy computation with initial precision $prec $basename."
    while prec ≤ max_prec
        try
            u = pseudonecromancy(F, prec; verbose=verbose, base=base)
            return u
        catch err
            verbose && @warn "Failed to compute pseudonecromancy for precision $prec, with error:\n$(err)\n Doubling precision and trying again."
        end
        prec *= 2
    end

    error("Failed to achieve pseudonecromancy with $max_prec $basename. Try increasing `max_prec`.")
end


@doc raw"""
    necromancy( F::AdmissibleTuple [; max_prec = 2^20, verbose = false])

Compute a numerical approximations to a SIC associated to the `AdmissibleTuple` `F`.
The maximum number of bits used in integer relation finding is set to `max_prec` (default of 1 Mb) and `verbose` can be toggled `true` or `false`.

# Examples

Check that the principal SIC in ``d=7`` satisfies the equiangularity conditions.

```
julia> d = 7; F = AdmissibleTuple(d)
AdmissibleTuple( d = 7, K = ℚ(√8), q = 2, Q = ⟨1,-6,1⟩, h = 1 )

julia> ψ = necromancy(F);

julia> all([ abs2(ψ'wh(p,ψ)) for p=1:d^2-1] .≈ 1/(d+1))
true
```
"""
function necromancy(F::AdmissibleTuple;
    max_prec::Integer=2^20,
    overlap_tol::Float64=1e-6,
    overlap_prec::Integer=256,
    base::Integer=2,
    verbose::Bool=false)

    @warn "This function is being updated and is presently nonfunctional."
    return nothing

    # Now try every shift in the Galois group until one of them gives a SIC
    verbose && println("Now searching through Galois shifts and using matrix completion.")
    for k = 0:n-1
        ψ = matrix_completion(circshift(x, radix(k, ords)), F)
        sft = sic_frame_test(ψ)
        if sft < overlap_tol
            verbose && println("Fiducial vector found with frame equations correct to ≤ $sft.")
            break
        end
    end

    # Finally, improve the precision again to standard BigFloat.
    z = re_im_proj(Complex{BigFloat}.(ψ))
    buffer_bits = 10
    precision_bump!(z, _sic_olp_func, overlap_prec + buffer_bits; base=base, verbose=verbose)
    ψ = re_im_proj(z)
    ψ = ψ / sqrt(ψ'ψ)
    setprecision(BigFloat, overlap_prec)
    z = reim(ψ)
    # round down the the desired precision so that the new precision of BigFloat is predictable
    return BigFloat.(z[1]) + im * BigFloat.(z[2])
end



@doc """
    _dualize( primal::Vector{T}, dual::Vector{T}, x::T) where T::AbstractFloat

Internal function that takes an `AbstractFloat` number `x`, rounds it into the `primal` basis, then expands it again in the `dual` basis.
If `x` is not faithfully represented in the `primal` basis then the result is unpredictable.
If `primal` and `dual` are related by a galois automorphism `g`, then ideally this outputs an approximation of `g(x)`.
"""
function _dualize(primal::Vector{T}, dual::Vector{T}, x::T) where {T<:AbstractFloat}
    t = guess_int_null_vec([primal; x])
    # println("‖t‖₁ = ", sum(abs.(t)))
    return -dot(dual, t[1:end-1]) / t[end]
end
