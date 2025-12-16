export necromancy, pseudonecromancy, ghost_orbit, class_field_bases


function ghost_orbit(
    F::AdmissibleTuple,
    ψ::Vector{Complex{BigFloat}}
)
    ords, porb = galois_order_orbit(F)

    prec = precision(real(ψ[1]))
    setprecision(BigFloat, prec) do
        ϕ = circshift(reverse(ψ), 1)
        ϕ .*= BigFloat(F.d + 1) / (ϕ'ψ)
        K = [real(ϕ' * wh(p, ψ)) for p in porb]
        return reshape(K, ords...)
    end
end


function _ensure_ghost!(
    F::AdmissibleTuple,
    ψ::Union{Nothing,Vector{Complex{BigFloat}}},
    target_prec::Int;
    base=10,
    verbose=false
)
    if ψ === nothing
        verbose && println("Computing an initial ghost at low precision.")
        # initial precision is set to 64 bits
        setprecision(BigFloat, 64) do
            ψ = ghost(F)
        end
    end

    ψ = precision_bump(ψ, target_prec; base=base, verbose=verbose)
    return ψ
end


# Inner loop version
function pseudonecromancy(
    F::AdmissibleTuple,
    prec::Int;
    ghost_cache::Union{Nothing,Vector{Complex{BigFloat}}}=nothing,
    recipes::Union{Nothing,Vector{NormalBasisRecipe}}=nothing,
    verbose=false,
    base=10 # possible inconsistency in default base
)

    Ψ_cache = _ensure_ghost!(F, ghost_cache, prec; verbose=verbose, base=base)
    u = ghost_orbit(F, Ψ_cache)
    verbose && println("Precisions of ghost orbit: $(precision(u[1])) and precision of Ψ_cache = $(precision(real(Ψ_cache[1])))")

    setprecision(precision(u[1])) do
        if recipes === nothing
            x, p, y, recipes = ghost_basis(u; verbose=verbose)
        else
            x, p, y, _ = ghost_basis(u; recipes=recipes, verbose=verbose)
        end
        verbose && println("Precisions of ghost basis: $(precision(x[1][1])) and $(precision(y[1][1]))")
        A = Zauner._multiplication_matrix(x, p, y)
        szu = size(u)

        # sign-switch A and construct z
        # Note: `prec` should be in bits here.
        Hprimal, Hdual = class_field_bases(F, prec)

        z = map(x -> similar(x, ComplexF64), p)
        for i = 1:length(z)
            B = Matrix{Float64}(undef, size(A[i])...)
            for t in CartesianIndices(A[i])
                B[t] = Float64(_dualize(Hprimal, Hdual, A[i][t], prec))
            end
            vals, vecs = eigen(B)
            z[i] = (vals[1] / vecs[1, 1]) .* vecs[:, 1]
        end

        Y = (length(y) == 1 ? y[1] : kron(reverse(y)...))
        s = zeros(Float64, szu)

        # This is a convolution and could be optimized by using FFTs
        for t = 0:length(u)-1
            shift = radix(t, szu)
            c = dot(circshift(u, shift)[:], Y)

            if log10.(abs.(c)) .< -prec
                continue # s is zero, no need to dualize
            end
            s[(1 .+ shift)...] = Float64(_dualize(Hprimal, Hdual, c, prec))
        end

        # This could be optimized by using FFTs
        Z = (length(z) == 1 ? z[1] : kron(reverse(z)...))
        uu = zeros(ComplexF64, szu)
        for t = 0:length(u)-1
            shift = radix(t, szu)
            uu[(1 .+ shift)...] = dot(circshift(s, shift)[:], Z)
        end

        @assert all(abs2.(uu) .≈ F.d + 1) "Reconstructed SIC units are not proportional to complex phase units, √(d+1)×exp(iθ)."

        return uu, recipes, Ψ_cache
    end
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
    prec = max(init_prec, ceil(Int, log(base, 10) * 2^6))

    ψ_cache = nothing
    recipes = nothing

    verbose && "Starting pseudonecromancy computation with initial precision $prec $basename."
    while prec ≤ max_prec
        try
            u, recipes, ψ_cache = pseudonecromancy(
                F,
                prec;
                ghost_cache=ψ_cache,
                recipes=recipes,
                verbose=verbose,
                base=base
            )
            return u
        catch err
            verbose && @warn "Failed to compute pseudonecromancy for precision $prec, with error:\n$(err)\n Doubling precision and trying again."
        end
        prec *= 2
    end

    error("Failed to achieve pseudonecromancy with $max_prec $basename. Try increasing `max_prec`.")
end


@doc raw"""
    _multiplication_matrix(
        x::Vector{Vector{BigFloat}},
        p::Vector{Vector{Complex{BigFloat}}},
        y::Vector{Vector{BigFloat}}
    )

Internal function to compute the multiplication matrix for the normal basis `x`.
If the dual basis `y` is not provided, it is computed using the inverse Fourier transform of the reciprocal of the Fourier transform of `x`.
"""
function _multiplication_matrix(
    x::Vector{Vector{T}},
    p::Vector{Vector{Complex{T}}},
    y::Vector{Vector{T}}
) where {T<:Number}
    r = length(x)
    @assert length(p) == r && length(y) == r "The dimensions of x, p, and y must match."
    A = Vector{Matrix{T}}(undef, r)
    for j = 1:r
        Ty = hcat([circshift(y[j], -k) for k = 0:length(y[j])-1]...)
        A[j] = real.(fft(p[j] .* fft(x[j] .* Ty, 1), 1))
    end

    return A
end

function _multiplication_matrix(x::Vector{<:Number})
    n = length(x)
    p = ifft(x)
    y = real.(ifft(1 ./ (n .* p)))
    Ty = hcat([circshift(y, -k) for k = 0:n-1]...)
    A = real.(fft(p .* fft(x .* Ty, 1), 1))
    return A
end



function class_field_bases(F::AdmissibleTuple, bits::Integer)
    # Make sure that the class field and sign-switching automorphism are initialized
    ghostclassfield(F)
    signswitch(F)
    hb = lll_basis(maximal_order(F.H)) # LLL-reduced basis for H
    gb = F.g.(hb) # the Galois-conjugate basis

    eH = real_embeddings(F.H)[1] # fix a real embedding

    # Evaluate at the embedding eH.
    _fH(x) = BigFloat.(real.(evaluation_function(eH, bits).(x)))
    primalbasis = _fH.(hb)
    dualbasis = _fH.(gb)

    # return numerical bases to precision bits.
    return primalbasis, dualbasis
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
    _dualize( primal::Vector{T}, dual::Vector{T}, x::T, prec::Integer) where T::AbstractFloat

Internal function that takes an `AbstractFloat` number `x`, rounds it into the `primal` basis, then expands it again in the `dual` basis.
If `x` is not reliably reconstructed in the `primal` basis at `prec` and `3 * prec ÷ 4`, then it throws an error.
If `primal` and `dual` are related by a galois automorphism `g`, then ideally this outputs an approximation of `g(x)`.
"""
function _dualize(primal::Vector{T}, dual::Vector{T}, x::T, prec::Integer) where {T<:AbstractFloat}
    R = guess_int_null_vec([primal; x], 3 * prec ÷ 4)
    R2 = guess_int_null_vec([primal; x], prec)

    if R == R2 || R == -R2
        return (dual' * R[1:end-1]) / (-R[end])
    else
        error("Unreliable integer relations: dualization failed.")
    end
end
