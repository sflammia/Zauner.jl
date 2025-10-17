export necromancy

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
        notj = NTuple{r-1,Int}(setdiff(1:r, j))
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
    # Ensure that we have initialized the class field for F
    verbose && println(F)
    d = Int(F.d)
    ghostclassfield(F)
    signswitch(F)
    hb = lll(maximal_order(F.H)).basis_nf # find an LLL-reduced basis for H
    gb = F.g.(hb) # the Galois-conjugate basis
    eH = real_embeddings(F.H)[1] # fix a real embedding

    # the normal form orders of the Galois group and a maximal p-orbit
    ords, porb = galois_order_orbit(F)
    verbose && println("Galois group with orders ", ords)
    r, n = length(ords), prod(ords)

    # get a low-precision ghost
    prec = 64
    setprecision(BigFloat, prec; base=2)
    verbose && println("Computing the ghost.")
    ψ = (verbose ? (@time ghost(F)) : ghost(F))
    # initialize for the (potentially shifted) sic overlaps
    x = Array{Complex{BigFloat}}(undef, Tuple(ords)...)

    while prec ≤ max_prec
        prec *= 2
        # define some test flags
        finite_invariants = true
        complex_phases = true
        unique_intersections = true

        verbose && println("\n$prec bit target precision.")

        # bump up the precision
        verbose && println("Current working precision = ", precision(real(ψ[1])), " bits.")
        ψ = precision_bump(ψ, prec; base=2, verbose=verbose)
        ϕ = circshift(reverse(ψ), 1)
        ϕ .*= (BigFloat(d + 1)) / ϕ'ψ # include normalization factors
        verbose && println("new precision = ", precision(real(ψ[1])), " bits")

        # compute the ghost overlaps
        verbose && println("Computing the high-precision ghost overlaps.")
        (verbose ? (@time K = [real(ϕ'wh(p, ψ)) for p in porb])
         : K = [real(ϕ'wh(p, ψ)) for p in porb])

        # compute the ghost invariants
        verbose && println("Computing the ghost invariants.")
        a, b, s = _ghost_invariants(K)

        # sign-switch to the SIC invariants.
        verbose && println("Sign-switching to the SIC invariants.")
        # create a high precision embedding map and sign-switching automorphism
        _fH(x) = BigFloat.(real.(evaluation_function(eH, prec).(x)))
        primalbasis = _fH.(hb)
        dualbasis = _fH.(gb)
        _dual(x) = _dualize(primalbasis, dualbasis, x)
        # sign-switch and test for finiteness
        for j = 1:r
            a[j] = _dual.(a[j])
            finite_invariants &= all(isfinite.(a[j]))
            !finite_invariants && break

            b[j] = _dual.(b[j])
            finite_invariants &= all(isfinite.(b[j]))
            !finite_invariants && break

            s[j] = _dual.(s[j])
            s[j] = reverse(pow_to_elem_sym_poly(s[j]))
            finite_invariants &= all(isfinite.(s[j]))
            !finite_invariants && break
        end

        setprecision(BigFloat, 320; base=2)

        if !finite_invariants
            verbose && println("Some SIC invariants were not finite.\n    ...Doubling precision.")
            continue
        end
        verbose && finite_invariants && println("All SIC invariants are finite.")
        # From this point on we are in SIC world.

        # Lower the precision back to standard BigFloat plus a 64-bit buffer.
        # NOTE: This line forces a lot of recomputation with ghosts.
        # Basically, we need to precision bump from this baseline each time instead of the previous level of precision.
        # However, if we keep the high precision then it is very slow to compute roots.
        # setprecision(BigFloat, 320; base=2)

        # Compute c', map to elementary symmetric polynomials
        # then find roots to get K'
        θ = map(x -> -roots(BigFloat.(x)), s)
        L = Matrix{Complex{BigFloat}}[]
        for j = 1:r
            θprime = [θ[j][1]]
            for k = 2:ords[j]
                push!(θprime, dot(b[j], [θprime[k-1]^n for n = 0:(ords[j]-1)]))
            end
            push!(L, Vandermonde(θprime) * a[j])
            for t = 1:ords[j]
                L[j][t, :] = -roots(reverse(pow_to_elem_sym_poly(L[j][t, :])))
            end
            L[j] = L[j] / sqrt(BigFloat(d + 1))
            # These should all be approximate phases,
            # otherwise break and try again
            for ol in L[j]
                # complex_phases &= all(abs.(abs.(L[j]) - 1) < 10^(-10))
                complex_phases &= abs(abs(ol) - 1) < 10^(-10)
            end
            !complex_phases && break
        end
        if !complex_phases
            verbose && println("Computed overlaps were ", L)
            verbose && println("Some SIC overlaps were not complex phases.\n    ...Doubling precision.")
            continue
        end
        verbose && println("All SIC overlaps are complex phases.")

        # now intersect to get x, which is nu up to an unknown Galois action.
        for k = 0:n-1
            t = radix(k, ords) .+ 1
            Kt = Tuple([L[j][t[j], :] for j = 1:r])
            # intersect at 128-bit precision by default
            a = reduce((x, y) -> _approx_complex_intersection(x, y; prec=128), Kt)
            unique_intersections &= (length(a) == 1)
            if unique_intersections
                x[t...] = a[1]
            else
                verbose && println("Intersection error...\n    Doubling precision.")
                break
            end
        end
        !unique_intersections && continue
        verbose && println("All intersections are unique.\n")
        # the only way we can get to here is if
        # finite_invariants && complex_phases && unique_intersections
        break
    end # while

    # Return to standard BigFloat precision
    setprecision(BigFloat, 256; base=2)
    (prec > max_prec) && error("max_prec exceeded without convergence.")
    x = complex.(BigFloat.(real.(x)), BigFloat.(imag.(x)))

    # if we ran out of precision, return to standard BigFloat precision and err.

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
    return -dot(dual, t[1:end-1]) / t[end]
end


@doc """
    _approx_complex_intersection(A::AbstractVector, B::AbstractVector; prec::Integer = 256, base::Integer = 2)

Internal function to compute the intersection of two complex lists using a default tolerance of 256 bit precision.
"""
function _approx_complex_intersection(A::AbstractVector, B::AbstractVector; prec::Integer=256, base::Integer=2)

    scale = BigInt(base)^prec

    # helper function to scale and round a Complex{BigFloat} to Complex{BigInt}
    function scale_round(z)
        re = round(BigInt, real(z) * scale)
        im = round(BigInt, imag(z) * scale)
        return complex(re, im)
    end

    # scale and round both lists
    rounded_A = Set{Complex{BigInt}}(scale_round.(A))
    rounded_B = Set{Complex{BigInt}}(scale_round.(B))

    # intersection of rounded lists
    intersection_rounded = intersect(rounded_A, rounded_B)

    # convert back to original scale
    intersection = Set{Complex{BigFloat}}()
    for z in intersection_rounded
        # scale back
        original = complex(real(z) / scale, imag(z) / scale)
        push!(intersection, original)
    end

    # convert back to an array before returning
    return collect(intersection)
end
