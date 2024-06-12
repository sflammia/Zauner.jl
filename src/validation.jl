export sic_overlap_test, ghost_overlap_test, sic_frame_test, ghost_frame_test

@doc raw"""
    sic_overlap_test(ψ::AbstractVector)

Check the SIC equiangularity conditions by returning
```math
\max_{\boldsymbol{p} \not=\boldsymbol{0}} \bigl|\nu_{\boldsymbol{p}} \nu_{-\boldsymbol{p}} - \tfrac{1}{d+1}\bigr|\,,
```
where ``\nu_{\boldsymbol{p}} = \langle\psi|D_{\boldsymbol{p}}|\psi\rangle`` are the SIC overlaps.
"""
function sic_overlap_test(ψ::AbstractVector)
    d = length(ψ)
    t = eltype(ψ)
    maximum( [ abs(abs2(ψ'*wh(p,ψ)) - one(t)/(d+1)) for p=1:d^2-1])
end


@doc raw"""
    ghost_overlap_test(ψ::AbstractVector)

Check the ghost overlap conditions.
If all ghost overlaps are approximately real, it returns
```math
\max_{\boldsymbol{p} \not=\boldsymbol{0}} \bigl|\nu_{\boldsymbol{p}} \nu_{-\boldsymbol{p}} - \tfrac{1}{d+1}\bigr|\,.
```
If they aren't approximately real it throws an error.
"""
function ghost_overlap_test(ψ::AbstractVector)
    d = length(ψ)
    t = eltype(ψ)
    ϕ = circshift(reverse(ψ),1)
    ov = [ϕ'*wh(p,q,ψ)*ϕ'*wh(-p,-q,ψ) for p=0:d-1, q=0:d-1]./(ϕ'ψ)^2

    # check approximate reality
    @assert all(ov .≈ real.(ov)) "Ghost overlaps must be real."
    ov = real.(ov)

    # quantify the deviation
    maximum( abs.(ov[2:end] .- one(t)/(d+1) ) )
end


@doc raw"""
    sic_frame_test(ψ::AbstractVector)

Return the absolute deviation of the pointwise conditions on the frame potential from eq. 8 of arXiv:0707.2071.
Let ``T(k,l) = \sum_j \psi_{j}\psi_{j+k}^* \psi_{j+l}^* \psi_{j+k+l}``, and recall that ``\psi`` is a SIC if and only if ``T(k,l) - \frac{\delta_{k,0}+\delta_{l,0}}{d+1} = 0``.
The function returns the maximum absolute deviation from these conditions.
"""
function sic_frame_test(ψ::AbstractVector)
    d = length(ψ)
    t = eltype(ψ)
    n = (ψ'ψ)^2/(d+one(t)) # normalizing factor

    # The function T below has an eight-fold symmetry generated by the transformations:
    #    (k,l) → ( l, k)
    #    (k,l) → (-k,-l)
    #    (k,l) → ( k,-l)*
    # where the last one means apply the symmetry then complex conjugate the result.
    # It therefore suffices to check a subset of the equations.
    function T(k,l)
        sum(ψ .* conj(circshift(ψ, k)) .* conj(circshift(ψ, l) ) .* circshift(ψ, k+l)) - ((k==0) + (l==0)) * n
    end

    dev = abs(zero(t))
    for j = 0:(d+1)÷2, k = 0:j
        dev = max( abs( T(j,k) ), dev)
    end
    dev
end


@doc """
    ghost_frame_test(ψ::AbstractVector)

Return the maximum absolute deviation of the pointwise frame conditions.
This is the ghost analog of eq. 8 of arXiv:0707.2071.
"""
function ghost_frame_test(ψ::AbstractVector)
    d = length(ψ)
    t = eltype(ψ)
    ϕ = circshift(reverse(ψ),1)
    n = (ϕ'ψ)^2/(d+one(t)) # normalizing factor

    function T(k,l)
        sum(ψ .* conj(circshift(ϕ, k)) .* conj(circshift(ϕ, l) ) .* circshift(ψ, k+l)) - ((k==0) + (l==0)) * n
    end

    dev = abs(zero(t))
    for j = 0:(d+1)÷2, k = 0:j
        dev = max( abs(T(j,k)), dev)
    end

    dev
end
