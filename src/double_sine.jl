export double_sine


@doc raw"""
    double_sine( x, ω1, ω2 [; kwargs...])

The double sine function, defined in terms of the double gamma function following Shintani's convention:
```math
S_2(x,\boldsymbol{\omega}) = \frac{\Gamma_2(x,\boldsymbol{\omega})}{\Gamma_2(\omega_1+\omega_2-x,\boldsymbol{\omega})}
```

This is currently implemented for real values only.

When ``0 < x < ω_1 + ω_2``, ``\log S_2`` has the integral representation
```math
    \log S_2(x,\boldsymbol{ω}) = \frac12 \int_0^\infty \left(\frac{\sinh\bigl((ω_1 + ω_2 - 2x)t\bigr)}{\sinh(ω_1 t)\sinh(ω_2 t)} - \frac{ω_1 + ω_2 - 2 x}{ω_1 ω_2 t}\right) \frac{\mathrm{d}t}{t}
```
Evaluation proceeds via numerical integration using adaptive Gauss quadrature using the optional keyword arguments for `QuadGK`.

If ``x`` is outside this fundamental domain, then the following period shift formulas can be applied until the integral formula can be applied.
We have
```math
    S_2(x,\boldsymbol{ω}) = S_2(x + ω_1,\boldsymbol{ω}) / \bigl(2 \sin(\pi x / ω_2) \bigr)
```
and similarly with ``ω_1`` and ``ω_2`` exchanged.
Note that ``S_2(x,\boldsymbol{ω})`` is symmetric with respect to ``ω_1, ω_2``.
In our implementation, we recursively shift by ``ω_2`` until the fundamental domain is reached.

*Note*: Many authors, including Koyama & Kurokawa, Kurokawa & Koyama, and Tangedal, use a convention which replaces ``S_2`` by ``1/S_2`` relative to our convention.
"""
function double_sine(w::Real, b1::Real, b2::Real; kwargs...)
    @assert b1 * b2 > 0 "Signs of b1 and b2 must agree."

    # ensure both periods are positive
    T = float(promote_type(typeof(w), typeof(b1), typeof(b2)))
    w, b1, b2 = T.(flipsign.((w, b1, b2), b1))

    # Find the smallest shifts that make
    #   w ∈ [Q/4, 3Q/4].
    # Only needs low precision.
    m, n = _best_shift(Float64(w), Float64(b1), Float64(b2))

    w, shiftcorr1, sgn1 = _shift_along(w, b1, b2, m)
    w, shiftcorr2, sgn2 = _shift_along(w, b2, b1, n)

    return exp(sgn1 * sgn2 * (_log_ds(w, b1, b2; kwargs...) + shiftcorr1 + shiftcorr2))
end



@doc """
    _best_shift(z::Float64, x::Float64, y::Float64; eps::Float64=0.25)

Helper function to find integers (m,n) nearly minimizing |m|+|n| subject to
    `|z - m*x - n*y - Q/2| ≤ eps*Q`
where Q = x + y. Note this only needs Float64 precision. This ensures the shifted value
is within the fundamental domain.

Assumes:
 - `x`, `y > 0`;
 - `z` is moderate-sized;
 - `1/2 ≥ eps > 0` is a large constant like 1/4,
 otherwise output is far from optimal and you should use continued fractions. Being away from the boundary is sufficient for the integral to have good convergence properties, so you should never need small `eps`.
"""
function _best_shift(z::Float64, x::Float64, y::Float64; eps::Float64=0.25)
    z /= (x + y)
    p = max(x, y) / (x + y)

    @assert 0 < p < 1
    @assert 0 < eps ≤ 0.5

    best_T = typemax(Int)
    best_m = 0
    best_n = 0

    k0 = round(Int, (z - 0.5) / p)

    # Test the neighborhood.
    # I think we might need radius 3 to avoid numerical edge cases, but I haven't seen any errors.
    for k in (k0-2):(k0+2)
        n = floor(Int, z - k * p)
        m = k + n
        r = z - m * p - n * (1 - p)

        if abs(r - 0.5) ≤ eps
            T = abs(m) + abs(n)
            if T < best_T
                best_T = T
                best_m = m
                best_n = n
            end
        end
    end

    return (x > y ? (best_m, best_n) : (best_n, best_m))
end

@doc """
Apply a shift by `m * δ1` to a double sine function with periods `δ1` and `δ2` and argument `w`.
"""
@inline function _shift_along(w, δ1, δ2, m)
    T = typeof(w)
    v = w - m * δ1
    corr = zero(T)
    sgn = one(T)
    if m ≥ 1
        @inbounds for k in 0:(m-1)
            x = sinpi((v + k * δ1) / δ2)
            sgn *= sign(x)
            corr += log(2 * abs(x))
        end
    elseif m ≤ -1
        @inbounds for k in 1:(-m)
            x = sinpi((v - k * δ1) / δ2)
            sgn *= sign(x)
            corr -= log(2 * abs(x))
        end
    end
    return v, corr, sgn
end


# some helper functions for the double sine integral
_g0(w, b1, b2, t) = sinh(((b1 + b2) / 2 - w) * t) / (2t * sinh(b1 * t / 2) * sinh(b2 * t / 2)) - (b1 + b2 - 2w) / (b1 * b2 * t^2)
_g1(w, b1, b2, t) = exp(-w) / (expm1(-b1 * (t / w + 1)) * expm1(-b2 * (t / w + 1)) * (t + w))
_g(w, b1, b2, t) = _g1(w, b1, b2, t) - _g1(b1 + b2 - w, b1, b2, t)


@doc """
Log of the double sine function inside the fundamental domain.
Uses numerical integration using Gauss-Kronrod quadrature via QuadGK.jl.
"""
function _log_ds(w, b1, b2; kwargs...)
    @assert 0 < w < b1 + b2 "Argument w must be in the fundamental domain."
    T = typeof(w)
    if b1 + b2 ≈ 2w
        return zero(T)
    end

    # integrate from [0,1]
    a = quadgk(t -> _g0(w, b1, b2, t), zero(T), one(T); kwargs...)[1]

    # boundary term
    b = -(b1 + b2 - 2w) / (b1 * b2)

    # integrate the rest, the change of variables makes it from [0,∞).
    c = quadgk(t -> exp(-t) * _g(w, b1, b2, t), zero(T), T(Inf); kwargs...)[1]

    return a + b + c
end
