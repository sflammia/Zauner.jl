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
    b1, b2 = sort((b1, b2)) # b1 < b2 now

    # Find small shifts to get into the fundamental domain, but not too close to the boundary.
    # Only needs low precision.
    n = _minimal_bulk_shift(Float64(w), Float64(b1), Float64(b2))
    w, shift, sgn = _shift_along(w, b2, b1, n)

    return exp(sgn * (_log_ds(w, b1, b2; kwargs...) + shift))
end


@doc """
    _minimal_bulk_shift(w::Float64, b1::Float64, b2::Float64)

Helper function to find an integer `n` that shifts `w` to inside the fundamental domain, but at least `b1/2` away from the boundary.
Assumes `0 < b1 < b2`.
"""
function _minimal_bulk_shift(w::Float64, b1::Float64, b2::Float64)
    Q = b1 + b2
    n_min = ceil((w - (Q - b1 / 2)) / b2)
    n_max = floor((w - b1 / 2) / b2)

    if n_min <= 0 <= n_max
        return 0
    elseif abs(n_min) < abs(n_max)
        return n_min
    else
        return n_max
    end
end


@doc """
Apply a shift by `-m * δ1` to a double sine function with periods `δ1` and `δ2` and argument `w`.
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


struct G0{T}
    w::T
    b1::T
    b2::T
end

@inline (g::G0)(t) = _g0(g.w, g.b1, g.b2, t)

struct Gexp{T}
    w::T
    b1::T
    b2::T
end

@inline (g::Gexp)(t) = exp(-t) * _g(g.w, g.b1, g.b2, t)


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
    g0 = G0(w, b1, b2)
    a = quadgk(g0, zero(T), one(T); kwargs...)[1]

    # boundary term
    b = -(b1 + b2 - 2w) / (b1 * b2)

    # integrate the rest, the change of variables makes it from [0,∞).
    gexp = Gexp(w, b1, b2)
    c = quadgk(gexp, zero(T), T(Inf); kwargs...)[1]

    return a + b + c
end
