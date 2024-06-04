export double_sine


@doc raw"""
    double_sine(x,ω1,ω2 [; points=21])

The double sine function, defined in terms of the double gamma function following Shintani's convention:
```math
S_2(x,\boldsymbol{\omega}) = \frac{\Gamma_2(x,\boldsymbol{\omega})}{\Gamma_2(\omega_1+\omega_2-x,\boldsymbol{\omega})}
```

This is currently implemented for real values only.

When ``0 < x < ω_1 + ω_2``, ``\log S_2`` has the integral representation
```math
    \log S_2(x,\boldsymbol{ω}) = \frac12 \int_0^\infty \left(\frac{\sinh\bigl((ω_1 + ω_2 - 2x)t\bigr)}{\sinh(ω_1 t)\sinh(ω_2 t)} - \frac{ω_1 + ω_2 - 2 x}{ω_1 ω_2 t}\right) \frac{\mathrm{d}t}{t}
```
Evaluation proceeds via numerical integration using adaptive Gauss quadrature using the optional keyword argument `points` set to 21 by default.

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
function double_sine( w, b1, b2; points=21)
    @assert b1 != 0 && b2 != 0 "Domain error, b1*b2 == 0."
    @assert all(is_real.([w, b1, b2])) "Only real values supported by this function."
    @assert sign(b1) == sign(b2) "Signs of b1 and b2 must agree."

    # return double sine with positive periods and BigFloat
    z, w1, w2 = BigFloat.(flipsign.([w,b1,b2],b1))
    return _ds_shift( z, w1, w2, points)
end



# domain shift recursively into the fundamental domain
function _ds_shift(z,w1,w2,points)
    if z <= 0
        return _ds_shift(z+w2, w1, w2, points)/(2*sin((pi*z)/w1))
    elseif z >= w1+w2
        return _ds_shift(z-w2, w1, w2, points)*(2*sin((pi*(z-w2))/w1))
    else
        return _ds_int_qgk(z, w1, w2, points)
    end
end



# some helper functions for the double sine integral
_g0(w,b1,b2,t) = sinh.(((b1+b2)/2-w)*t)./(2t.*sinh.(b1*t/2).*sinh.(b2*t/2)).-(b1+b2-2w)./(b1*b2*t.^2)
_g1(w,b1,b2,t) = exp.(-w)./(expm1.(-b1.*(t./w .+ 1)).*expm1.(-b2.*(t./w.+1)).*(t.+w))
_g(w,b1,b2,t)  = _g1(w,b1,b2,t) .- _g1(b1+b2-w,b1,b2,t)



# numerical integration using gauss-kronrod quadrature
function _ds_int_qgk(w, b1, b2, pts)
    if b1 + b2 ≈ 2w; return one(BigFloat); end

    # integrate from [0,1]
    a = quadgk(t -> _g0(w,b1,b2,t), BigFloat(0), BigFloat(1), order = pts)[1]

    # boundary term
    b = -(b1+b2-2w)/(b1*b2)

    # integrate the rest, the change of variables makes it from [0,∞).
    c = quadgk(t -> exp(-t).*_g(w,b1,b2,t), BigFloat(0), BigFloat(Inf), order = pts)[1]

    exp(a+b+c)
end
