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
    return dsShift( z, w1, w2, points)
end

function dsShift(z,w1,w2,points)
    if z <= 0
        return dsShift(z+w2, w1, w2, points)/(2*sin((pi*z)/w1))
    elseif z >= w1+w2
        return dsShift(z-w2, w1, w2, points)*(2*sin((pi*(z-w2))/w1))
    else
        return dsIntQGK(z, w1, w2, points)
    end
end


# some helper functions for the double sine integral
_g0(w,b1,b2,t) = sinh.(((b1+b2)/2-w)*t)./(2t.*sinh.(b1*t/2).*sinh.(b2*t/2)).-(b1+b2-2w)./(b1*b2*t.^2)
_g1(w,b1,b2,t) = exp.(-w)./(expm1.(-b1.*(t./w .+ 1)).*expm1.(-b2.*(t./w.+1)).*(t.+w))
_g(w,b1,b2,t)  = _g1(w,b1,b2,t) .- _g1(b1+b2-w,b1,b2,t)


function dsIntQGK(w, b1, b2, pts)
    if b1 + b2 ≈ 2w; return one(BigFloat); end

    # integrate from [0,1]
    a = quadgk(t -> _g0(w,b1,b2,t), BigFloat(0), BigFloat(1), order = pts)[1]

    # boundary term
    b = -(b1+b2-2w)/(b1*b2)

    # integrate the rest, the change of variables makes it from [0,∞).
    c = quadgk(t -> exp(-t).*_g(w,b1,b2,t), BigFloat(0), BigFloat(Inf), order = pts)[1]

    exp(a+b+c)
end





# Experimenting with directly enforcing the Zauner symmetry...


# _g0(w,x,t) = sinh.(((1+x)/2-w)*t)./(2t.*sinh.(t/2).*sinh.(x*t/2)).-(1+x-2w)./(x*t.^2)
# _g1(w,x,t) = exp.(-w)./(expm1.(-(t./w .+ 1)).*expm1.(-x.*(t./w.+1)).*(t.+w))
# _g(w,x,t)  = _g1(w,x,t) .- _g1(1+x-w,x,t)


# @doc raw"""
#     zine(p,d;points=21)

# Zauner-symmetrized double sine function for the principal form. Numerical integration is done via adaptive Gauss quadrature with optional keyword argument `points` set to 21 by default.
# """
# function zine( p, d; points=21)
#     @assert p[1] ≥ 0 && p[2] ≥ 0 && p[1] < d && p[2] < d && d > 3 "p should be in 0:d-1"

#     if p[1] == 0 && p[2] == 0; return one(BigFloat); end

#     p1, p2 = p
#     p0 = mod(-p1-p2,d)
#     x =  (d-1 + sqrt(BigInt((d-3)*(d+1))))/2

#     w0 = 1 + (p2*x-p1)/d
#     w1 = 1 + (p0*x-p2)/d
#     w2 = 1 + (p1*x-p0)/d
#     W = (w0+w1+w2)/3

#     a = quadgk(t -> _g0(w0,x,t)+_g0(w1,x,t)+_g0(w2,x,t), BigFloat(0), BigFloat(1), order = points)[1]

#     # boundary term
#     b = -3(1+x-2W)/x

#     # integrate the rest
#     c = quadgk(t -> exp(-t).*(_g(w0,x,t)+_g(w1,x,t)+_g(w2,x,t)), BigFloat(0), BigFloat(Inf), order = points)[1]

#     exp(a+b+c)
# end






# Experimenting with other numerical integration routines...



# _g12(w,b1,b2,t) = exp(((b1+b2-w)*(t*(1-w)+w))/((1-t)*w))/(expm1(b1+(b1*t)/(w-t*w))*expm1(b2+(b2*t)/(w-t*w))*(1-t)*(t*(1-w)+w))
# _g2(w,b1,b2,t) = _g12(w,b1,b2,t) .- _g12(b1+b2-w,b1,b2,t)



# do the integral with QuadGK.jl, and count the number of function evals.
# function dsIntQGK(w, b1, b2, pts)
#     # integrate from [0,1]
#     a = quadgk_count(t -> _g0(w,b1,b2,t), BigFloat(0), BigFloat(1), order = pts)

#     # boundary term
#     b = -(b1+b2-2w)/(b1*b2)

#     # integrate the rest
#     t = @elapsed c = quadgk_count(t -> exp(-t).*_g(w,b1,b2,t), BigFloat(0), BigFloat(Inf), order = pts)
#     # println("  ",c[3]," evaluations for STANDARD, ",c[1])

#     t2 = @elapsed c2 = quadgk_count(t -> _g2(w,b1,b2,t), BigFloat(0), BigFloat(1), order = pts)
#     # println("  ",c2[3]," evaluations for t-TRANS., ",c[1])

#     # println("t transform with QuadGK")
#     # @time c4 = quadgk_count(t -> 0.5.*sin(t).*_g(w,b1,b2,-log.((1 .+cos.(t))./2)), BigFloat(0), BigFloat(1), order = pts)
#     # println("  ",c4[3]," evaluations")


#     exp(a[1]+b+c[1])


# end



# do the integral with GaussQuadrature.jl
# slower, but if we do it multiple times then we can reuse the nodes and weights
# function dsIntGQ(w, b1, b2, pts)
#     # compute the nodes and weights for the [0,1] piece
#     t0, w0 = legendre(BigFloat, pts)

#     # integrate using gauss-legendre quadrature
#     a = BigFloat(1)/2 * w0' * _g0.(w,b1,b2,(t0.+1)./2)

#     # boundary term
#     b = -(b1+b2-2w)/(b1*b2)

#     # compute the nodes and weights for the [1,∞) piece
#     t1, w1 = laguerre(pts, BigFloat(0))

#     # integrate using gauss-laguerre quadrature
#     c = w1' * _g(w,b1,b2,t1)

#     exp(a+b+c)
# end


# function dsIntGQ2(w, b1, b2, t0, w0, t1, w1)
#     # integrate using gauss-legendre quadrature
#     a = BigFloat(1)/2 * w0' * _g0.(w,b1,b2,(t0.+1)./2)
#     # boundary term
#     b = -(b1+b2-2w)/(b1*b2)
#     # integrate using gauss-laguerre quadrature
#     c = w1' * _g(w,b1,b2,t1)
#     exp(a+b+c)
# end
