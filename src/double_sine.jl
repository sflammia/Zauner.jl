export double_sine


@doc raw"""
double_sine(w,b1,b2,pts=21)

Real double sine function. Numerical integration is done via adaptive Gauss quadrature with `pts` points, 21 by default. 
"""
function double_sine( w, b1, b2, pts=21)
    if b1 == 0 || b2 == 0
        error("Domain error, b1*b2 == 0.")
    elseif !all(is_real.([w, b1, b2]))
        error("Only real values supported by this function.")
    elseif sign(b1) != sign(b2)
        error("Signs of b1 and b2 must agree.")
    end

    # return double sine with positive periods and BigFloat
    z, w1, w2 = BigFloat.(flipsign.([w,b1,b2],b1))
    return dsShift( z, w1, w2, pts)
end

function dsShift(z,w1,w2,pts)
    if z <= 0
        return 2*sin((pi*z)/w1)*dsShift(z+w2, w1, w2, pts)
    elseif z >= w1+w2
        return dsShift(z-w2, w1, w2, pts)/(2*sin((pi*(z-w2))/w1))
    else
        return dsIntQGK(z, w1, w2, pts)
    end
end

# some helper functions for the double sine integral
_g0(w,b1,b2,t) = sinh.(((b1+b2)/2-w)*t)./(2t.*sinh.(b1*t/2).*sinh.(b2*t/2)).-(b1+b2-2w)./(b1*b2*t.^2)
_g1(w,b1,b2,t) = exp.(-w)./(expm1.(-b1.*(t./w .+ 1)).*expm1.(-b2.*(t./w.+1)).*(t.+w))
_g(w,b1,b2,t)  = _g1(w,b1,b2,t) .- _g1(b1+b2-w,b1,b2,t)


# do the integral with QuadGK.jl, and count the number of function evals.
# function dsIntQGK(w, b1, b2, pts)
#     # integrate from [0,1]
#     a = quadgk_count(t -> _g0(w,b1,b2,t), BigFloat(0), BigFloat(1), order = pts)
        
#     # boundary term
#     b = -(b1+b2-2w)/(b1*b2)

#     # integrate the rest
#     c = quadgk_count(t -> exp(-t).*_g(w,b1,b2,t), BigFloat(0), BigFloat(Inf), order = pts)

#     println((a[3],c[3],w,b1))
#     exp(a[1]+b+c[1])
# end


function dsIntQGK(w, b1, b2, pts)
    # integrate from [0,1]
    a = quadgk(t -> _g0(w,b1,b2,t), BigFloat(0), BigFloat(1), order = pts)[1]

    # boundary term
    b = -(b1+b2-2w)/(b1*b2)

    # integrate the rest
    c = quadgk(t -> exp(-t).*_g(w,b1,b2,t), BigFloat(0), BigFloat(Inf), order = pts)[1]

    exp(a+b+c)
end


# do the integral with GaussQuadrature.jl
# slower, but if we do it multiple times then we can reuse the nodes and weights
function dsIntGQ(w, b1, b2, pts)    
    # compute the nodes and weights for the [0,1] piece
    t0, w0 = legendre(BigFloat, pts)
    
    # integrate using gauss-legendre quadrature
    a = BigFloat(1)/2 * w0' * _g0.(w,b1,b2,(t0.+1)./2)
        
    # boundary term
    b = -(b1+b2-2w)/(b1*b2)

    # compute the nodes and weights for the [1,∞) piece
    t1, w1 = laguerre(pts, BigFloat(0))
    
    # integrate using gauss-laguerre quadrature
    c = w1' * _g(w,b1,b2,t1)

    exp(a+b+c)
end


function dsIntGQ2(w, b1, b2, t0, w0, t1, w1)
    # integrate using gauss-legendre quadrature
    a = BigFloat(1)/2 * w0' * _g0.(w,b1,b2,(t0.+1)./2)
    # boundary term
    b = -(b1+b2-2w)/(b1*b2)
    # integrate using gauss-laguerre quadrature
    c = w1' * _g(w,b1,b2,t1)
    exp(a+b+c)
end

