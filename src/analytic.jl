export q_pochhammer, q_pochhammer_exp, e, nu, ghost

@doc """
q_pochhammer(a, q, n)

Finite q-Pochhammer symbol.
"""
q_pochhammer(a, q, n) = prod([one(a)-a*q^k for k=0:n-1])


@doc """
q_pochhammer_exp(z, τ, n)

Finite exponentiated q-Pochhammer symbol, extended to ``n < 0``.
"""
q_pochhammer_exp(z, tau, n) = ( n ≥ 0 ? q_pochhammer(e(z),e(tau),n) : (1-e(z))/q_pochhammer(e(z),e(-tau),1-n) )


@doc raw"""
e(z)

Normalized exponential function, ``e(z) = \exp(2 \pi i z)``.
"""
e(z) = cispi(2*z)



function sds(z,β)
    n = floor( -z) + floor(β/2)
    a = q_pochhammer_exp( z/β, -1/β, -n)
    b = e( (6*(z+n)^2+6*(1-β)*(z+n)+β^2-3*β+1)/(24*β) )
    c = double_sine( z+n+1, β, 1)
    a*b*c
end


function shin(A,d,p,β)
    B = BigFloat.(A)
    
    z = (-(β*B[2,1]+B[2,2])*p[1] + (β*B[1,1]+B[1,2])*p[2])/d
    W = psl2word(A)
    n = length(W)-1
    
    vals = Vector{Tuple{typeof(z),typeof(z)}}(undef,n)
    
    for j=1:n
        B = [0 1; -1 W[j]]*B
        J = (B[2,1]*β+B[2,2])
        vals[j] = ( z/J, (B[1,1]*β+B[1,2])/J )
    end
    m = Int((-A[2,1]*p[1]+(A[1,1]-1)*p[2])/d)
    S = mapreduce( x -> sds(x...), *, vals)
    S / q_pochhammer_exp( (p[2]*β-p[1])/d, β, m )
end


@doc """
  nu((A,d,p,β,q=[0, 0])

Ghost overlaps.
"""
function nu(A,d,p,β,q=[0, 0])
    if rem.(p,d) == [0, 0]
        return(BigFloat(1))
    end
    
    ζ = -e(BigFloat(1)/(2*d))
    QA = BigFloat(-QuadBin(A[2,1],A[2,2]-A[1,1],-A[1,2])(p...)/(d*(d-2)))
    s = BigFloat( d%2 == 1 ? 1 : (1+p[1])*(1+p[2])+q[1]*p[2]-q[2]*p[1] )
    f = ζ^QA * (-1)^s * e(-BigFloat(rademacher(A))/24) / sqrt(BigFloat(d+1))

    real(f * shin(A, d, p, β))
end



@doc """
  ghost((A,d,β,q=[0, 0])

Ghost.
"""
function ghost(A,d,β,q=[0,0])
    GG = [ nu(A,d,[p1,p2],β,q) for p1=0:d-1, p2=0:d-1]
    G = zeros(Complex{BigFloat},d,d)
    for p1=0:d-1, p2=0:d-1
        G += GG[p1+1,p2+1]*WH([p1,p2],d,Complex{BigFloat})
    end
    G /= d
end

# function ghost(d,Q,q=[0,0])
    
#     nothing
# end
    