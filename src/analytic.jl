export q_pochhammer, q_pochhammer_exp, e, nu, ghost

@doc """
q_pochhammer(a, q, n)

Finite q-Pochhammer symbol.
"""
q_pochhammer(a, q, n) = prod([one(a)-a*q^k for k=0:n-1])
# function q_pochhammer(a, q, n) 
#     x = one(a)
#     for k=0:n-1
#         x *= one(a)-a*q^k
#     end
#     x
# end
    


@doc """
q_pochhammer_exp(z, τ, n)

Finite exponentiated q-Pochhammer symbol, extended to ``n < 0``.
"""
q_pochhammer_exp(z, tau, n) = ( n ≥ 0 ? q_pochhammer(e(z),e(tau),n) : (1-e(z))/q_pochhammer(e(z),e(-tau),1-n) )
# function q_pochhammer_exp(z, tau, n) 
#     println(round(Int,max(n,1-n)))
#     ( n ≥ 0 ? q_pochhammer(e(z),e(tau),n) : (1-e(z))/q_pochhammer(e(z),e(-tau),1-n) )
# end

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


# function ghost(d,Q::QuadBin,q=[0,0])
#     D = discriminant(Q)
    
#     # test compatability of arguments
#     t =  coredisc( (d+1)*(d-3) ) .// coredisc(D)
#     @assert t[1] == 1 "Fundamental discriminant of the quadratic form must equal the fundamental discriminant of Q(√a) where a=(d+1)*(d-3)."
#     @assert isinteger(t[2]) "Conductor of the quadratic form must divide the conductor of Q(√a) where a=(d+1)*(d-3)."
    
#     L = stabilizer(Q)
#     A = L^sl2zorder(L,d)
#     β = (-BigFloat(Q.b)+sqrt(BigFloat(D)))/(2Q.a)
#     ζ = -cispi(BigFloat(1)/d)
#     QQ = QuadBin(A[2,1],A[2,2]-A[1,1],-A[1,2])
#     c = e(-BigFloat(rademacher(A))/24) / sqrt(BigFloat(d+1))
#     W = psl2word(A)
#     n = length(W)-1
    
#     B = zeros(eltype(A),n+2,2)
#     B[1:2,1:2] = A
#     for j=1:n
#         B[j+2,:] = [-1 W[j]]*B[j:j+1,:]
#     end
#     ω = BigFloat.(B)*[β; BigFloat(1)]
#     r = ω ./ circshift(ω,-1)

#     χ = zeros(Complex{BigFloat},d,2)
#     χ[1,1] = 1
#     for j = 1:2*d-1
#         p = radix(j,[d,d])
#         z = (ω[1]*p[2]-ω[2]*p[1])/d
#         m = Int((-A[2,1]*p[1]+(A[1,1]-1)*p[2])/d)

#         QA = BigFloat(-QQ(p...)/(d*(d-2)))
#         s = ( d%2 == 1 ? 1 : (1+p[1])*(1+p[2])+q[1]*p[2]-q[2]*p[1] )
#         nu = ζ^QA * (-1)^s * c / q_pochhammer_exp( (p[2]*β-p[1])/d, β, m )
#         for i=1:n
#             nu *= sigma_S(z/ω[i+2],r[i+1])
#         end
    
#         χ[p[2]+1,p[1]+1] = ζ^(p[2]*p[1])*real(nu)
#     end
#     χ = ifft(χ,1)
#     sqrt(abs(χ[1,1]))*circshift(cumprod(χ[:,2]./χ[:,1]), 1)
#     # to obtain Ghost projector:
#     # ψ = ghost(d,q) 
#     # ϕ = circshift(reverse(ψ),1)
#     # G = ψ*ϕ'/ϕ'ψ
# end


function ghost(d,Q::QuadBin,q=[0,0])
    D = discriminant(Q)
    
    # test compatability of arguments
    t =  coredisc( (d+1)*(d-3) ) .// coredisc(D)
    @req t[1] == 1 "Fundamental discriminant of the quadratic form must equal the fundamental discriminant of Q(√a) where a=(d+1)*(d-3)."
    @req isinteger(t[2]) "Conductor of the quadratic form must divide the conductor of Q(√a) where a=(d+1)*(d-3)."
    
    L = stabilizer(Q)
    A = L^sl2zorder(L,d)
    β = (-BigFloat(Q.b)+sqrt(BigFloat(D)))/(2Q.a)
    ζ = -cispi(BigFloat(1)/d)
    QQ = QuadBin(A[2,1],A[2,2]-A[1,1],-A[1,2])
    c = e(-BigFloat(rademacher(A))/24) / sqrt(BigFloat(d+1))

    ω = _get_periods(A,β)
    r = ω ./ circshift(ω,-1)
    
    χ = zeros(Complex{BigFloat},d,2)
    χ[1,1] = 1
    for j = 1:2*d-1
        p = radix(j,[d,d])
        z = (ω[1]*p[2]-ω[2]*p[1])/d
        m = Int((-A[2,1]*p[1]+(A[1,1]-1)*p[2])/d)

        QA = BigFloat(-QQ(p...)/(d*(d-2)))
        s = ( d%2 == 1 ? 1 : (1+p[1])*(1+p[2])+q[1]*p[2]-q[2]*p[1] )
        nu = ζ^QA * (-1)^s * c / q_pochhammer_exp( (p[2]*β-p[1])/d, β, m )
        for i=1:(length(ω)-2)
            nu *= sigma_S(z/ω[i+2],r[i+1])
        end
    
        χ[p[2]+1,p[1]+1] = ζ^(p[2]*p[1])*real(nu)
    end
    χ = ifft(χ,1)
    sqrt(abs(χ[1,1]))*circshift(cumprod(χ[:,2]./χ[:,1]), 1)
    # to obtain Ghost projector:
    # ψ = ghost(d,q) 
    # ϕ = circshift(reverse(ψ),1)
    # G = ψ*ϕ'/ϕ'ψ
end


function _get_periods(A,β)
    W = psl2word(A)
    n = length(W)-1
    
    B = zeros(eltype(A),n+2,2)
    B[1:2,1:2] = A
    for j=1:n
        B[j+2,:] = [-1 W[j]]*B[j:j+1,:]
    end
    BigFloat.(B)*[β; BigFloat(1)]
end




function sigma_S(z,β)
    n = floor(Int, -z) + floor(Int, β/2)
    a = q_pochhammer_exp( z/β, -1/β, -n)
    b = e( (6*(z+n)^2+6*(1-β)*(z+n)+β^2-3*β+1)/(24*β) )
    c = dsIntQGK( z+n+1, β, BigFloat(1), 21)
    a*b*c
end

