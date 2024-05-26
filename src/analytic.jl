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

    z = (-(β*B[2,1]+B[2,2])*BigFloat(p[1]) + (β*B[1,1]+B[1,2])*BigFloat(p[2]))/BigFloat(d)
    
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
    S / q_pochhammer_exp( (BigFloat(p[2])*β-BigFloat(p[1]))/d, β, m )
end


@doc raw"""
    nu(F::AdmissibleTuple) → Matrix{BigFloat}
    nu(F::AdmissibleTuple,p::Vector{Integer}) → BigFloat

Calculate all the ghost overlaps, or one specific ghost overlap specified by `p`.  
"""
function nu(F::AdmissibleTuple)
    M = zeros(BigFloat,F.d,F.d)
    M[1,1] = BigFloat(1)
    
    ζ = -e(BigFloat(1)/(2*F.d))
    t = e(-BigFloat(rademacher(F.A))/24)
    for p in [radix(pp,[F.d,F.d]) for pp=1:(F.d^2-1)]
        QA = BigFloat(-F.Q(p...)*(F.r*F.f//conductor(F.Q)))
        s = BigFloat( F.d%2 == 1 ? 1 : (1+p[1])*(1+p[2]) )
        f = ζ^QA * (-1)^s * t
        M[p[1]+1,p[2]+1] = real(f * shin(F.A, F.d, p, F.x))
    end
    M
end

function nu(F::AdmissibleTuple,p)
    if rem.(p,F.d) == [0, 0]
        return(BigFloat(1))
    end
    
    ζ = -e(BigFloat(1)/(2*F.d))
    QA = BigFloat(-F.Q(p...)*(F.r*F.f//conductor(F.Q)))
    s = BigFloat( F.d%2 == 1 ? 1 : (1+p[1])*(1+p[2]) )
    f = ζ^QA * (-1)^s * e(-BigFloat(rademacher(F.A))/24)
    
    real(f * shin(F.A, F.d, p, F.x))
end



@doc raw"""
    ghost(F:AdmissibleTuple) → Matrix{Complex{BigFloat}}

Compute a ghost as a d × d matrix from the admissible tuple `F`. 
"""
ghost(F::AdmissibleTuple) = ( F.r == 1 ? _rank_1_ghost(F) : _general_ghost(F) )

# Compute all d^2 values of nu for computing the ghost
function _general_ghost(F::AdmissibleTuple)
    M = nu(F)/(sqrt(BigFloat(F.n))*BigFloat(F.d))
    G = zeros(Complex{BigFloat},F.d,F.d)
    for p in [radix(pp,[F.d,F.d]) for pp=1:(F.d^2-1)]
        G += M[p[1]+1,p[2]+1].*wh(p,F.d,Complex{BigFloat})
    end
    G += BigFloat(F.r)/BigFloat(F.d)*wh([0,0],F.d,Complex{BigFloat})
end


# compute the rank-1 ghosts in two cases
_rank_1_ghost(F::AdmissibleTuple) = ( F.Q.a == 1 && F.Q.b == 1-F.d && F.Q.c == 1 ? _principle_ghost(F) : _generic_rank_1_ghost(F) )


# symmetrized double sine for principle ghosts
# NOTE: This function averages over the action of [-1 -1; 1 0]
# The original choice of Zauner is [0 -1; 1 -1], and we could use this instead
function _triple_double_sine(p,q,F::AdmissibleTuple)
    r = mod(-p-q,F.d)
    (-1)^( F.d * (p+q) + p*q + min( F.d, p+q) ) * 
        double_sine( 1 + (q*F.x-p)/F.d, F.x, 1) * 
        double_sine( 1 + (p*F.x-r)/F.d, F.x, 1) * 
        double_sine( 1 + (r*F.x-q)/F.d, F.x, 1)
end


# compute the principle ghost with the simplified algorithm
function _principle_ghost(F::AdmissibleTuple)
    dsp = zeros(BigFloat,2,F.d)
    dsp[1,1] = sqrt(F.d+one(BigFloat))
    k = div(F.d,2)
    for p2=1:k
        dsp[0+1,p2+1] = _triple_double_sine( 0, p2, F)
        dsp[0+1,F.d-p2+1] = one(BigFloat)/dsp[0+1,p2+1]
    end
    for p2=0:F.d-3
        dsp[1+1,p2+1] = _triple_double_sine( 1, p2, F)
    end
    dsp[1+1,end-1] = dsp[1+1,1+1]
    dsp[1+1,end] = dsp[0+1,1+1]
    
    ζ = -cispi(BigFloat(1)/F.d)
    χ = [ ζ^(p*q) for p=0:1, q = 0:F.d-1] .* dsp
    χ = ifft(χ,2)
    χ = circshift(cumprod(χ[2,:]./χ[1,:]), 1)
    χ./χ[1]
    # This uses a "projective" normalization instead of 2-norm
end

# use special features of the rank-1 case to avoid calculating all nu.
function _generic_rank_1_ghost(F::AdmissibleTuple)
    ζ = -cispi(BigFloat(1)/F.d)
    QQ = QuadBin(F.A[2,1],F.A[2,2]-F.A[1,1],-F.A[1,2])
    c = e(-BigFloat(rademacher(F.A))/24) / sqrt(BigFloat(F.d+1))

    ω = _get_periods(F.A,F.x)
    r = ω ./ circshift(ω,-1)
    
    χ = zeros(Complex{BigFloat},F.d,2)
    χ[1,1] = 1
    for j = 1:2*F.d-1
        p = radix(j,[F.d,F.d])
        z = (ω[1]*p[2]-ω[2]*p[1])/F.d
        m = Int((-F.A[2,1]*p[1]+(F.A[1,1]-1)*p[2])/F.d)

        QA = BigFloat(-QQ(p...)/(F.d*(F.d-2)))
        s = ( F.d%2 == 1 ? 1 : (1+p[1])*(1+p[2]) )
        nu = ζ^QA * (-1)^s * c / q_pochhammer_exp( (p[2]*F.x-p[1])/F.d, F.x, m )
        for i=1:(length(ω)-2)
            nu *= sigma_S(z/ω[i+2],r[i+1])
        end
    
        χ[p[2]+1,p[1]+1] = ζ^(p[2]*p[1])*real(nu)
    end
    χ = ifft(χ,1)
    sqrt(abs(χ[1,1]))*circshift(cumprod(χ[:,2]./χ[:,1]), 1)
    # to obtain Ghost projector, replace last line with:
    # ψ = sqrt(abs(χ[1,1]))*circshift(cumprod(χ[:,2]./χ[:,1]), 1)
    # ϕ = circshift(reverse(ψ),1)
    # G = ψ*ϕ'/ϕ'ψ
end



function _ourchi(F::AdmissibleTuple)
    ζ = -cispi(BigFloat(1)/F.d)
    QQ = QuadBin(F.A[2,1],F.A[2,2]-F.A[1,1],-F.A[1,2])
    c = e(-BigFloat(rademacher(F.A))/24)

    ω = _get_periods(F.A,F.x)
    r = ω ./ circshift(ω,-1)
    
    # χ = zeros(Complex{BigFloat},F.d,2)
    χ = zeros(Complex{BigFloat},F.d,F.d)
    χ[1,1] = sqrt(BigFloat(F.d+1))
    for j = 1:F.d^2-1 #2*F.d-1
        p = radix(j,[F.d,F.d])
        z = (ω[1]*p[2]-ω[2]*p[1])/F.d
        m = Int((-F.A[2,1]*p[1]+(F.A[1,1]-1)*p[2])/F.d)

        QA = BigFloat(-QQ(p...)/(F.d*(F.d-2)))
        s = ( F.d%2 == 1 ? 1 : (1+p[1])*(1+p[2]) )
        nu = ζ^QA * (-1)^s * c / q_pochhammer_exp( (p[2]*F.x-p[1])/F.d, F.x, m )
        for i=1:(length(ω)-2)
            nu *= sigma_S(z/ω[i+2],r[i+1])
        end
    
        χ[p[1]+1,p[2]+1] = ζ^(p[2]*p[1])*real(nu)
    end
    χ = ifft(χ,2)
    # sqrt(abs(χ[1,1]))*circshift(cumprod(χ[:,2]./χ[:,1]), 1)
    # to obtain Ghost projector, replace last line with:
    # ψ = sqrt(abs(χ[1,1]))*circshift(cumprod(χ[:,2]./χ[:,1]), 1)
    # ϕ = circshift(reverse(ψ),1)
    # G = ψ*ϕ'/ϕ'ψ
end


function _chi(F::AdmissibleTuple)
    ζ = -cispi(BigFloat(1)/F.d)
    QQ = QuadBin(F.A[2,1],F.A[2,2]-F.A[1,1],-F.A[1,2])
    c = e(-BigFloat(rademacher(F.A))/24) / sqrt(BigFloat(F.d+1))

    ω = _get_periods(F.A,F.x)
    r = ω ./ circshift(ω,-1)
    
    # χ = zeros(Complex{BigFloat},F.d,2)
    χ = zeros(Complex{BigFloat},F.d,F.d)
    χ[1,1] = 1
    for j = 1:F.d^2-1 #2*F.d-1
        p = radix(j,[F.d,F.d])
        z = (ω[1]*p[2]-ω[2]*p[1])/F.d
        m = Int((-F.A[2,1]*p[1]+(F.A[1,1]-1)*p[2])/F.d)

        QA = BigFloat(-QQ(p...)/(F.d*(F.d-2)))
        s = ( F.d%2 == 1 ? 1 : (1+p[1])*(1+p[2]) )
        nu = ζ^QA * (-1)^s * c / q_pochhammer_exp( (p[2]*F.x-p[1])/F.d, F.x, m )
        for i=1:(length(ω)-2)
            nu *= sigma_S(z/ω[i+2],r[i+1])
        end
    
        χ[p[2]+1,p[1]+1] = ζ^(p[2]*p[1])*real(nu)
    end
    χ = ifft(χ,1)
    # sqrt(abs(χ[1,1]))*circshift(cumprod(χ[:,2]./χ[:,1]), 1)
    # to obtain Ghost projector, replace last line with:
    # ψ = sqrt(abs(χ[1,1]))*circshift(cumprod(χ[:,2]./χ[:,1]), 1)
    # ϕ = circshift(reverse(ψ),1)
    # G = ψ*ϕ'/ϕ'ψ
end



# Base.@kwarg 


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
