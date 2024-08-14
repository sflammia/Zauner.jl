export q_pochhammer, q_pochhammer_exp, e, ghost, shin

@doc """
    q_pochhammer(a, q, n)

Finite q-Pochhammer symbol, ``\\prod_{k=0}^{n-1} \\bigl(1-a q^k\\bigr)``.
"""
function q_pochhammer(a, q, n)
    prod(one(a) - a * q^k for k = 0:n-1)
end


@doc raw"""
    q_pochhammer_exp(z, τ, n)

Finite exponential-variant q-Pochhammer symbol, extended to include ``n < 0``.
Defined as
```math
\varpi_n(z,\tau) =
\begin{cases}
        \prod_{j=0}^{n-1}\bigl(1-\mathrm{e}^{2 \pi i (z+j \tau)}\bigr)
        \qquad & n>0
        \\
        1 \qquad & n=0
        \\
        \prod_{j=n}^{-1}\bigl(1-\mathrm{e}^{2 \pi i (z+j \tau)}\bigr)^{-1}
        \qquad & n<0
\end{cases}\,.
```
"""
q_pochhammer_exp(z, tau, n) = (n ≥ 0 ? q_pochhammer(e(z), e(tau), n) : (1 - e(z)) / q_pochhammer(e(z), e(-tau), 1 - n))



@doc raw"""
    e(z)

Normalized exponential function, ``e(z) = \exp(2 \pi i z)``.
"""
function e(z)
    cispi(2 * z)
end



@doc raw"""
    ghost(F:AdmissibleTuple)

Compute a ghost as a d × d matrix from the admissible tuple `F`.

Only rank-1 ghosts are supported at this time.
"""
ghost(F::AdmissibleTuple) = (F.r == 1 ? _rank_1_ghost(F) : _general_ghost(F))

# No support for rank-r ghosts yet.
function _general_ghost(F::AdmissibleTuple)
    error("Only rank-1 ghosts are supported at this time.")
end


# compute the rank-1 ghosts in two cases
# Ideally this should test if F.Q is equivalent to principle
# (using 'is_equivalent( Q1, Q2; proper = true)')
# rather than using only the reduced principle form
_rank_1_ghost(F::AdmissibleTuple) = (F.Q.a == 1 && F.Q.b == 1 - F.d && F.Q.c == 1 ? _principle_ghost(F) : _generic_rank_1_ghost(F))


# symmetrized double sine for principle ghosts
# NOTE: This function averages over the action of [-1 -1; 1 0]
# The original choice of Zauner is [0 -1; 1 -1], and we could use this instead
function _triple_double_sine(p, q, F::AdmissibleTuple)
    d = Int(F.d)
    r = mod(-p - q, d)
    (-1)^(d * (p + q) + p * q + min(d, p + q)) *
    double_sine(1 + (q * F.x - p) / d, F.x, 1) *
    double_sine(1 + (p * F.x - r) / d, F.x, 1) *
    double_sine(1 + (r * F.x - q) / d, F.x, 1)
end


# compute the principle ghost with the simplified algorithm
function _principle_ghost(F::AdmissibleTuple)
    d = Int(F.d)
    dsp = zeros(BigFloat, 2, d)
    dsp[1, 1] = sqrt(d + one(BigFloat))
    k = div(d, 2)
    for p2 = 1:k
        dsp[0+1, p2+1] = _triple_double_sine(0, p2, F)
        dsp[0+1, d-p2+1] = one(BigFloat) / dsp[0+1, p2+1]
    end
    for p2 = 0:d-3
        dsp[1+1, p2+1] = _triple_double_sine(1, p2, F)
    end
    dsp[1+1, end-1] = dsp[1+1, 1+1]
    dsp[1+1, end] = dsp[0+1, 1+1]

    ζ = -cispi(BigFloat(1) / d)
    χ = [ζ^(p * q) for p = 0:1, q = 0:d-1] .* dsp
    χ = ifft(χ, 2)
    χ = circshift(cumprod(χ[2, :] ./ χ[1, :]), 1)
    χ ./ χ[1]
    # This uses a "projective" normalization instead of 2-norm
end

# use special features of the rank-1 case to avoid calculating all nu.
function _generic_rank_1_ghost(F::AdmissibleTuple)
    d = Int(F.d)
    ζ = -cispi(BigFloat(1) / d)
    QQ = QuadBin(F.A[2, 1], F.A[2, 2] - F.A[1, 1], -F.A[1, 2])
    c = e(-BigFloat(rademacher(F.A)) / 24) / sqrt(BigFloat(d + 1))

    ω = _get_periods(F.A, F.x)
    r = ω ./ circshift(ω, -1)

    χ = zeros(Complex{BigFloat}, d, 2)
    χ[1, 1] = 1
    for j = 1:2*d-1
        p = radix(j, [d, d])
        z = (ω[1] * p[2] - ω[2] * p[1]) / d
        m = Int((-F.A[2, 1] * p[1] + (F.A[1, 1] - 1) * p[2]) / d)

        QA = BigFloat(-QQ(p...) / (d * (d - 2)))
        s = (isodd(d) ? 1 : (1 + p[1]) * (1 + p[2]))
        nu = ζ^QA * (-1)^s * c / q_pochhammer_exp((p[2] * F.x - p[1]) / d, F.x, m)
        for i = 1:(length(ω)-2)
            nu *= _sigma_s(z / ω[i+2], r[i+1])
        end

        χ[p[2]+1, p[1]+1] = ζ^(p[2] * p[1]) * real(nu)
    end
    χ = ifft(χ, 1)
    sqrt(abs(χ[1, 1])) * circshift(cumprod(χ[:, 2] ./ χ[:, 1]), 1)
    # to obtain Ghost projector, replace last line with:
    # ψ = sqrt(abs(χ[1,1]))*circshift(cumprod(χ[:,2]./χ[:,1]), 1)
    # ϕ = circshift(reverse(ψ),1)
    # G = ψ*ϕ'/ϕ'ψ
end


function _get_periods(A, β)
    W = psl2word(A)
    n = length(W) - 1

    B = zeros(eltype(A), n + 2, 2)
    B[1:2, 1:2] = A
    for j = 1:n
        B[j+2, :] = [-1 W[j]] * B[j:j+1, :]
    end
    BigFloat.(B) * [β; BigFloat(1)]
end


function _sigma_s(z, β)
    n = floor(Int, -z) + floor(Int, β / 2)
    a = q_pochhammer_exp(z / β, -1 / β, -n)
    b = e((6 * (z + n)^2 + 6 * (1 - β) * (z + n) + β^2 - 3 * β + 1) / (24 * β))
    c = _ds_int_qgk(z + n + 1, β, BigFloat(1), 21)
    a * b * c
end



@doc """
     shin()

 Shintani-Faddeev modular cocycle.
 Not yet implemented.
 """
function shin()
    return nothing
end
