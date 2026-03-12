# Tangedal double sine products with lambda phase
# Removes the need for any finite q-Pochhammer symbols in the computation of the ghost overlaps

# I plan to make some of these functions internal after some testing.
export lambda_kopp, gamma_kopp, u_tangedal, shin_of_tuple, shin_rm, sf_phase, rank_1_ghost_var, ghost_proj_var

@doc """
    _sympt(p::Vector, x)

Symplectic inner product. Computes `p[2] * x - p[1]`.
"""
function _sympt(p::Vector, x)
    p[2] * x - p[1]
end


"""
    lambda_kopp(A::Matrix, x::BigFloat, p::Vector, d::Int64)
    lambda_kopp(t::AdmissibleTuple, p::Vector)

Function ``\\lambda_(\\mathbf{p}/d)(A)`` describing the ``\\mathbf{r}``-dependent phase factor as defined in Proposition 7.20 of Kopp, "The Shintani--Faddeev modular cocycle: Stark units from ``q``-Pochhammer ratios".
This is known to be a rational number but is not currently computed as such. May modify later.
"""
function lambda_kopp(A::Matrix, x::BigFloat, p::Vector, d::Int64)
    hjcd = hj_cycle_data(A, x, p, d)
    rhos = hjcd[1]
    ps = hjcd[2]
    n = length(rhos)
    ws = typeof(x)[]
    for j = 1:n
        push!(ws, _sympt(ps[j], rhos[j]) / BigInt(d))
    end
    sum = typeof(x)(0)
    for j = 1:n
        sum += (rhos[j] - ws[j]) * (1 - ws[j]) / rhos[j]
    end
    sum
end
function lambda_kopp(t::AdmissibleTuple, p::Vector)
    lambda_kopp(t.A, t.x, p, t.d)
end


@doc """
    gamma_kopp(A::Matrix, x::BigFloat, d::Int64)
    gamma_kopp(t::AdmissibleTuple)

Function ``\\gamma(A)`` describing the global phase factor as defined in Proposition 7.20 of Kopp, "The Shintani--Faddeev modular cocycle: Stark units from ``q``-Pochhammer ratios".
Should be identical to Rademacher function. Implemented for testing purposes.
"""
function gamma_kopp(A::Matrix, x::BigFloat, d::Int64)
    hjcd = hj_cycle_data(A, x, [BigInt(0); BigInt(0)], d)
    rhos = hjcd[1]
    n = length(rhos)
    sum = typeof(x)(0)
    for j = 1:n
        sum += rhos[j] - 3 + 1 / rhos[j]
    end
    sum
end
function gamma_kopp(t::AdmissibleTuple)
    gamma_kopp(t.A, t.x, t.d)
end


@doc """
    u_tangedal(A::Matrix, x::BigFloat, p::Vector, d::Int64)
    u_tangedal(t::AdmissibleTuple, p::Vector)

Stark--Tangedal--Yamamoto invariant as defined in Definition 7.13 of Kopp, "The Shintani--Faddeev modular cocycle: Stark units from ``q``-Pochhammer ratios".
"""
function u_tangedal(A::Matrix, x::BigFloat, p::Vector, d::Int64)
    hjcd = hj_cycle_data(A, x, p, d)
    rhos = hjcd[1]
    ps = hjcd[2]
    n = length(rhos)
    ws = typeof(x)[]
    for j = 1:n
        push!(ws, _sympt(ps[j], rhos[j]) / BigInt(d))
    end
    prod = typeof(x)(1)
    for j = 1:n
        prod *= double_sine(ws[j], rhos[j], typeof(x)(1))
    end
    prod
end
function u_tangedal(t::AdmissibleTuple, p::Vector)
    u_tangedal(t.A, t.x, p, t.d)
end


"""
    shin_of_tuple(t::AdmissibleTuple, p::Vector)

Formula for shin is from Proposition 7.20 of Kopp, "The Shintani--Faddeev modular cocycle: Stark units from ``q``-Pochhammer ratios".
Valid for admissible tuples.
"""
function shin_of_tuple(t::AdmissibleTuple, p::Vector)
    e(gamma_kopp(t) / typeof(t.x)(24) + lambda_kopp(t, p) / typeof(t.x)(4)) * u_tangedal(t, p)
end

"""
    shin_rm(A::Matrix, x::BigFloat, p::Vector, d::Int64)
    shin_rm(t::AdmissibleTuple, p::Vector)

Formula for shin is from Proposition 7.20 of Kopp, "The Shintani--Faddeev modular cocycle: Stark units from ``q``-Pochhammer ratios".
Valid for powers of real multiplication values.
"""
function shin_rm(A::Matrix, x::BigFloat, p::Vector, d::Int64)
    e(gamma_kopp(A, x, d) / typeof(x)(24) + lambda_kopp(A, x, p, d) / typeof(x)(4)) * u_tangedal(A, x, p, d)
end
function shin_rm(t::AdmissibleTuple, p::Vector)
    shin_of_tuple(t, p)
end


"""
    sf_phase(t::AdmissibleTuple)

DRAFT FUNCTION
The Shintani--Faddeev phase function, from Definition 1.30 of AFK.
ATTENTION: Doesn't preserve precision yet!
"""
function sf_phase(t::AdmissibleTuple, p::Vector)
    d = t.d
    s = (d + (d + 1) * (p[1] + 1) * (p[2] + 1)) % 2
    r = exp(-im * pi * BigFloat(rademacher(t.A)) / 12)
    ξ = -exp(im * pi / d)
    Qp = BigFloat(t.Q.a * p[1]^2 + t.Q.b * p[1] * p[2] + t.Q.c * p[2]^2)
    rjm = BigFloat(t.f)/BigFloat(t.q) # = f_{jm}/f, where f = t.q, for when we generalize to higher rank.

    return (-1)^s * r * ξ^(-rjm * Qp)
end

"""
    rank_1_ghost_var(t::AdmissibleTuple)

DRAFT FUNCTION
Alternative computation of a rank 1 ghost
Passes tests for all SICs in both even and odd dimension d <= 8
"""
function rank_1_ghost_var(t::AdmissibleTuple)
    d = t.d
    ζ = -e(BigFloat(1) / (2 * d))
    χ = zeros(Complex{BigFloat}, d, 2)
    χ[1, 1] = 1
    for j = 1:2*d-1
        p = radix(j, [d, d])
        nu = sf_phase(t,p)*shin_of_tuple(t,p)/sqrt(BigFloat(d+1))
        χ[p[2]+1, p[1]+1] = ζ^(p[2] * p[1]) * real(nu)
    end
    χ = ifft(χ, 1)
    sqrt(abs(χ[1, 1])) * circshift(cumprod(χ[:, 2] ./ χ[:, 1]), 1)
end

"""
    ghost_proj_var(t::AdmissibleTuple)

DRAFT FUNCTION
Alternative computation of ghost projector
Currently just for rank 1
Here for testing purposes
Later will be useful for higher rank, with modifications
"""
function ghost_proj_var(t::AdmissibleTuple)
    d = t.d
    sum = zeros(Complex{BigFloat}, d, d)
    for p = 0:d-1
        for q = 0:d-1
            if [p,q] == [0,0]
                sum += wh(0,0,d)/BigFloat(d) # Multiply by r for higher rank
            else
                sum += sf_phase(t,[p,q])*shin_of_tuple(t,[p,q])*wh(p,q,d)/(d*sqrt(BigFloat(d+1))) # Modify for higher rank
            end
        end
    end
    sum
end
