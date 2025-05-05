# Tangedal double sine products with lambda phase
# Removes the need for any finite q-Pochhammer symbols in the computation of the ghost overlaps

# I plan to make some of these functions internal after some testing.
export lambda_kopp, gamma_kopp, u_tangedal, shin_of_tuple, shin_rm

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
    lambda_kopp(t.A, t.x,p, t.d)
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
