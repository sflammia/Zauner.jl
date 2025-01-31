# Tangedal double sine products with lambda phase
# Removes the need for any finite q-Pochhammer symbols in the computation of the ghost overlaps

# I plan to make some of these functions internal after some testing.
export lambda_kopp, gamma_kopp, u_tangedal, shin_of_tuple

function _sympt(p::Vector, x)
    p[2] * x - p[1]
end

# Function $\lambda_(\mathbf{p}/d)(A)$ describing the "$\mathbf{r}$-dependent phase factor" as defined in Proposition 7.20 of Kopp, "The Shintani--Faddeev modular cocycle: Stark units from $q$-Pochhammer ratios"
# This is known to be a rational number but is not currently computed as such. May modify later.
function lambda_kopp(t::AdmissibleTuple, p::Vector)
    hjcd = hj_cycle_data(t, p)
    rhos = hjcd[1]
    ps = hjcd[2]
    n = length(rhos)
    ws = typeof(t.x)[]
    for j = 1:n
        push!(ws, _sympt(ps[j], rhos[j]) / BigInt(t.d))
    end
    sum = typeof(t.x)(0)
    for j = 1:n
        sum += (rhos[j] - ws[j]) * (1 - ws[j]) / rhos[j]
    end
    sum
end

# Function $\gamma(A)$ describing the "global phase factor" as defined in Proposition 7.20 of Kopp, "The Shintani--Faddeev modular cocycle: Stark units from $q$-Pochhammer ratios"
# Should be identical to Rademacher function. Implemented for testing purposes.
function gamma_kopp(t::AdmissibleTuple)
    hjcd = hj_cycle_data(t, [BigInt(0); BigInt(0)])
    rhos = hjcd[1]
    n = length(rhos)
    sum = typeof(t.x)(0)
    for j = 1:n
        sum += rhos[j] - 3 + 1 / rhos[j]
    end
    sum
end

# Stark--Tangedal--Yamamoto invariant as defined in Definition 7.13 of Kopp, "The Shintani--Faddeev modular cocycle: Stark units from $q$-Pochhammer ratios"
function u_tangedal(t::AdmissibleTuple, p::Vector)
    hjcd = hj_cycle_data(t, p)
    rhos = hjcd[1]
    ps = hjcd[2]
    n = length(rhos)
    ws = typeof(t.x)[]
    for j = 1:n
        push!(ws, _sympt(ps[j], rhos[j]) / BigInt(t.d))
    end
    prod = typeof(t.x)(1)
    for j = 1:n
        prod *= double_sine(ws[j], rhos[j], typeof(t.x)(1))
    end
    prod
end

# Formula for shin is from Proposition 7.20 of Kopp, "The Shintani--Faddeev modular cocycle: Stark units from $q$-Pochhammer ratios"
function shin_of_tuple(t::AdmissibleTuple, p::Vector)
    e(gamma_kopp(t) / typeof(t.x)(24) + lambda_kopp(t, p) / typeof(t.x)(4)) * u_tangedal(t, p)
end
