# Should compute the HJ cycle data associated to a pair (t,p)

export hj_cycle_data

# Returns the list of matrices $A_{n,0}$ from Definition 7.4 of Kopp, "The Shintani--Faddeev modular cocycle: Stark units from $q$-Pochhammer ratios"
function _hj_cycle_matrices(t::AdmissibleTuple)
    bs = psl2word(t.A)
    n = length(bs)-1
    if bs[n+1] != 0
        # User should never see this error.
        error("The function _hj_cycle_matrices is not implemented when t.Q is not strictly HJ-reduced. (A strictly HJ-reduced form is one for which the last entry of psl2word(t.A) is zero.)")
    end
    S = matrix_S()
    T = matrix_T()
    Ainvs = typeof(t.A)[matrix_I()]
    for j = 1:(n-1)
        push!(Ainvs, Ainvs[j]*T^BigInt(bs[j])*S)
    end
    As = typeof(t.A)[]
    for j = 1:n
        push!(As, sl2z_inverse(Ainvs[j]))
    end
    As
end

# Shifted vector modulus equal to $d\{\mathbf{p}/d\}$ in the notation of Definition 7.3 of Kopp, "The Shintani--Faddeev modular cocycle: Stark units from $q$-Pochhammer ratios"
function _vec_mod(p::Vector, d)
    [mod(p[1],d)-d; mod(p[2],d)]
end

# Returns the lists of $\beta_n$ (called rhos here) and $d*\mathbf{r}_n$ (called ps here) from Definition 7.4 of Kopp, "The Shintani--Faddeev modular cocycle: Stark units from $q$-Pochhammer ratios"
function hj_cycle_data(t::AdmissibleTuple, p::Vector)
    As = _hj_cycle_matrices(t)
    rhos = typeof(t.x)[]
    ps = typeof(p)[]
    for A in As
        push!(rhos, sl2z_act(A, t.x))
        push!(ps, _vec_mod(A*p, t.d))
    end
    [rhos, ps]
end
