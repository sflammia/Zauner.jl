# Should compute the HJ cycle data associated to a pair (t,p)

export hj_cycle_data

@doc """
    _hj_cycle_matrices(t::AdmissibleTuple)
    _hj_cycle_matrices(A::Matrix)

Returns the list of matrices ``A_{n,0}`` from Definition 7.4 of Kopp, "The Shintani--Faddeev modular cocycle: Stark units from ``q``-Pochhammer ratios".
"""
function _hj_cycle_matrices(A::Matrix)
    bs = psl2word(A)
    if bs[end] != 0
        # User should never see this error.
        error("The function _hj_cycle_matrices is not implemented when t.Q is not strictly HJ-reduced. (A strictly HJ-reduced form is one for which the last entry of psl2word(t.A) is zero.)")
    end
    _hj_cycle_matrices_from_psl2word(bs)
end
function _hj_cycle_matrices(t::AdmissibleTuple)
    _hj_cycle_matrices(t.A)
end

function _hj_cycle_matrices_from_psl2word(W::Vector)
    # The matrix [0 1; -1 n] is (T^n * S)^(-1), but
    # this avoids powers and inverses and
    # easily ensures type consistency
    U = pushfirst!(map(x -> [0 1; -1 x], W[1:end-2]), eltype(W).([1 0; 0 1]))
    accumulate((A, B) -> B * A, U)
end


@doc """
    _vec_mod(p::Vector, d)

Shifted vector modulus equal to ``d\\{\\mathbf{p}/d\\}`` in the notation of Definition 7.3 of Kopp, "The Shintani--Faddeev modular cocycle: Stark units from ``q``-Pochhammer ratios".
"""
function _vec_mod(p::Vector, d)
    [mod(p[1], d) - d; mod(p[2], d)]
end

@doc """
    hj_cycle_data(t::AdmissibleTuple, p::Vector)
    hj_cycle_data(A:Matrix, x::BigFloat, p::Vector, d::Int64)

Returns the lists of ``\\beta_n`` (called rhos here) and ``d \\mathbf{r}_n`` (called ps here) from Definition 7.4 of Kopp, "The Shintani--Faddeev modular cocycle: Stark units from ``q``-Pochhammer ratios".
"""
function hj_cycle_data(A::Matrix, x::BigFloat, p::Vector, d::Int64)
    As = _hj_cycle_matrices(A)
    rhos = typeof(x)[]
    ps = typeof(p)[]
    for Aj in As
        push!(rhos, sl2z_act(Aj, x))
        push!(ps, _vec_mod(Aj * p, d))
    end
    [rhos, ps]
end
function hj_cycle_data(t::AdmissibleTuple, p::Vector)
    hj_cycle_data(t.A,t.x,p,t.d)
end
