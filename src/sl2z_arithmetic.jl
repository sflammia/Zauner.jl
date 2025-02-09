# Basic 2 x 2 matrices and operations
# Note that matrices can be added (A+B), multiplied (A*B), and raised to positive integer powers (A^Int64(n))

export sl2z_I, sl2z_P, sl2z_S, sl2z_T, sl2z_inverse, sl2z_act

@doc """
    sl2z_I

The matrix [1 0; 0 1] with entries of type `ZZRingElem`.
"""
const sl2z_I = [ZZ(1) ZZ(0); ZZ(0) ZZ(1)]

@doc """
    sl2z_P

The matrix [-1 0; 0 -1] with entries of type `ZZRingElem`.
"""
const sl2z_P = [-ZZ(1) ZZ(0); ZZ(0) -ZZ(1)]

@doc """
    sl2z_S

The matrix [0 -1; 1 0] with entries of type `ZZRingElem`.
"""
const sl2z_S = [ZZ(0) -ZZ(1); ZZ(1) ZZ(0)]

@doc """
    sl2z_T

The matrix [1 1; 0 1] with entries of type `ZZRingElem`.
"""
const sl2z_T = [ZZ(1) ZZ(1); ZZ(0) ZZ(1)]


@doc """
    sl2z_inverse(A::Matrix)

Compute the inverse of the ``\\mathrm{SL}_2(\\mathbb{Z})`` matrix `A`.
"""
function sl2z_inverse(A::Matrix)
    [A[2, 2] -A[1, 2]; -A[2, 1] A[1, 1]]
end


@doc """
    sl2z_act(A::Matrix, x)

Action of an ``\\mathrm{SL}_2(\\mathbb{Z})`` matrix `A` as a linear fractional transformation on the number `x`.
"""
function sl2z_act(A::Matrix, x)
    (A[1, 1] * x + A[1, 2] * one(x)) / (A[2, 1] * x + A[2, 2] * one(x))
end
