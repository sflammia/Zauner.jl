# Basic 2 x 2 matrices and operations
# Note that matrices can be added (A+B), multiplied (A*B), and raised to positive integer powers (A^Int64(n))

export matrix_I, matrix_P, matrix_S, matrix_T, sl2z_inverse, sl2z_act

const matrix_I = [ZZ(1) ZZ(0); ZZ(0) ZZ(1)]

const matrix_P = [-ZZ(1) ZZ(0); ZZ(0) -ZZ(1)]

const matrix_S = [ZZ(0) -ZZ(1); ZZ(1) ZZ(0)]

const matrix_T = [ZZ(1) ZZ(1); ZZ(0) ZZ(1)]

function sl2z_inverse(A::Matrix)
    [A[2, 2] -A[1, 2]; -A[2, 1] A[1, 1]]
end

function sl2z_act(A::Matrix, x)
    (A[1, 1] * x + A[1, 2]) / (A[2, 1] * x + A[2, 2])
end
