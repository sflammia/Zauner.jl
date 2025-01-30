# Basic 2 x 2 matrices and operations
# Note that matrices can be added (A+B), multiplied (A*B), and raised to positive integer powers (A^Int64(n))

export matrix_I, matrix_P, matrix_S, matrix_T, sl2z_inverse, sl2z_act

function matrix_I()
    [ZZ(1) ZZ(0); ZZ(0) ZZ(1)]
end

function matrix_P()
    [-ZZ(1) ZZ(0); ZZ(0) -ZZ(1)]
end

function matrix_S()
    [ZZ(0) -ZZ(1); ZZ(1) ZZ(0)]
end

function matrix_T()
    [ZZ(1) ZZ(1); ZZ(0) ZZ(1)]
end

function sl2z_inverse(A::Matrix)
    [A[2,2] -A[1,2]; -A[2,1] A[1,1]]
end

function sl2z_act(A::Matrix, x)
    (BigInt(A[1,1])*x+BigInt(A[1,2]))/(BigInt(A[2,1])*x+BigInt(A[2,2]))
end
