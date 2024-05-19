export QuadBin, quadbinid, qmat

# Extend multiplication and whole number exponents to binary quadratic forms
import Base.*, Base.^

*(Q1::QuadBin,Q2::QuadBin) = reduction(compose(Q1,Q2))

function ^(Q::QuadBin,n::Integer)
    @req n ≥ 0 "Negative powers not supported yet for type QuadBin."
    n == 0 ? quadbinid(discriminant(Q)) : reduction(reduce(*,[Q for k=1:n]))
end
function ^(Q::QuadBin,n::ZZRingElem)
    @req n ≥ 0 "Negative powers not supported yet for type QuadBin."
    n == 0 ? quadbinid(discriminant(Q)) : reduction(reduce(*,[Q for k=1:n]))
end


# Convert an ideal to the associated quadratic form
function QuadBin(x::AbsSimpleNumFieldOrderIdeal) 
    x1, x2 = basis(x)
    n1, n2 = norm(x1), norm(x2)
    a, b, c = ZZ.( (n1, n1*x2//x1 + n2*x1//x2, n2) .//norm(x) )
    binary_quadratic_form(a,b,c)
end



@doc raw"""
quadbinid(D)

The principal reduced form with discriminant D. 
"""
function quadbinid(D)
    @req is_discriminant(D) "D must be a valid discriminant."
    
    if D > 0
        s = isqrt(D)
        b = s - (D%2)⊻(s%2)
    else
        b = D%2
    end
    
    binary_quadratic_form( one(D), b, div(b^2-D,4) )
end


@doc raw"""
qmat(Q::QuadBin)

The matrix of a quadratic form `Q`.
"""
qmat(Q::QuadBin) = [Q.a Q.b//2; Q.b//2 Q.c]


# Convert a matrix back to a quadratic form.
function QuadBin(M::Matrix{QQFieldElem})
    @req size(M) == (2,2) "Size of M must be (2,2)"
    @req M[2,1] == M[1,2] "M must be symmetric"
    @req isinteger(2M[1,2]) "M[1,2] must be a half-integer"
    @req isinteger(M[1,1]) && isinteger(M[2,2]) "M[1,1] and M[2,2] must be integers"
    binary_quadratic_form(ZZ.([M[1,1],2M[2,1],M[2,2]])...)
end
function QuadBin(M::Matrix{ZZRingElem})
    @req size(M) == (2,2) "Size of M must be (2,2)"
    @req M[2,1] == M[1,2] "M must be symmetric"
    @req isinteger(2M[1,2]) "M[1,2] must be a half-integer"
    @req isinteger(M[1,1]) && isinteger(M[2,2]) "M[1,1] and M[2,2] must be integers"
    binary_quadratic_form(M[1,1],2M[2,1],M[2,2])
end