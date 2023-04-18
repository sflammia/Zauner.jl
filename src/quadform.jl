export QuadBin, quadbinid, qmat

# Extend multiplication and whole number exponents to binary quadratic forms
import Base.*, Base.^

*(Q1::QuadBin,Q2::QuadBin) = reduction(compose(Q1,Q2))


function ^(Q::QuadBin,n::Integer)
    if n > 0
        reduction(reduce(*,[Q for k=1:n]))
    elseif n == 0
        quadbinid(discriminant(Q))
    else
        error("negative powers not supported yet for type QuadBin.")
    end
end


# Convert an ideal to the associated quadratic form
function QuadBin(x::NfOrdIdl) 
    x1, x2 = basis(x)
    n1, n2 = norm(x1), norm(x2)
    a, b, c = ZZ.( (n1, n1*x2//x1 + n2*x1//x2, n2) .//norm(x) )
    QuadBin(a,b,c)
end




@doc raw"""
quadbinid(D)

The principal reduced form with discriminant D. 
"""
function quadbinid(D)
    if !is_discriminant(D)
        error("D must be a valid discriminant.")
    end
    
    if D > 0
        s = isqrt(D)
        b = s - (D%2)⊻(s%2)
    else
        b = D%2
    end
    
    QuadBin( 1, b, div(b^2-D,4) )
end


@doc raw"""
qmat(Q::QuadBin)

The matrix of a quadratic form `Q`.
"""
qmat(Q::QuadBin) = [Q.a Q.b//2; Q.b//2 Q.c]

