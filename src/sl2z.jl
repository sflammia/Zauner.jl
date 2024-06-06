export is_sl2z, sl2zorder, rademacher, psl2word, stabilizer

@doc raw"""
    is_sl2z(M)

Returns `true` if the matrix `M` is a 2x2 integer matrix with determinant 1.
"""
function is_sl2z(M::Matrix{T}) where T<:Integer
    (size(M) == (2,2)) && (M[1,1]*M[2,2]-M[1,2]*M[2,1] == 1)
end
function is_sl2z(M::Matrix{ZZRingElem})
    (size(M) == (2,2)) && (M[1,1]*M[2,2]-M[1,2]*M[2,1] == 1)
end



@doc raw"""
    sl2zorder(L, d)

Find the order of ``L \bmod d`` where ``L`` is in ``\mathrm{SL}(2,\mathbb{Z})``.
Brute force algorithm.
"""
function sl2zorder(L,d)
    @assert is_sl2z(L) error("First input must be in SL(2,Z).")
    @assert d > 1 error("d must be > 1.")

    k=1
    A = L .% d
    while (A - [1 0; 0 1]) .% d != [0 0; 0 0]
        A = L*A .% d
        k += 1
    end
    k
end

@doc raw"""
    rademacher(M)

Rademacher invariant of the ``\mathrm{SL}(2,\mathbb{Z})`` matrix `M`.
Reference:
Hans Rademacher
Zur Theorie der Dedekindschen Summen
Mathematische Zeitschriften vol 63, pp. 445--463 (1955).
"""
function rademacher(M)
    @assert is_sl2z(M) error("Input matrix should be in SL(2,Z).")

    a, c, b, d = M
    t = a + d

    ( c == 0 ? b//d : t//c - 3*sign(c*t) - 12*sign(c)*dedekind_sum(a, abs(c)) )
end



@doc raw"""
    psl2word(v::Vector)
    psl2word(A::Matrix)

Return the product ``T^{v_1} S T^{v_2} S...S T^{v_n}`` where ``S`` and ``T`` are the standard generators of ``\mathrm{SL}(2,\mathbb{Z})``.

Decompose a matrix `A` in ``\mathrm{SL}(2,\mathbb{Z})`` into a product of ``S`` and ``T`` generators, modulo ``-I``.
Reduction is done using rounding up with ceiling (Hirzebruch-Jung or negative regular continued fraction reduction) and returning a product strictly in terms of ``S`` and ``T`` except for the first or final element, which might be negative.
"""
psl2word(A::Vector) = reduce(*, [ [A[k] -1;1 0] for k=1:length(A)] ) * [0 1;-1 0]

function psl2word(A::Matrix)
    @req is_sl2z(A) "Input should be in SL(2,Z)."

    # reduce
    B = ZZ.(A)
    w = ZZRingElem[]
    while abs.(B) != ZZ.([1 0; 0 1])
        if B[2,1] == 0
            push!(w, sign(B[1,1])*B[1,2])
            return(w)
        else
            n = max(ZZ(0),ceil(ZZRingElem,B[1,1]//B[2,1]))::ZZRingElem
            push!(w, n)
            B = [B[2,1]  B[2,2]; n*B[2,1] - B[1,1]  n*B[2,2] - B[1,2]]
        end
    end
    push!(w,ZZ(0))
    w
end



@doc raw"""
    stabilizer(Q::QuadBin [, u::AbsSimpleNumFieldOrderElem])

Compute the stabilizer of `Q` in ``\mathrm{SL}(2,\mathbb{Z})``, that is, compute the matrix ``L`` such that ``L^TQL=Q``, where ``Q`` is the matrix of the quadratic form.
The first input `Q` is a binary quadratic form and the second (optional) input `u` is a fundamental totally positive unit with norm 1 in the ring of integers ``\mathbb{Z}[\omega]``, where ``\omega = \bigl(\Delta\bmod 4 + \sqrt{\Delta}\bigr)/2``.
If the second argument is not specified, the unit is computed.
"""
function stabilizer(Q::QuadBin,u::AbsSimpleNumFieldOrderElem)
    x = trace(u)//2
    y = coordinates(u)[2]
    ZZ.(x.*[1 0; 0 1] + y.*[0 -1; 1 0]*qmat(Q))
end
stabilizer(Q::QuadBin) = stabilizer(Q,pell(Q))
