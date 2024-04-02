# find an integer relation between the elements of x
function guess_int_null_vec( x::Vector{BigFloat}; warn::Bool = false)
    prec = precision(x[1])
    t = ZZ(2)^prec
    n = length(x)
    # maxnorm = 2maximum(ZZ.(ceil.(x)))*t+1
    v = ZZ.(round.(t .* x))
    L = matrix_space(ZZ,n,n+1)
    T = L([ ZZ(j==k) + (k==n+1)*v[j] for j=1:n, k=1:n+1])
    # B = dropdims(BigInt.(lll_with_removal(T,maxnorm)[2][1,:]),dims=1)[1:n]
    # B = BigInt.(lll_with_removal(T,maxnorm)[2][1,:])[1:n]
    # println(lll(T)[1,1:end-1])
    B = BigInt.(lll(T)[1,:])[1:n]
    
    if warn && -log2( abs( dot(x,B) ) ) - prec + 2n < 0
        @warn "Precision of dot product is low; integer relation may be spurious."
    end
    
    return B
end



# Round into the class field H and then Galois conjugate
function round_conj( F::AdmissibleTuple, V::Vector{BigFloat})
    hb = lll(maximal_order(F.H)).basis_nf # find a good basis
    eH = real_embeddings(F.H)[1] # need to pick a real embedding
    prec = precision(V[1])
    fH = x -> BigFloat.(real.(evaluation_function( eH, prec).(x)))
    primalbasis = fH.(hb)
    dualbasis = fH.(F.g.(hb))
    W = copy(V)

    k = 1
    for v in V
        t = Zauner.guess_int_null_vec( [ primalbasis; v] )
        x = -t[1:end-1]/t[end]
        W[k] = dot( dualbasis, x )
        k += 1
    end
    W
end