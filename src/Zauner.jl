using Hecke, QuadGK#, LaTeXStrings

### Shintani-Fadeev Cocycles, Double Sine, q-Pochhammer, and Shin 

"""
qPochhammer(a, q, n=Inf[; maxterms = 100])

q-Pochhammer symbol.
"""
function qPochhammer(a, q, n=Inf; maxterms = 100)
    if n < 0; @error "n must be nonnegative."; end
    if n == 0; return 1; end
    infinite = (n == Inf)
    same = (a == q)
    if infinite
        if abs(q) >= 1
            if same && (q == -1 || q == 1)
                return 0
            end
            @error "q-Pochhammer only defined for |q| < 1."
        elseif q == 0
            return 1 - a
        end
    end
    if infinite && same
        # Euler's pentagonal theorem
        terms = [1.0]
        x1 = q
        x2 = q^2
        for k=1:maxterms
            push!(terms, (-1)^k * x1)
            push!(terms, (-1)^k * x2)
            x1 *= q^(3*k+1)
            x2 *= q^(3*k+2)
        end
        return sum(terms)
    else
        factors = [1.0]
        k = 0
        r = 1
        while k < n && k ≤ maxterms
            push!(factors, 1 - a*r)
            r *= q
            k += 1
        end
        return prod(factors)
    end
end


"""
e(z)

Normalized exponential function, e(z) = exp(2*pi*im*z).
"""
e(z) = exp(2*pi*im*z)


DS0(w,b1,b2,t) = (exp(-t*w)-exp(-t*(b1+b2-w)))/((1-exp(-b1*t))*(1-exp(-b2*t))*t)

DS1(w,b1,b2,t) = DS0(w,b1,b2,t) - (b1+b2-2w)/(b1*b2*t^2)

DS2(w,b1,b2) =  quadgk(t->DS1(w,b1,b2,t),0,1).+
                quadgk(t->DS0(w,b1,b2,t),1,Inf).-
                ((b1+b2-2w)/(b1*b2), 0.0)

DoubleSine0(w,b1,b2) = let (x,y) = DS2(w,b1,b2); return exp(-x), y; end

function DoubleSine1(w,b1,b2)
    if real(b1) < 0 && real(b2) < 0
        return DoubleSine1(-w, -b1, -b2)
    elseif real(b1) > real(b2)
        return DoubleSine1(w, b2, b1)
    elseif 0 < real(w) ≤ real(b1)
        return DoubleSine0(w, b1, b2)
    elseif real(w) ≤ 0
        return (2*sin(pi*w*b1/b2))^(-1) .* DoubleSine1(w + b1, b1, b2)
    else
        return 2*sin(pi*(w-b1)*b1/b2) .* DoubleSine1(w - b1, b1, b2)
    end
end

"""
DoubleSine
"""
function DoubleSine(w, b1, b2)
    if b1*b2 == 0 || (imag(b1/b2) == 0 && real(b1/b2) < 0)
        @error "Domain error on b1, b2."
    end
    ϕ = exp((-im/2)*(angle(b1) + angle(b2)))
    return DoubleSine1(ϕ*w, ϕ*b1, ϕ*b2)
end


"""
sfcocycle
"""
function sfcocycle(M, z, τ, NN)
    if M == [1 0; 0 1]; return 1; end
    if M == [-1 0; 0 -1]; @error "M cannot equal -I."; end
    
    a, c, b, d = M
    S = [0 -1; 1 0]
    if c < 0
        sfcocycle(round.(Int,M^-1), z/(c*τ+d), (a*τ+b)/(c*τ+d), NN)^(-1)
    elseif M == S
        sfs(z, τ, NN)
    else
        k = (c == 0 ? 0 : ceil(Int,a/c) )
        δ = [0 1; -1 k]*M # S^-1*T^-k*M
        s, u, t, v = δ
        sfcocycle(S, z/(u*τ+v), (s*τ+t)/(u*τ+v), NN)*sfcocycle(δ, z, τ, NN)
    end
end


sfs(z, τ, NN) = e((τ-3+τ^-1)/24 + (τ-z)*(1-z)/(4*τ) ) * (1-e(z/τ))*DoubleSine(z, τ, 1, NN)


"""
shin
"""
shin(r, A, τ, NN) = qPochhammer(e(r*[τ; 1]), e(τ), r*([1; 0]-A*[1; 0]))*
                    sfcocycle(A, r*A*[τ; 1], τ, NN)


DedekindPart(x) = x-floor(x)-1/2

Dedekinds(a, c) = sum( [ DedekindPart(k/c)*DedekindPart(k*a/c) for k = 1:(abs(c)-1) ] )

function Rademacher(M)
    a, c, b, d = M
    (c == 0 ? b/d : (a + d)/c - 12*sign(c)*Dedekinds(a, c) )
end

Meyer(A) = 3 - Rademacher(A)

ssep(M) = e(Rademacher(M)/24)



### SIC algebraic data

"""
Fundamental discriminant for the SIC base field.
"""
Δ(d::Integer) = fundamental_discriminant((d+1)*(d-3))

"""
Return the base field of the SIC in dimension d as a tuple K, a
where K is the real quadratic field Q(a) and a = √Δ
"""
sic_base_field(d::Integer,r::Integer=1) = quadratic_field(Δ(d))


"""
The class number of the SIC base field
"""
sic_class_number(K::AnticNumberField) = order(class_group(K)[1])
sic_class_number(d::Integer) = order(class_group(sic_base_field(d)[1])[1])

