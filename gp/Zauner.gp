/*
print("\n");
print("---------------------------------------------------------------------");
print("       _/_/_/_/_/                                                    ");
print("            _/      _/_/_/  _/    _/  _/_/_/      _/_/    _/  _/_/   ");
print("         _/      _/    _/  _/    _/  _/    _/  _/_/_/_/  _/_/        ");
print("      _/        _/    _/  _/    _/  _/    _/  _/        _/           ");
print("   _/_/_/_/_/    _/_/_/    _/_/_/  _/    _/    _/_/_/  _/            ");
print("---------------------------------------------------------------------");
print("\n\nZauner.gp: A PARI/GP package for computing the conjectural properties of SICs and related objects.\n\n");
print("To display the help message for a function f in the command line, type ?f. To view a list of all user-defined functions in PARI (not just these) type \\u.\n");
*/


\\ -----------------------------------------------------------------------------
\\ Functions for generating valid quadratic forms with a given discriminant.
\\ -----------------------------------------------------------------------------




addhelp(D0,"D0(d):\n\n Outputs a row vector [Delta0,f] of a fundamental discriminant Delta0 and integer f such that f^2*Delta0 = (d+1)*(d-3).")

D0(d) = coredisc((d+1)*(d-3),1)




addhelp(ghostbasis,"ghostbasis(d,DeltaQ):\n\n Let Delta0 be the fundamental discriminant Delta0 = quaddisc((d+1)*(d-3)). If DeltaQ divides (d+1)*(d-3), then this function generates a basis for all the quadratic forms with discriminant DeltaQ = fQ^2*Delta0. The output is a vector [Q, p] where Q[k] is a binary quadratic form and p[k] is the period of Q[k]. Thus, p is the cycle decomposition of the class group and Q are explicit generators under composition.\n\n Example:\n\n ? d = 23;\n ? [qb,p] = ghostbasis(d,(d+1)*(d-3)) \n%1 [[Qfb(5, 20, -4, 0.E-38), Qfb(7, 16, -8, 0.E-38)], [2, 2]]\n\n These ghosts have Z2 x Z2 class group. We can generate a representative form from each class:\n\n ? [qb[1]^0*qb[2]^0, qb[1]^0*qb[2]^1, qb[1]^1*qb[2]^0, qb[1]^1*qb[2]^1]\n %2 [Qfb(1, 20, -20, 0.E-38), Qfb(7, 16, -8, 0.E-38), Qfb(5, 20, -4, 0.E-38), Qfb(3, 18, -13, 0.92936294776622502376290462060643844164)]")

ghostbasis(d,DeltaQ) = 
{
    my( dd = (d+1)*(d-3), ff, rem, v);
    if( dd % DeltaQ != 0, error("DeltaQ must divide (d+1)*(d-3)."));
    
    [ff,rem] = divrem(DeltaQ, D0(d)[1]);
    if( rem != 0, error("Ghost fundamental discriminant must divide DeltaQ."));
    if( !issquare(ff), error("DeltaQ must be a square times the ghost fundamental discriminant."));
  
    v = quadclassunit(DeltaQ);
    \\ if the class number is 1, insert an explicit generator.
    if(v[1]==1, 
        v[2] = [1];
        v[3] = [prinredqfb(DeltaQ)];
    );
    [v[3],v[2],v[1]]
}




addhelp(prinredqfb,"prinredqfb(D):\\ A reduced principal form for discriminant D. Returns an error if D is a square or is not 0, 1 mod 4.")

prinredqfb(D) = 
{
    my(s, b, r = D%4); 
    if( issquare(D) || r >= 2 , error("D must be a valid discriminant.") );
    if( D > 0,
            s = sqrtint(D);
            b = s - bitxor( r%2, s%2 );,
        D < 0,
            b = r%2;
    );
    return( Qfb(1,b,(b^2-D)/4) )    
}




addhelp(e, "e(z):\n\n The normalized exponential function exp(2*Pi*I*z).")

e(z) = exp(2*Pi*I*z)




addhelp(WH, "WH(m,n,d):\n\n The Weyl-Heisenberg displacement operator with p[1]=m, p[2]=n in d dimensions using Appleby's 2005 convention for the polarization.")

WH(m,n,d) = matrix(d,d,j,k, (-e(1/(2*d)))^(m*n)*((j-k) % d == m)*e((k-1)*n/d) )




addhelp(radix,"radix(n,r):\n\n The #r least significant digits of the integer n in mixed radix r = [r1,r2,...,rk].\n\n Example:\n\n radix(n,[2,2,2]) computes the three least significant binary digits of n. \n\n Example: \n\n radix(n,[365,24,60,60]) converts seconds to days, hours, minutes, seconds, modulo one year. ")

radix(n,r) = vector(#r, k, (n\vecprod(r[k+1..#r])) % r[k] )




addhelp(hjcontfrac,"hjcontfrac(x,k):\n\n The first k integers in the Hirzebruch-Jung (or negative regular) continued fraction expansion of x.")

hjcontfrac(x,k) = -contfrac(x, vector(k,j,-1));




addhelp(stabilizer,"stabilizer(Q):\n\n Compute a generator L of the stability group of the integral binary quadratic form Q. This is chosen to have det(L) = 1, and the full stabilizer group (including torsion) is given by elements of the form +/- L^k for integer k.\n\n Example:\n\n ? Q = Qfb(1, 3, -1); \n\n ? L = stabilizer(Q); \n\n ? L~*Mat(Q)*L == Mat(Q) && matdet(L) == 1 \n\n %1 1")

stabilizer(Q) = 
{
    my(DeltaQ = -4*matdet(Mat(Q)), u, x, y, S = [0,-1;1,0]);
    u = quadunit(DeltaQ);

    if( norm(u) == -1, u = u^2);
    
    x = trace(u)/2;
    y = imag(u);
    
    x * [1, 0; 0, 1] + y * S * Mat(Q)
}




addhelp(towerh,"towerh(d): Compute the height h in the dimension tower for d_h (usually h is called r).")

towerh(d) = 
{
    my(Delta0 = D0(d)[1], u);
    u = quadunit(Delta0);
    if( norm(u) == -1, u = u^2);
    round(acosh((d-1)/2)/log(u))
}




addhelp(stabilizersl2order,"stabilizersl2order(L,d):\n\n Compute the order of the SL(2,Z) stabilizer matrix L modulo d. No check is done that the matrix is indeed a valid stabilizer, so if the order is not a multiple of 3 (as it must be for valid input) then you may encounter an infinite loop.")
\\ This is a brute force algorithm and surely a much faster approach exists. 
\\ But since k ≤ 11 for d ≤ 10^6, this suffices in practice.

stabilizersl2order(L,d) = 
{
    if( issl2z(L), ,"First input is not an element of SL(2,Z).");
    my(k=1, A, M = Mod(L*L*L,d));
    A = M;
    while(A!=matid(2),
        A = A*M;
        k = k+1;
    );
    3*k
}




addhelp(issl2z,"issl2z(M):\n\n Test if the input M is an element of SL(2,Z).")

issl2z(A) = ( type(A) == "t_MAT" ) && ( matdet(A) == 1 ) && ( length(A) == 2 ) && ( type(A[1,1]) == "t_INT" ) && ( type(A[1,2]) == "t_INT" ) && ( type(A[2,1]) == "t_INT" ) && ( type(A[2,2]) == "t_INT" )




addhelp(psl2word,"psl2word(A):\n\n Decompose a matrix M in SL(2,Z) into a product of S and T generators, modulo -I. Reduces using rounding up with ceiling (Hirzebruch-Jung or negative regular continued fraction reduction) and returning a product strictly in terms of S and T except for the first or final element, which might be negative.\n\n If a vector v is input, then return the product T^v[1]*S*T^v[2]*S...*S*T^v[#v]. ")

psl2word(A) = 
{
    \\ If A is a vector, multiply T^A[1]*S*T^A[2]*S...*S*T^A[#A].
    if( type(A) == "t_VEC", 
        return( fold( (a,b)->a*b, vector(#A, k, [A[k],-1;1,0])) * [0,1;-1,0] )
    );
    
    \\ Otherwise, input must be in SL2(Z).
    if( !issl2z(A), error("A matrix input should be in SL(2,Z)."));
    
    \\ reduce
    my(B = A, a, b, c, d, w = List(), n);
    while(abs(B) != [1, 0; 0, 1],
        a = B[1,1]; b = B[1,2]; c = B[2,1]; d = B[2,2];
        if( c == 0,
            listput(w, sign(a)*b);
            return( Vec(w) );
            ,
            n = max(0,ceil(a/c)); \\ `round`, or `ceil` for HJ reduction
            listput(w, n);
            B = [c, d; n*c - a, n*d - b];
        );
    );
    listput(w,0);
    return( Vec(w) );
}




addhelp(rademacher,"rademacher(A):\n\n Compute the integer Rademacher invariant of an element of SL(2,Z). This is a class invariant function.")
\\ For hyperbolic elements, Meyer(A) = -1/3*Rademacher(A).
rademacher(A) = 
{
    if( !issl2z(A), error("Input should be in SL(2,Z)."));
    
    my(a = A[1,1], b = A[1,2], c = A[2,1], d = A[2,2], t=trace(A));
    
    if(c == 0,
            b/d
        ,
            t/c - 12*sign(c)*sumdedekind(a,abs(c)) - 3*sign(c*t)
    )
}




addhelp(ds,"Double sine function, ds(w,b1,b2).")
ds(w,b1,b2) = 
{
    my(z,w1,w2);
    
    if( b1*b2 == 0 , 
        error("Domain error, b1*b2 == 0.")
    );
    if( (imag(b1/b2) == 0) && (real(b1/b2) < 0) ,
        error("Domain error, imag(b1/b2)==0 and real(b1/b2) < 0.");
    );
    
    \\ Rotate so that w1 and w2 have arg(w1) == -arg(w2)
    [z,w1,w2] = exp((-I/2)*(arg(b1) + arg(b2)))*[w,b1,b2];
    
    \\ Real parts of w1 and w2 have the same sign now, so flip arguments if necessary.
    if( real(w1) < 0, 
        [z,w1,w2] = [-z,-w1,-w2];
    );
    
    \\ Now real(w1) > 0 and real(w2) > 0
    return( dsShift(z,w1,w2) );
}


\\ Shift real(z) until 0 < real(z) < real(w1+w2)
dsShift(z,w1,w2) =
{
    if(real(z) <= 0, 
        return( 2*sin(Pi*z/w1)*dsShift(z+w2, w1, w2) ),
    real(z) >= real(w1+w2),
        return( dsShift(z-w2, w1, w2)/(2*sin(Pi*(z-w2)/w1)) ),
    1,
        return( dsInt(z, w1, w2) )
    );
}


\\ Integral formula for the double sine function
dsInt(w,b1,b2) =
{    
    my( a, b, c, d, EPS);

    EPS = 1/10; \\ hard code the boundary for now. In principle, this choice shouldn't matter.
    
    \\ integrate the power series around zero, up to EPS
    a = dsPowerSeries(w,b1,b2,EPS);
    
    \\ boundary term
    b = - (b1+b2-2*w)/(b1*b2*EPS); 

    \\ integrate EPS to oo
    c = intnum(t=EPS, [+oo, real(w)], \
            exp(-t*w)/(t*(expm1(-b1*t))*(expm1(-b2*t))) ); 

    d = intnum(t=EPS, [+oo, real(b1+b2-w)], \
            -exp(-t*(b1+b2-w))/(t*(expm1(-b1*t))*(expm1(-b2*t))) ); 

    exp(+(a+b+c+d))
}



\\ Integrate the Double sine power series in the main domain from zero to eps.
dsPowerSeries(w,b1,b2,eps) = 
{
    my( isp = default(seriesprecision), csp = default(seriesprecision), sa, sb, a = 1, b = 0);
    
    \\ keep doubling seriesprecision until convergence
    while( abs(a-b) > 0,
        my( t='t);
        default(seriesprecision, csp);
        \\ integrate up to csp terms.
        sa = intformal(serchop(Ser( sinh(((b1+b2)/2-w)*t)/ \
            (2*t*sinh(b1*t/2)*sinh(b2*t/2)) - (b1+b2-2*w)/(b1*b2*t^2), t)));
        \\ explicitly enforce odd parity
        sa = (sa - subst(sa,t,-t))/2;
        \\ for comparison, truncate to order csp-4.
        sb = Ser(sa,t,csp-4);
        a = subst(truncate(sa),t,eps);
        b = subst(truncate(sb),t,eps);
        csp = 2*csp;
    );
    \\restore original seriesprecision and return a
    default(seriesprecision, isp);
    a
}



\\ Compute the nodes and weights for n-point Gauss-Legendre quadrature
gausslegendre(n) = 
{
    my( t, ww, w);

    t = polrootsreal(pollegendre(n));
    ww = vector(#t, k, pollegendre(n,t[k],1));
    w = vector(#t, k , 2*(1-t[k]^2)/(n^2*(ww[k][1] - t[k]*ww[k][2])^2) );
    [t,w]
}

\\ Compute the nodes and weights for n-point Gauss-Laguerre quadrature
gausslaguerre(n) = 
{
    my( t, w);

    t = polrootsreal(pollaguerre(n));
    w = vector(#t, k , t[k]/((n+1)^2*(pollaguerre(n+1,0,t[k]))^2) );
    [t,w]
}




addhelp(qp,"qp(a,q,n): the finite q-Pochhammer symbol (a,q)_n.")

qp(a, q, n) = prod(k=0,n-1, 1-a*q^k)




addhelp(qpe,"qpe(z,tau,n): the finite q-Pochhammer symbol with exponential variables e(z) and e(tau), extended for n < 0.")

qpe(z,tau,n) = if( n >= 0, qp(e(z),e(tau),n) , (1-e(z))/qp(e(z),e(-tau),-n+1) )




\\ sigma_S as a factor times the double sine integral formula.
sds(z,beta) = 
{
    my( n = floor(-z), a, b, c );
    
    a = qpe(z/beta,-1/beta,-n);
    b = e((6*(z+n)^2+6*(1-beta)*(z+n)+beta^2-3*beta+1)/(24*beta));
    c = ds(z+n+1,beta,1);
    a*b*c
}




shin(A,d,p,beta) = 
{
    my( W, z, m, B = A, J, zvals = List([]), bvals = List([]), S);
    
    z = ((A^-1*[p[1];p[2]])~ * [0,-1;1,0] * [beta;1])[1,1] / d;
    W = psl2word(A);
    for( j=2, #W,
        B = [0,1;-1,W[j-1]]*B;
        J = (B[2,1]*beta+B[2,2]);
        listput(zvals, z/J );
        listput(bvals, (B[1,1]*beta+B[1,2])/J );
    );
    zvals = Vec(zvals);
    bvals = Vec(bvals);
    m = (-A[2,1]*p[1]+(A[1,1]-1)*p[2])/d;
    
    S = prod(k=1, #zvals, sds(zvals[k],bvals[k]) );
    
    S / qpe( (p[2]*beta-p[1])/d, beta, m ) 
}




nu(A,d,p,beta,q=[0,0]) = 
{
    if( p % d  == [0,0], return(1) );
    
    my( zeta = -e(1/(2*d)), s, QA, f );
    
    QA = Mat(Qfb(A[2,1],A[2,2]-A[1,1],-A[1,2]))/(d*(d-2));
    s = if( d%2 == 1, 1, (1+p[1])*(1+p[2])+q[1]*p[2]-q[2]*p[1]);
    f = zeta^(-(p*QA*p~)) * (-1)^s * e(-rademacher(A)/24) / sqrt(d+1);
    return( real(f * shin(A, d, p, beta)) )
}




\\ from (d,Q), compute (A,beta).
qdata(d,Q) = {
    my(L, A, a, b, c, w);

    L = stabilizer(Q);
    a = stabilizersl2order(L,d);
    A = L^a;

    [a,b,c] = Vec(Q);
    w = quadgen(b^2-4*a*c,'w); \\ quadratic generator for disc_Q
    beta = (-b+2*w-(b^2%4))/(2*a); \\ a positive root of Q
    [A,beta]
}




\\ For a given d, find a (reduced) Q for each allowed equivalence class.
allQ(d) = {
    my( Delta0, f, qb, p, c, r, Q = List());

    [Delta0,f] = D0(d);
    f = divisors(f);
    
    for( j=1, #f,
        [qb, p, c] = ghostbasis(d,f[j]^2*Delta0);
        for( k=0, c-1,  
            r = radix(k,p);
            listput(Q,prod( n=1, #r, qb[n]^(r[n]) ) );
        );
    );
    Vec(Q)
}




\\ create the ghost associated to (d,Q).
ghost(d,Q,q=[0,0]) = 
{
    my(A, beta, go);
    [A,beta] = qdata(d,Q);
    go = matrix(d,d,p1,p2, nu(A,d,[p1-1,p2-1],beta,q) );
    sum(p1=0,d-1, sum(p2=0,d-1, go[p1+1,p2+1]*WH(p1,p2,d) ))/d
}




\\ create the ghost associated to (d,Q).
/*
ghost(d,Q,q=[0,0]) = 
{
    my( Delta0, f, Dq, fq, A, B, J, beta, W);
    my( z=List(), b=List(), m, p1, p2, nu, S);
    my( zeta = -e(1/(2*d)), s, QA, g, h, z0);
    
    [Delta0,f] = D0(d);
    [Dq,fq] = coredisc(-4*matdet(Mat(Q)),1);
    
    if( (Dq == Delta0) && (f % fq == 0), 
        Dq = fq^2*Dq; \\ If disc(Q) is valid, redefine Dq = disc(Q).
        ,
        error("Discriminant of Q should divide (d+1)(d-3).");
    );
    
    [A,beta] = qdata(d,Q);
    W = psl2word(A);
    B = A;

    nu = matrix(d,d);
    nu[1,1] = 1;
    
    QA = Mat(Qfb(A[2,1],A[2,2]-A[1,1],-A[1,2]))/(d*(d-2));
    h = e(-rademacher(A)/24) / (d*sqrt(d+1));
    
    for( j=2, #W,
        B = [0,1;-1,W[j-1]]*B;
        J = (B[2,1]*beta+B[2,2]);
        listput(z, 1/J );
        listput(b, (B[1,1]*beta+B[1,2])/J );
    );
    z = Vec(z);
    b = Vec(b);

    
    for(p=1,d^2-1,
        [p1,p2] = radix(p,[d,d]);
        z0 = ((A^-1*[p1;p2])~ * [0,-1;1,0] * [beta;1])[1,1] / d;
        m = (-A[2,1]*p1+(A[1,1]-1)*p2)/d;
        
        S = prod(k=1, #z, sds(z0*z[k],b[k]) ) / qpe( (p2*beta-p1)/d, beta, m );

        s = if( d%2 == 1, 1, (1+p1)*(1+p2)+q[1]*p2-q[2]*p1);
        g = zeta^(-([p1,p2]*QA*[p1;p2])[1]) * (-1)^s;
        nu[p1+1,p2+1] = real(g * h * S);
    );

    sum(p1=0,d-1, sum(p2=0,d-1, nu[p1+1,p2+1]*WH(p1,p2,d) ))
    
}



vals(A,W,beta,d,p1,p2) = 
{
    my( z, m, B = A, J, zvals = List([]), bvals = List([]), S);
    
    z = ((A^-1*[p1;p2])~ * [0,-1;1,0] * [beta;1])[1,1] / d;
    for( j=2, #W,
        B = [0,1;-1,W[j-1]]*B;
        J = (B[2,1]*beta+B[2,2]);
        listput(zvals, z/J );
        listput(bvals, (B[1,1]*beta+B[1,2])/J );
    );
    zvals = Vec(zvals);
    bvals = Vec(bvals);
    m = (-A[2,1]*p1+(A[1,1]-1)*p2)/d;
    [zvals, bvals, m]
}


shin2(A,d,beta,W,p) = 
{
    my( z, m, B = A, J, zvals = List([]), bvals = List([]), S);
    
    z = ((A^-1*[p[1];p[2]])~ * [0,-1;1,0] * [beta;1])[1,1] / d;
    for( j=2, #W,
        B = [0,1;-1,W[j-1]]*B;
        J = (B[2,1]*beta+B[2,2]);
        listput(zvals, z/J );
        listput(bvals, (B[1,1]*beta+B[1,2])/J );
    );
    zvals = Vec(zvals);
    bvals = Vec(bvals);
    m = (-A[2,1]*p[1]+(A[1,1]-1)*p[2])/d;
    
    S = prod(k=1, #zvals, sds(zvals[k],bvals[k]) );
    
    S / qpe( (p[2]*beta-p[1])/d, beta, m ) 
}

*/




print("Zauner.gp  v0.0.1");
