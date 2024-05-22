
addhelp(D,"D(d): the fundamental discriminant for a SIC in dimension d.")
D(d) = quaddisc((d+1)*(d-3))


addhelp(prinredqfb,"prinredqfb(D): the reduced principal form for discriminant D.")
prinredqfb(D) = 
{
  my(s = sqrtint(D), b);
  b = s - bitxor( D%2, s%2 );
  Qfb(1,b,(b^2-D)/4)
}

addhelp(radix,"radix(n,r): the #r least significant digits of the integer n in radix r.\n\nExample: radix(n,[2,2,2]) computes the three least significant binary digits of n.\n\nExample: radix(n,[365,24,60,60]) converts seconds to days, hours, minutes, seconds, modulo one year.")
radix(n,r) = vector(#r, k, (n\vecprod(r[k+1..#r])) % r[k] )


addhelp(hjcontfrac,"hjcontfrac(x,k): the first k integers in the Hirzebruch-Jung (or negative regular) continued fraction expansion of x.")
hjcontfrac(x,k) = -contfrac(x, vector(k,j,-1));










addhelp(sichilbert,"sichilbert(d): the minimal polynomial for the Hilbert class field over the SIC field, Q(D(d)).")
sichilbert(d) = quadhilbert(D(d))




\\ Lists the SIC class group information for all SICs in dimension d.
\\ Note that we augment class number 1 with an explicit qfb for the class.
\\
sicclass(d) = 
{
  my(D,s,f,ff,Delta,v);
  [D,s] = coredisc( (d+1)*(d-3), 1);
  f = divisors(s);
  ff = apply(sqr,f);
  Delta = ff*D;
  v = vector( #f, k, quadclassunit(Delta[k]) );
  for(k = 1, #v,
    if(v[k][1]==1, 
      v[k][2] = [1];
      v[k][3] = [prinredqfb(Delta[k])];
    );
  );
  v
}


\\ Lists the SIC multiplet numbers in dimension d.
\\ 
sicmult(d) = apply(x->x[1],sicclass(d))


\\{
\\  my(D,DD,f,Delta,s);
\\  D = (d+1)*(d-3);
\\  [DD,s] = coredisc(D,1);
\\  f = divisors(s);
\\  Delta = apply(sqr,f)*DD;
\\  vector( length(f), k, quadclassunit(Delta[k]).no )
\\} 


\\ Total number of SICs in dimension d
\\ 
numsic(d) = 
{
  my(s);
  s = sicmult(d);
  sum( k=1, #s, s[k] )
}


\\ Compute the number of SICs in dimensions n–m inclusive
\\ 
allsicnums(n,m) = parvector(m-n+1, k, numsic(k+n-1) )


/* For a given dimension d, we can compute all of the essential information about the (ghost) SIC that would let us compute the ghost overlaps. 
(1)  d is the dimension,
(2)  D is the fundamental discriminant,
(3)  r is position in tower,
(4)  f is the conductor,
(5)  c is the class number,
(6)  Q is a choice for the primitive HJ reduced form,
(7)  tau is the larger of the two roots of Q,
(8)  HJ is the minimal period of the HJ continued fraction expansion of tau,
(9)  LQ is the generator of the stability group of Q,
(10) n is the order of LQ in SL(2,Z)/Γ(d) (so G=LQ^n) 
?? Also include: Meyer invariant, rmin, anti-unitary symmetry, and total number of SICs, and tau in the integral basis? */
\\ 
sic(d) = 
{
  my(dd, D, s, f, DD, U, m);
  dd = (d+1)*(d-3);
  [D,s] = coredisc(dd,1); \\ fundamental discriminant D, where dd = D*s^2.
\\  db = d*(3+(-1)^d)/2; \\ d-bar
\\  w = quadgen(D); \\ the "canonical" quadratic generator for Q[√D]
\\  H = quadhilbert(D); \\ minimal polynomial of the Hilbert class field
  f = divisors(s);
  DD = apply(sqr,f)*D; \\ discriminants of each order
  U = apply(quadunit,DD); \\ fundamental units for each order
  s = vector( #f, k, quadclassunit(DD[k]) ); \\ all class info about sic multiplets
  m = apply(x -> x[1], s); \\ sic multiplet numbers
}

\\ The Meyer invariant associated to the SL2(Z) matrix A = [a, b; c, d]
\\ Uses the formula 6.17 from Atiyah, "The Logarithm of the Dedekind Eta Function".
\\ ***Currently assumes c ≠ 0 and 𝜖 = -1.***
\\ 
meyer(A) = -(A[1,1]+A[2,2])/(3*A[2,1]) + 4*sign(A[2,1])*sumdedekind(A[1,1],A[2,1]) - (-1)*sign(A[2,1]*(A[1,1]+A[2,2]))



\\ Given a (standard) reduced binary quadratic form, compute the associated set of HJ reduced forms
\\ 
hjSub(q) = 
{
  my(a,b,c,nmax,v);
  [a,b,c] = Vec(q);
  nmax = floor(abs((2*a)/(-b + sqrt(b^2 - 4*a*c))));

  v = vector(nmax,n,[a + sign(a)*n * b + n^2 * c,
   -sign(a)*2*a - (2*n-1) * b - sign(a)*2*(n-1)*n*c, 
   a + sign(a)*(n-1)*b + (n-1)^2*c])
}


\\ Given a (standard) reduced binary quadratic form, compute an associated HJ reduced form and the (negative) Meyer invariant
\\ 
hjSub2(q) = 
{
  my(a,b,c,m,v);
  [a,b,c] = Vec(q);
  m = floor(abs((2*a)/(-b + sqrt(b^2 - 4*a*c))));
  v = Qfb( a, -2*abs(a)-b, a+sign(a)*b+c);
  [m,v]
}



\\ Compute a reduced quadratic form for every class from quadclassunit. 
\\ Here q2 is the list of the orders of q3, the generators. 
\\ 
allQfb(q2,q3) = 
{
  my(L, z, r);
  L = vecprod(q2);
  z = length(q3);
  r = vector(L, k, radix(k-1,q2) );
  vector( L, j, vecprod( vector( z, k, q3[k]^r[j][k] ) ) )
}


/* Let q be a quadratic form with positive root r, and suppose that r has a purely periodic HJ continued fraction expansion c. You could speed this up by iterating each HJ continued fraction expansion step instead of recomputing it each time. */
\\ 
hj(q) = 
{
  my( a, b, c, cf, r, n=1, Q=Mat(q), LQ=[1,0;0,1], p );
  [a,b,c] = Vec(q);
  r = (-b + sqrt(b^2 - 4*a*c) )/(2*a);
  cf = hjcontfrac(r,n);
  LQ = LQ*[cf[n],-1;1,0];
  while( LQ~*Q*LQ != Q , 
    n++; 
    cf = hjcontfrac( r, n ); 
    LQ = LQ*[cf[n],-1;1,0];
  );
  [#cf,cf,LQ]
}


/* Print all (conjectural) information about a SIC in dimension d that appears in Marcus' table.*/
sic(d) = 
{
  my(dd = (d+1)*(d-3), f, s, v = vector(13,k,1));

  v[1] = d; \\ Dimension d
  [DD,f] = coredisc(dd,1);
  v[2] = DD; \\ fundamental discriminant 
  v[3] = 0; \\ position in the dimension tower
  v[4] = divisors(f); \\ list of conductors
  s = sicclass(d); 
  v[5] = apply( x -> x[1], s); \\ class num of order f^2*DD, i.e. size of multiplet.
  v[6] = 0; \\ red. qfb for each class
  v[13] = s;
  v
}


\\ Given a reduced Qfb, find all the 



sic(d) = 
{
  my(s = vector(15,k,0), D, f, ff, Delta, v, h, Q);
  [D,f] = coredisc( (d+1)*(d-3), 1);
  f = divisors(f);
  ff = apply(sqr,f);
  Delta = ff*D;
  v = vector( #f, k, quadclassunit(Delta[k]) );
  for(k = 1, #v,
    if(v[k][1]==1, 
      v[k][2] = [1];
      v[k][3] = [prinredqfb(Delta[k])];
    );
  );
  h = apply(x->x[1],v); \\ class numbers
  Q = apply( x -> allQfb(x[2],x[3]), v); \\ a std. reduced Q for each class
  QQ = Q
  Q = vector(#Q,k,apply(hjSub, Q[k])); \\ all HJ reduced Qs for each class
\\  for(j=1,#h, \
\\    for(k=1,h[j], \
\\      QQ[j][k] = apply(x->hj(call(Qfb,x)),Q[j][k]); \
\\    ); \
\\  );
  s[1]  = d; \\ dimension d
  s[2]  = D; \\ fundamental discriminant Δ
  s[3]  = 0; \\ position in dimension tower r
  s[4]  = f; \\ conductors f
  s[5]  = h; \\ class numbers h
  s[6]  = Q; \\ HJ reduced quadratic forms Q (with minimal period) for each class
  s[7]  = 0; \\ tau, the positive root of each Q
  s[8]  = 0; \\ HJ periodic expansion of tau
  s[9]  = 0; \\ LQ, the SL(2,Z) matrix
  s[10] = 0; \\ n, order of LQ so that G = LQ^n
  s[11] = 0; \\ -Meyer invariant m, (length of HJ)
  s[12] = 0; \\ rmin
  s[13] = 0; \\ anti-unitary? 0 (no) / 1 (yes)
  s[14] = 0; \\ #sics, the number of inequivalent SICs in dimension d
  s[15] = QQ; \\ tau in the integral basis for the order
  s
}




\\ fix the output of quadclassunit so that class number 1 also returns a reduced form


\\------------------------------- Commented for now  -------------------------------\\
/* 
\\ HJ continued fraction expansion to n steps
hjcfstep(x,n) =
{ 
  my( a= vector(n), t );
  x *= -1;
  for (i = 1, n,
    a[i] = floor(x);
    t = x - a[i]; if (!t || i == n, break);
    x = -1 / t;
  ); 
  a; 
}


hjcfstep(x) = 
{ 
  my( a = floor(-x) );
  [1/(x+a), a];
}



\\ Let q be a quadratic form with positive root r, and suppose that r has a purely periodic HJ continued fraction expansion c. 
\\ You could speed this up by iterating each HJ continued fraction expansion step instead of recomputing it each time. 
\\ 
hjtest(q) = 
{
  my( a, b, c, cf, r, n=1, Q=Mat(q), LQ=[1,0;0,1], p );
  [a,b,c] = Vec(q);
  r = (-b + sqrt(b^2 - 4*a*c) )/(2*a);
  cf = hjcfexpand(r,n);
  LQ = LQ*[cf[n],-1;1,0];
  while( LQ~*Q*LQ != Q , 
    n++; 
    cf = hjcfexpand( r, n ); 
    LQ = LQ*[cf[n],-1;1,0];
  );
  [cf,LQ]
}




  
\\ relative equation for the ray class SIC field
\\ sicray(d) = quadray( quaddisc( (d+1)*(d-3) ), [d*(3+(-1)^(d))/2,[2,0]] )


\\ compute the number of SICs in that list of dimensions.
for( d = 4, 101, print("d = "d" , n(d) = "numsic(d)) )

dmax = 2^12+1;
t = vector( dmax-3 );
for( d = 4, dmax, t[d-3] = [d ,numsic(d) ]  )
write("~/Desktop/numsics.txt",t)






dmax=2^13;
t = vector( dmax-3 );
for( d = 4, dmax, t[d-3] = [d ,cn(d) ]  )
write("~/Desktop/cn.txt",t)





d = 6; 

\\ degree of the RCF for the SIC:
rcftest(d) = 
{
  dd = quaddisc((d+1)*(d-3)); 
  db = d*(3+(-1)^d)/2;
  bnf = bnfinit( y^2 - dd );
  bnr = bnrinit( bnf, [db,[1,1]], 1 );
  bnr.clgp[1]
}
*/
