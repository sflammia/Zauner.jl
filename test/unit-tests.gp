\\ Double sine tests
{
my(N = 40, t = (5+sqrt(21))/2 , b1, b2, L = List() );
localprec(N+20);

\\ Some known values from K & K
listput(~L, ds(1,4,6) - 1         );
listput(~L, ds(2,4,6) - sqrt(2)   );
listput(~L, ds(3,4,6) - sqrt(2)   );
listput(~L, ds(4,4,6) - sqrt(3/2) );
listput(~L, ds(5,4,6) - 1         );
listput(~L, ds(6,4,6) - sqrt(2/3) );
listput(~L, ds(7,4,6) - 1/sqrt(2) );
listput(~L, ds(8,4,6) - 1/sqrt(2) );
listput(~L, ds(9,4,6) - 1         );

\\ Some generic identities from K & K
b1 = random(1.0);
b2 = random(1.0);

listput(~L, ds( (b1+b2)/2, b1, b2) - 1             );
listput(~L, ds( b1/2, b1, b2) - sqrt(2)            );
listput(~L, ds( b2/2, b1, b2) - sqrt(2)            );
listput(~L, ds( b1 + b2/2, b1, b2) - sqrt(2)/2     );
listput(~L, ds( b2 + b1/2, b1, b2) - sqrt(2)/2     );
listput(~L, ds( b1, b1, b2) - sqrt(b2/b1)          );
listput(~L, ds( b2, b1, b2) - sqrt(b1/b2)          );

\\ from Shintani
listput(~L, ds(1/2,1,t) - sqrt(2) );
listput(~L, ds(t/2,1,t) - sqrt(2) );
listput(~L, ds(1/3,1,t)*ds(1+t/3,1,t)*ds((2+2*t)/3,1,t) - sqrt((1+sqrt(21))/4 - sqrt((3+sqrt(21))/2)/2) );

\\ from K & W
listput(~L, ds(2-I,1,I) - (-1)^(1/4) );
listput(~L, ds(2-sqrt(2),1,sqrt(2)) + 2^(5/4) * cos(Pi/sqrt(2)) );
L = apply( x->abs(x) < 10^-N, L)
}