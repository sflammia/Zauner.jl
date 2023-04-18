using Test
using Zauner

#####
#=

TODO list:
 - add tests for `dsShift` part of `double_sine`
 - add tests for `quadclassunit`, or whatever replaces it
=#
####

@testset "Zauner algebraic tests" begin
    
    @test conductor(QuadBin(1,1,-1)) == 1
    @test conductor(QuadBin(2,4,8)) == 4
    
    @test coredisc(5) == (5,1)
    @test coredisc(20) == (5,2)
    @test coredisc(QuadBin(1,1,-1)) == (5,1)
    @test coredisc(QuadBin(2,4,8)) == (-3,4)
    
    x = pell(5)
    @test 1//x+x == trace(x)
    @test norm(x) == 1
    @test (x-1//x)^2 == 5

    y=pell(QuadBin(1,1,-1))
    # note: it if false that x == y (!)
    # presumably this is due to partial type initialization
    @test x.elem_in_nf == y.elem_in_nf
    
    @test towerh.([4,8,19]) == [1,2,3]
    @test towerh.([5,15]) == [1,2]
    
    
    # stabilizer tests for d=4:16
    L =[[3 -1; 1 0], [4 -1; 1 0], [5 -1; 1 0], [7 -2; 4 -1], [6 -1; 1 0], 
        [3 -1; 1 0], [7 -1; 1 0], [8 -1; 1 0], [9 -5; 2 -1], [9 -1; 1 0], 
        [11 -2; 6 -1], [10 -1; 1 0], [11 -4; 3 -1], [13 -3; 9 -2], [11 -1; 1 0],
         [12 -1; 1 0], [13 -7; 2 -1], [13 -1; 1 0], [14 -5; 3 -1], [4 -1; 1 0],
         [15 -2; 8 -1], [14 -1; 1 0], [16 -11; 3 -2], [15 -1; 1 0], [18 -11; 5 -3]]
    q = [1 -3 1; 1 -4 1; 1 -5 1; 2 -4 1; 1 -6 1; 1 -3 1; 1 -7 1; 1 -8 1; 2 -10 5; 
    1 -9 1; 3 -6 1; 1 -10 1; 3 -12 4; 3 -5 1; 1 -11 1; 1 -12 1; 2 -14 7; 1 -13 1; 
    3 -15 5; 1 -4 1; 4 -8 1; 1 -14 1; 3 -18 11; 1 -15 1; 5 -21 11];
    q = [QuadBin(q[k,:]...) for k=1:25]
    for k=1:25
        @test stabilizer(q[k]) == L[k]
    end
    
    @test WH(0,1,3,BigFloat)[2,2] ≈ (-1+sqrt(3)*im)/2
    @test WH(0,1,3)[3,3] ≈ (-1-sqrt(3)*im)/2
    @test WH(1,0,3) == [ 0  0  1; 1  0  0; 0  1  0]
    @test WH([1,0],3) == [ 0  0  1; 1  0  0; 0  1  0]
end


@testset "Zauner analytic tests" begin
    
    @test e(1) == 1
    @test e(0) == 1
    @test e(1/2) ≈ -1
    @test e(1/3) ≈ -1/2+sqrt(3)*im/2
    @test e(1/4) ≈ im
    @test e(1/8) ≈ (im+1)/sqrt(2)
    
    @test q_pochhammer(1,1,1) == 0
    @test q_pochhammer(1,1,0) == 1
    @test q_pochhammer(rand(),rand(),0) == 1
    @test q_pochhammer(1/2,rand(),1) == 1/2
    @test q_pochhammer(1/2,1,5) == (1/2)^5
    
    @test q_pochhammer_exp(rand(),rand(),0) == 1
    @test q_pochhammer_exp(sqrt(2),1,-3) ≈ (1-e(sqrt(2)))^(-3)
end


@testset "Zauner quadform tests" begin
    
    a,b,c = [1,2,3]
    Q = QuadBin(a,b,c)
    @test qmat(Q) == [1 1; 1 3]
    
    @test quadbinid(5)  == QuadBin(1,1,-1)
    @test quadbinid(8)  == QuadBin(1,2,-1)
    @test quadbinid(12) == QuadBin(1,2,-2)
    @test quadbinid(13) == QuadBin(1,3,-1)
    @test quadbinid(17) == QuadBin(1,3,-2)
    @test quadbinid(20) == QuadBin(1,4,-1)
    @test quadbinid(21) == QuadBin(1,3,-3)
    
    Q = QuadBin(1,3,-3)
    @test Q*Q == Q
    @test Q^2 == Q
    
    # example with class number 2
    Q1 = QuadBin(-3,6,2)
    Q0 = QuadBin(1,6,-6)
    @test Q0*Q0 == Q0
    @test Q0^2 == Q0
    @test Q1*Q0 == Q1
    @test Q0*Q1 == Q1
    @test Q1*Q1 == Q0
    @test Q1^2 == Q0
end


@testset "Zauner SL(2,Z) tests" begin
    S = [0 -1; 1 0]
    T = [1  1; 0 1]
    U = [1  0; 1 1]

    # Test is_sl2z function
    @test is_sl2z( [ 1  1; 1 1]) == false
    @test is_sl2z( T) == true
    @test is_sl2z( S) == true
    @test is_sl2z( S*T*S*T^4) == true
    @test is_sl2z( [-1 -3; 2 5]) == true
    
    # Tests from ghost tables up to d = 15
    @test sl2zorder( [ 3 -1; 1  0],  4) == 3
    @test sl2zorder( [ 4 -1; 1  0],  5) == 3
    @test sl2zorder( [ 5 -1; 1  0],  6) == 3
    @test sl2zorder( [ 7 -2; 4 -1],  7) == 3
    @test sl2zorder( [ 6 -1; 1  0],  7) == 3
    @test sl2zorder( [ 3 -1; 1  0],  8) == 6
    @test sl2zorder( [ 7 -1; 1  0],  8) == 3
    @test sl2zorder( [ 8 -1; 1  0],  9) == 3
    @test sl2zorder( [ 9 -5; 2 -1],  9) == 3
    @test sl2zorder( [ 9 -1; 1  0], 10) == 3
    @test sl2zorder( [11 -2; 6 -1], 11) == 3
    @test sl2zorder( [10 -1; 1  0], 11) == 3
    @test sl2zorder( [11 -4; 3 -1], 11) == 3
    @test sl2zorder( [13 -3; 9 -2], 12) == 3
    @test sl2zorder( [11 -1; 1  0], 12) == 3
    @test sl2zorder( [12 -1; 1  0], 13) == 3
    @test sl2zorder( [13 -7; 2 -1], 13) == 3
    @test sl2zorder( [13 -1; 1  0], 14) == 3
    @test sl2zorder( [14 -5; 3 -1], 14) == 3
    @test sl2zorder( [ 4 -1; 1  0], 15) == 6
    @test sl2zorder( [15 -2; 8 -1], 15) == 3
    @test sl2zorder( [14 -1; 1  0], 15) == 3
    @test sl2zorder( [16 -11;3 -2], 15) == 3
    
    # Rademacher invariant    
    @test rademacher( [1 0; 0 1]) == 0
    @test rademacher( T) == 1
    @test rademacher( S) == 0
    @test rademacher( U) == -1
    @test rademacher( [15 -47; 8 -25]) == 7
    @test rademacher( -[15 -47; 8 -25]) == 7
    @test rademacher( [-25 47; -8 15]) == -7
    @test rademacher( -S*U*U*[15 -47; 8 -25]*T*T*S) == 7
    
    
    # words in PSL(2,Z).
    # vector input
    @test psl2word([3]) == T^3
    @test psl2word([3,3]) == T^3*S*T^3
    @test psl2word([0,0]) == S
    @test psl2word([-3,0,3]) == [-1 0; 0 -1]
    # matrix input
    @test psl2word([15 -47; 8 -25]) == [2,8,-3]
    @test psl2word([15 -47; 8 -25]*T^3) == [2,8,0]
    
    r = rand(-10:10,10)
    A = psl2word(r)
    @test is_sl2z(A)
    # word decomposition should be nonnegative except on the ends
    w = psl2word(A)
    @test all(w[2:end-1] .>= 0)
    # should reproduce the original matrix up to an overall ±1
    AA = psl2word(w)
    @test AA == A || AA == -A
    
    # test stabilizer
    
    
end # SL(2,Z) testset



# Unit tests for the real double sine function
@testset "Zauner double sine tests" begin
    # Our convention for doulble sine follows Shintani (and is reciprical to K & K).
    t0 = sqrt(BigFloat(21))
    sq2 = sqrt(BigFloat(2))
    sq3 = sqrt(BigFloat(3))
    t = (5+t0)/2
    
    # Some known values from K & K
    b1 = 4
    b2 = 6
    
    @test 1/double_sine(1,b1,b2) ≈ 1
    @test 1/double_sine(2,b1,b2) ≈ sq2
    @test 1/double_sine(3,b1,b2) ≈ sq2
    @test 1/double_sine(4,b1,b2) ≈ sq3/sq2
    @test 1/double_sine(5,b1,b2) ≈ 1
    @test 1/double_sine(6,b1,b2) ≈ sq2/sq3
    @test 1/double_sine(7,b1,b2) ≈ 1/sq2
    @test 1/double_sine(8,b1,b2) ≈ 1/sq2
    @test 1/double_sine(9,b1,b2) ≈ 1
    
    # Some generic identities from K & K
    for k = 1:3
        b1 = rand(BigFloat)
        b2 = rand(BigFloat)
        
        @test 1/double_sine( (b1+b2)/2, b1, b2) ≈ 1
        @test 1/double_sine( b1/2, b1, b2) ≈ sq2
        @test 1/double_sine( b2/2, b1, b2) ≈ sq2
        @test 1/double_sine( b1 + b2/2, b1, b2) ≈ sq2/2
        @test 1/double_sine( b2 + b1/2, b1, b2) ≈ sq2/2
        @test 1/double_sine( b1, b1, b2) ≈ sqrt(b2/b1)
        @test 1/double_sine( b2, b1, b2) ≈ sqrt(b1/b2)
    end
    
    # from Shintani
    @test 1/double_sine(1//2,1,t) ≈ sq2
    @test 1/double_sine(t/2,1,t) ≈ sq2
    @test 1/(double_sine(1//3,1,t)*double_sine(1+t/3,1,t)*double_sine((2+2*t)/3,1,t)) ≈ sqrt((1+t0)/4 - sqrt((3+t0)/2)/2)
    
    # from K & W
    @test 1/double_sine(2-sq2,1,sq2) ≈ -2^(BigFloat(5)/4) * cos(pi/sq2)
    
    # need to add a few more to test the shift formulas.
end # double sine testset



@testset "Zauner utils tests" begin
    
    # Test radix function
    @test radix(5,[2,2,2]) == [1,0,1]
    @test radix(5,[2,2,2,2,2]) == [0,0,1,0,1]
    @test radix(86400,[24,60,60]) == [0,0,0]
    @test radix(86400,[7,24,60,60]) == [1,0,0,0]
    
end # utils testset

@testset "Zauner ghost tests" begin
    setprecision(100)
    q = map(x->QuadBin(x...),
      [[1 -3 1], [1 -4 1], [1 -5 1], [2 -4 1], [1 -6 1], [1 -3 1], [1 -7 1], [1 -8 1], [2 -10 5], 
      [1 -9 1], [3 -6 1], [1 -10 1], [3 -12 4], [3 -5 1], [1 -11 1], [1 -12 1], [2 -14 7], [1 -13 1], 
      [3 -15 5], [1 -4 1], [4 -8 1], [1 -14 1], [3 -18 11], [1 -15 1], [5 -21 11]])
    d = [4; 5; 6; 7; 7; 8; 8; 9; 9; 10; 11; 11; 11; 12; 12; 13; 13; 14; 14; 15; 15; 15; 15; 16; 16] 
    L = stabilizer.(q)
    n = sl2zorder.(L,d)
    A = L.^n
    β = map(x->(-BigFloat(x.b)+sqrt(BigFloat(discriminant(x))))/(2x.a),q)
    for k=1:10
        G = ghost(A[k],d[k],β[k])
        @test norm(G*G-G) <= 1e-19
    end
end