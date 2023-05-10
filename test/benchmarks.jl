using Test, ProfileSVG, BenchmarkTools, Serialize
using Zauner

#=
This code is here to benchmark the performance of the double sine integral rules. 
We run through various rules for the double sine integral formula, and check for correctness and speed. 
=#

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




@testset "Zauner ghost tests II" begin
    setprecision(100)
    q = deserialize("../data/Q")
    Q = map(x -> QuadBin(x...),q)
    d = deserialize("../data/d")
    L = stabilizer.(Q)
    n = sl2zorder.(L,d)
    A = L.^n
    β = map(x->(-BigFloat(x.b)+sqrt(BigFloat(discriminant(x))))/(2x.a),Q)
    for k=1:10
        G = ghost(A[k],d[k],β[k])
        @test norm(G*G-G) <= 1e-19
    end
end
