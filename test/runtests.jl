using Test,LinearAlgebra,StaticArrays
using CifToSpinFRGLattice
import CifToSpinFRGLattice as Cif
##
@testset "parseInequivSites" begin
    @test getInequivalentSites("SimpleCubic.cif") == [[0.0,0.0,0.0]]
end

@testset "parseSymmetries" begin
    symlist = [
        "'   '-z, -x, y''"
        "'   '-z+1/4, -x-0.1, y*2''"
    ]
    syms = Cif.modifySymmetries!(symlist)
    @test syms == [
        "r-> SA[-r[3],-r[1],r[2]]"
        "r-> SA[-r[3]+1/4,-r[1]-0.1,r[2]*2]"
    ]
    syms = getSymmetries(syms)
    x,y,z = r = [1,2,3]
    @test syms[1](r) == [-3,-1,2]
    @test syms[2](r) == [-z+1/4, -x-0.1, y*2]
end

@testset "reading Symmetries" begin
    syms = getSymmetries("SimpleCubic.cif")
    @test length(syms) == 12
    
    x,y,z = r = [10.5,11.2,-10.4]
    r´ = [-y, -z, x]

    @test syms[1](r) == r # identity
    @test syms[end](r) == r´
end
##
@testset "convert lattice angles to unit vectors" failfast = true begin
    a = SA[1,0,0]
    b_(x) = LinearAlgebra.normalize(SA[1,x,0])
    c_(x,y) = LinearAlgebra.normalize(SA[1,x,y])

    range = LinRange(-2,2,30)
    for x in range
        b = b_(x)
        for y in range
            c = c_(x,y)
            alpha,beta,gamma = Cif.latticeparameters(a,b,c)

            a2,b2,c2 =  Cif.getUnitVectors(alpha,beta,gamma)
            alpha2,beta2,gamma2 = Cif.latticeparameters(a2,b2,c2)

            @test alpha ≈ alpha2
            @test beta ≈ beta2
            @test gamma ≈ gamma2
        end
    end
end
##
@testset "Vesta readLatticevecs"  begin
    
    a,b,c = Cif.readLatticeVectors("SimpleCubic.vesta")
    @test a ≈ SA[1,0,0]
    @test b ≈ SA[0,1,0]
    @test c ≈ SA[0,0,1]
end
##
@testset "vesta readSites" begin
    sites = Cif.readVestaSites("test.vesta")
    @test sites == [
        SA[0.,0.,0.],
        SA[0.25,0.25,0.25],
    ]
end