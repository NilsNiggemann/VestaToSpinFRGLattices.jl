using CifToSpinFRGLattice

using Test,LinearAlgebra,StaticArrays

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
    @test length(syms) == 48
    
    x,y,z = r = [10.5,11.2,-10.4]
    r´ = [-y, -z, x]

    @test syms[1](r) == r # identity
    @test syms[23](r) == r´
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
    @test getproperty.(sites,:position) == [
        SA[0.,0.,0.],
        SA[0.25,0.25,0.25],
    ]
end

##
@testset "reading Vesta Symmetries" begin
    syms = getSymmetries("SimpleCubic.cif")
    syms2 = Cif.readVestaSymops("SimpleCubic.vesta")
    @test length(syms) == length(syms2)

    x,y,z = r = SA[10.5,11.2,-10.4]
    for (s1,s2) in zip(syms,syms2)
        @test s1(r) == s2(r)

    end

end
##

@testset "Centered Pyrochlore Unit Cell" begin
    
    UCSites = sort!([
        SA[0.000000, 0.000000, 0.000000],
        SA[0.250000, 0.000000, 0.250000],
        SA[0.250000, 0.250000, 0.000000],
        SA[0.000000, 0.250000, 0.250000],
        SA[0.000000, 0.500000, 0.500000],
        SA[0.250000, 0.500000, 0.750000],
        SA[0.250000, 0.750000, 0.500000],
        SA[0.000000, 0.750000, 0.750000],
        SA[0.500000, 0.000000, 0.500000],
        SA[0.750000, 0.000000, 0.750000],
        SA[0.750000, 0.250000, 0.500000],
        SA[0.500000, 0.250000, 0.750000],
        SA[0.500000, 0.500000, 0.000000],
        SA[0.750000, 0.500000, 0.250000],
        SA[0.750000, 0.750000, 0.000000],
        SA[0.500000, 0.750000, 0.250000],
        SA[0.125000, 0.125000, 0.125000],
        SA[0.375000, 0.375000, 0.875000],
        SA[0.875000, 0.875000, 0.875000],
        SA[0.125000, 0.625000, 0.625000],
        SA[0.375000, 0.875000, 0.375000],
        SA[0.625000, 0.125000, 0.625000],
        SA[0.875000, 0.375000, 0.375000],
        SA[0.625000, 0.625000, 0.125000]
    ])
    Sites = Cif.getUnitCell("CentredPyrochlore.vesta").positions |> sort!
    @test length(UCSites) == length(Sites)

    smallerArray, largerArray = length(UCSites) > length(Sites) ? (UCSites,Sites) : (Sites,UCSites)

    uniqueSites = [s for (i,s) in enumerate(smallerArray) if s ∉ largerArray]

    @test isempty(uniqueSites)
end
##
using CifToSpinFRGLattice.SpinFRGLattices
@testset "transforms" begin
    
    syms = Cif.readVestaSymops("SimpleCubic.vesta")
    Basis = Cif.getBasis("SimpleCubic.vesta")

    T1 = Cif.gettransform(syms[1],Basis)
    T2 = Cif.gettransform(syms[3],Basis)

    R1 = Rvec(1,2,3,1)
    R2 = Rvec(-2,2,3,1)
    @test T1(R1) == R1
    @test T1(R2) == R2
    @test T2(R1) == Rvec(-1,-2,3,1)
    @test T2(R2) == Rvec(2,-2,3,1)
end
##
a = Cif.generateSystem(7,"SimpleCubic.vesta",test = true)
##
@testset "Read Cubic" begin
    b = SimpleCubic.getCubic(7)
    @test length(a.PairList) == length(b.PairList)
end
a= Cif.generateSystem(7,"CentredPyrochlore.vesta",test = true)
##
@testset "PairDict IO" begin
    file = tempname()
    a = Cif.generateReducedLattice(6,"CentredPyrochlore.vesta").pairNumberDict
    Cif.h5writeReducedLattice(file,"/test",a)
    b = Cif.h5readPairNumberDict(file,"/test")
    @test a == b
end