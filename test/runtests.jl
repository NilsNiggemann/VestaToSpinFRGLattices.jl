using VestaToSpinFRGLattices

using Test,LinearAlgebra,StaticArrays

import VestaToSpinFRGLattices as Vesta
##
@testset "convert lattice angles to unit vectors" begin
    a = SA[1,0,0]
    b_(x) = LinearAlgebra.normalize(SA[1,x,0])
    c_(x,y) = LinearAlgebra.normalize(SA[1,x,y])

    range = LinRange(-2,2,10)
    for x in range
        b = b_(x)
        for y in range
            c = c_(x,y)
            alpha,beta,gamma = Vesta.latticeparameters(a,b,c)

            a2,b2,c2 =  Vesta.getUnitVectors(alpha,beta,gamma)
            alpha2,beta2,gamma2 = Vesta.latticeparameters(a2,b2,c2)

            @test alpha ≈ alpha2
            @test beta ≈ beta2
            @test gamma ≈ gamma2
        end
    end
end
##
@testset "Vesta readLatticevecs"  begin
    
    a,b,c = Vesta.readLatticeVectors("SimpleCubic.vesta")
    @test a ≈ SA[1,0,0]
    @test b ≈ SA[0,1,0]
    @test c ≈ SA[0,0,1]
end
@testset "Vesta generateSymms"  begin
    
    transl = float(SA[0,0,0])
    xmirr = Vesta.FractionalTransformation(transl,float(SA[-1 0 0; 0 1 0; 0 0 1]) )
    ymirr = Vesta.FractionalTransformation(transl,float(SA[1 0 0; 0 -1 0; 0 0 1]) )
    zmirr = Vesta.FractionalTransformation(transl,float(SA[1 0 0; 0 1 0; 0 0 -1]) )
    xyRot = Vesta.FractionalTransformation(transl,float(SA[0 -1 0; 1 0 0; 0 0 1]) )
    xzRot = Vesta.FractionalTransformation(transl,float(SA[0 0 -1; 0 1 0; 1 0 0]) )
    yzRot = Vesta.FractionalTransformation(transl,float(SA[1 0 0; 0 0 -1; 0 1 0]) )
    
    testIrrs = [xmirr,ymirr,zmirr,xyRot,xzRot,yzRot]
    generatedSymms = Vesta.generateSymms(testIrrs)

    readSymms = Vesta.readSymops("SimpleCubic.vesta")

    @test Set(generatedSymms) == Set(readSymms)
end
##

##
@testset "vesta readSites" begin
    sites = Vesta.readUniqueSites("test.vesta")
    @test getproperty.(sites,:position) == [
        SA[0.,0.,0.],
        SA[0.25,0.25,0.25],
    ]
end

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
    Sites = Vesta.getUnitCell("CentredPyrochlore.vesta").positions |> sort!
    @test length(UCSites) == length(Sites)

    smallerArray, largerArray = length(UCSites) > length(Sites) ? (UCSites,Sites) : (Sites,UCSites)

    uniqueSites = [s for (i,s) in enumerate(smallerArray) if s ∉ largerArray]

    @test isempty(uniqueSites)
end
##
using VestaToSpinFRGLattices.SpinFRGLattices
@testset "transforms" begin
    
    syms = Vesta.getSymmetries("SimpleCubic.vesta")
    Basis = Vesta.getBasis("SimpleCubic.vesta")

    T1 = syms.refSyms[1][1]
    T2 = syms.refSyms[1][3]

    R1 = Rvec(1,2,3,1)
    R2 = Rvec(-2,2,3,1)
    @test T1(R1) == R1
    @test T1(R2) == R2
    @test T2(R1) == Rvec(-1,-2,3,1)
    @test T2(R2) == Rvec(2,-2,3,1)
end
##
@testset "Read Cubic" begin
    a = Vesta.generateSystem(7,"SimpleCubic.vesta",test = true)
    b = SimpleCubic.getCubic(7)
    @test length(a.PairList) == length(b.PairList)
end

Vesta.generateSystem(7,"CentredPyrochlore.vesta",test = true)
##
@testset "reduced Geometry IO" begin
    file = tempname()
    a = Vesta.generateReducedLattice(6,"CentredPyrochlore.vesta")
    Vesta.h5saveReducedLattice(file,a)
    b = Vesta.h5readPairNumberDict(file)

    @testset "pairNumberDict" begin
        @test a.pairNumberDict == b
    end
    @testset "Basis" begin
        Basis = Vesta.h5readBasis(file)
        for f in fieldnames(Basis_Struct_3D)
            @test getfield(Basis,f) == getfield(a.Basis,f)
        end

    end
end