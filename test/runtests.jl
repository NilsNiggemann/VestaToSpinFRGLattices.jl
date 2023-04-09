using CifToSpinFRGLattice,Test,LinearAlgebra,StaticArrays

@testset "parseInequivSites" begin
    @test getInequivalentSites("SimpleCubic.cif") == [[0.0,0.0,0.0]]
end

@testset "parseSymmetries" begin
    symlist = [
        "'   '-z, -x, y''"
        "'   '-z+1/4, -x-0.1, y*2''"
    ]
    syms = CifToSpinFRGLattice.modifySymmetries!(symlist)
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
@testset "convert lattice angles to unit vectors" begin
    a = SA[1,0,0]
    b_(x) = LinearAlgebra.normalize(SA[1,x,0])
    c_(x,y) = LinearAlgebra.normalize(SA[1,x,y])

    range = LinRange(-2,2,30)
    for x in range
        b = b_(x)
        for y in range
            c = c_(x,y)
            alpha,beta,gamma = CifToSpinFRGLattice.latticeparameters(a,b,c)

            a2,b2,c2 =  CifToSpinFRGLattice.getUnitVectors(alpha,beta,gamma)
            alpha2,beta2,gamma2 = CifToSpinFRGLattice.latticeparameters(a2,b2,c2)

            @test alpha ≈ alpha2
            @test beta ≈ beta2
            @test gamma ≈ gamma2
        end
    end
end