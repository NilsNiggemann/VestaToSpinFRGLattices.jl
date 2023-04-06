using CifToSpinFRGLattice,Test

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