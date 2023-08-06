
"""writes crucial Information about reduced Lattice to an HDF5 file
The data layout includes a list of all sites appearing in pairNumbersDict and a vector of all pairNumbers. If a pair does not appear in the dict, write 0 instead.
"""
function h5saveReducedLattice(filename::AbstractString,reducedLattice,key = standardKey())
    pairNumberDict = reducedLattice.pairNumberDict
    R1,R2,pairNumber = SpinFRGLattices.convertDictToArrays(pairNumberDict)
    AllSites = unique(R2)

    UC = SpinFRGLattices.getUnitCell(reducedLattice.Basis)
    sort!(UC,by = x -> x.b)

    allPairNumbers = [pairNumberDict[(R1,R2)] for R1 in UC, R2 in AllSites]
    
    h5saveRvecs(filename,key*"/AllSites",AllSites)
    h5open(filename,"cw") do file
        file[key*"/pairNumbers",blosc =9] = allPairNumbers
    end
    h5saveBasis(filename,reducedLattice.Basis)

end


function h5saveBasis(filename,Basis::Bas,key = standardKeyBasis()) where Bas <: Basis_Struct
    fields = fieldnames(Bas)
    
    conv(x) = x
    conv(x::StaticArray) = Array(x)
    conv(x::AbstractArray{<:AbstractArray}) = Array(hcat(x...))
    BasisDict = Dict(string(k) => conv(getfield(Basis,k)) for k in fields)
    
    delete!(BasisDict,"bLatt")
    delete!(BasisDict,"T")

    h5write(filename,key,BasisDict)
end

function h5readBasis3D(filename,key = standardKeyBasis())
    fields = h5read(filename,key)
    b = SVector{3}.(collect(eachcol(fields["b"])))
    refSites = [Rvec(x.n1,x.n2,x.n3,x.b) for x in fields["refSites"]]
    
    fieldsTup = (Symbol(k) => v for (k,v) in fields)
    return Basis_Struct_3D(;fieldsTup...,b,refSites)
end

function h5readBasis2D(filename,key = standardKeyBasis())
    fields = h5read(filename,key)
    b = SVector{2}.(collect(eachcol(fields["b"])))
    refSites = [Rvec(x.n1,x.n2,x.b) for x in fields["refSites"]]
    
    fieldsTup = (Symbol(k) => v for (k,v) in fields)
    return Basis_Struct_3D(;fieldsTup...,b,refSites)
end

function h5readBasis(filename,key = standardKeyBasis())
    Dim = h5read(filename,key*"/b")
    if size(Dim,1) == 3
        return h5readBasis3D(filename,key)
    else
        return h5readBasis2D(filename,key)
    end
end

standardKey() = "ReducedLattice"
standardKeyBasis() = "ReducedLattice/Basis"

function h5saveRvecs(filename::AbstractString,key::AbstractString,RV::AbstractArray{R}) where R <:Rvec
    h5open(filename,"w") do file
        for k in fieldnames(R)
            file[key*"/"*string(k)] = getfield.(RV,k)
        end
    end
end

function h5readRvecs(filename::AbstractString,key::AbstractString)
    fields = h5read(filename,key)
    return convertDictToRvecs(fields)
end

function convertDictToRvecs(D)
    if "n3" in keys(D)
        return Rvec.(D["n1"],D["n2"],D["n3"],D["b"])
    else
        return Rvec.(D["n1"],D["n2"],D["b"])
    end

end

function h5readPairNumberDict(filename::AbstractString,key::AbstractString = standardKey())
    fields = h5read(filename,key)
    AllSites = h5readRvecs(filename,key*"/AllSites")
    # UC = SpinFRGLattices.getUnitCell(Basis)
    UC = filter(isInUnitCell,AllSites)
    sort!(UC,by = x -> x.b)

    pairNumber = fields["pairNumbers"]
    
    iterations = Iterators.product(UC,AllSites)

    D = Dict(R1R2 => Int(pairNumber[i]) for (i,R1R2) in enumerate(iterations))
    filter!(x->x.second != 0,D)
    return PairNumbersDict(D)
end

function h5readReducedLattice(filename::AbstractString,key::AbstractString = standardKey())
    pairNumberDict = h5readPairNumberDict(filename,key)
    Basis = h5readBasis(filename,key*"/Basis")
    return (;pairNumberDict,Basis)
end