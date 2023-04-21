
"""writes crucial Information about reduced Lattice to an HDF5 file
The data layout includes a list of all sites appearing in pairNumbersDict and a vector of all pairNumbers. If a pair does not appear in the dict, write 0 instead.
"""
function h5writeReducedLattice(filename::AbstractString,key,pairNumbersDict)
    R1,R2,pairNumber = SpinFRGLattices.convertDictToArrays(pairNumbersDict)
    AllSites = unique(R1)
    getPNumber(pair) = pair in keys(pairNumbersDict) ? pairNumbersDict[pair] : 0

    allPairNumbers = [getPNumber((R1,R2)) for R1 in AllSites, R2 in AllSites]
    
    h5writeRvecs(filename,key*"/AllSites",AllSites)
    h5write(filename,key*"/pairNumbers",allPairNumbers)
end

function h5writeRvecs(filename::AbstractString,key::AbstractString,RV::AbstractArray{R}) where R <:Rvec
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

function h5readPairNumberDict(filename::AbstractString,key::AbstractString)
    fields = h5read(filename,key)
    AllSites = h5readRvecs(filename,key*"/AllSites")
    pairNumber = fields["pairNumbers"]
    
    iterations = Iterators.product(AllSites,AllSites)

    D = Dict(R1R2 => pairNumber[i] for (i,R1R2) in enumerate(iterations))
    filter!(x->x.second != 0,D)
    return D
end