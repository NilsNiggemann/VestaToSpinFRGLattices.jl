module VestaToSpinFRGLattices
    
    using StaticArrays,HDF5

    abstract type CrystallographicFile end
    struct VestaFile <:CrystallographicFile end

    """given a line from a file returns true if it contains any word. If it contains an empty space or a selection of numbers return false"""
    containsData(line::String,::CrystallographicFile) = !(occursin("_",line) || line == "")

    """reading cif file returns information under key as a vector of strings"""
    function readFileInfo(filename::String,key::String,FileTypye::CrystallographicFile)
        positions = String[]
        open(filename) do file
            for line in eachline(file)
                if occursin(key,line)
                    for line in eachline(file)
                        # if occursin("_",line) || line == ""
                        if !containsData(line,FileTypye)
                            break
                        end
                        push!(positions,line)
                    end
                    break
                end
            end
        end
        return positions
    end

    """returns whether is symmetry leaves reference site invariant"""
    isRefSymmetry(sym,refSite) = sym(refSite) == refSite

    """returns whether a symmetry leaves any reference site invariant"""
    function isAnyRefSymmetry(sym,refSites::AbstractArray)
        for refSite in refSites
            if isRefSymmetry(sym,refSite)
                return true
            end
        end
        return false
    end

    getSymmetriesForRefSite(refSite,syms) = [s for s in syms if isRefSymmetry(refSite,s)]

    function getFields(line::String)
        fields = split(line," ")
        filter!(!=(""),fields)
        return fields
    end
    
    parseAsSVector(V::AbstractVector{<:AbstractString}) = SVector{3,Float64}([parse(Float64,s) for s in V])

    include("readVesta.jl")
    include("generateSymmetriesFromIrrep.jl")
    include("IO.jl")

    export getBasis,generateReducedLattice,generateSymmetries,readBonds,generateSystem
end # module SpinLatticeFromCif
