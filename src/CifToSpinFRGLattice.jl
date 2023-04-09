module CifToSpinFRGLattice
    
    using StaticArrays

    """given a line from a file returns true if it contains any word. If it contains an empty space or a selection of numbers return false"""
    function containsWord(line::String)
        if occursin(" ",line) || occursin("[0-9]",line)
            return false
        end
        return true
    end

    """reading cif file returns information under key as a vector of strings"""
    function readCifInfo(filename::String,key::String)
        positions = String[]
        open(filename) do file
            for line in eachline(file)
                if occursin(key,line)
                    for line in eachline(file)
                        # if occursin("_",line) || line == ""
                        if containsWord(line)
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

    getSymmetriesForRefSite(refSite,syms) = [s for s in syms if s(refSite) == refSite]

    """reading cif file returns positions of symmetry inequivalent sites"""
    readCifPositions(filename::String) = readCifInfo(filename,"_atom_site_type_symbol")

    function getFields(line::String)
        fields = split(line," ")
        filter!(!=(""),fields)
        return fields
    end
    
    function getInequivalentSites(poslist::Vector{String})
        spl = getFields.(poslist)
        return [SVector{3,Float64}([parse(Float64,s) for s in p[3:5]]) for p in spl]
    end
    getInequivalentSites(filename::String) = filename |> readCifPositions |> getInequivalentSites

    """reading cif file returns the symmetries specified under _space_group_symop_operation_xyz as a string"""
    readCifSymmetries(filename::String) = readCifInfo(filename,"_space_group_symop_operation_xyz")
    
    """given a vector of symmetric positions as a strings modify the strings such that they can be evaluated by julia"""
    function modifySymmetries!(symmetries::Vector{String})
        for (i,s) in enumerate(symmetries)
            s = replace(s," " => "")
            s = replace(s,"'" => "")
            s = replace(s,"x" => "r[1]")
            s = replace(s,"y" => "r[2]")
            s = replace(s,"z" => "r[3]")
            symmetries[i] = join(("r-> SA[",s,"]"))
            # symmetries[i] = join(s)
        end
        return symmetries
    end
    modifySymmetries(symmetries::Vector{String}) = modifySymmetries!(copy(symmetries))

    fncFromString(s) = eval(Meta.parse(s))

    getSymmetries(symlist::Vector{String}) = fncFromString.(symlist)

    getSymmetries(filename::String) = filename |> readCifSymmetries |> modifySymmetries! |> getSymmetries

    include("readVesta.jl")
    export getInequivalentSites,getSymmetries
end # module SpinLatticeFromCif
