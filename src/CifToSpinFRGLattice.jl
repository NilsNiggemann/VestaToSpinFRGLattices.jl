module CifToSpinFRGLattice
    
    using StaticArrays

    """reading cif file returns information under key as a vector of strings"""
    function readCifInfo(filename::String,key::String)
        positions = String[]
        open(filename) do file
            for line in eachline(file)
                if occursin(key,line)
                    for line in eachline(file)
                        if occursin("_",line) || line == ""
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

    function getInequivalentSites(poslist::Vector{String})
        spl = split.(poslist," ")
        filter!.(!=(""),spl)
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

    export getInequivalentSites,getSymmetries
end # module SpinLatticeFromCif
