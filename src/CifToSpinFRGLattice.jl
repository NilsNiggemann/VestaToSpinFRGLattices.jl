module CifToSpinFRGLattice
    
    using StaticArrays

    abstract type CrystallographicFile end
    struct CifFile <:CrystallographicFile end
    struct VestaFile <:CrystallographicFile end

    """given a line from a file returns true if it contains any word. If it contains an empty space or a selection of numbers return false"""
    containsData(line::String,::CifFile) = !(occursin("_",line) || line == "")

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

    getSymmetriesForRefSite(refSite,syms) = [s for s in syms if s(refSite) == refSite]

    """reading cif file returns positions of symmetry inequivalent sites"""
    readCifPositions(filename::String) = readFileInfo(filename,"_atom_site_type_symbol",CifFile())

    function getFields(line::String)
        fields = split(line," ")
        filter!(!=(""),fields)
        return fields
    end
    
    parseAsSVector(V::AbstractVector{<:AbstractString}) = SVector{3,Float64}([parse(Float64,s) for s in V])

    function getInequivalentSites(poslist::Vector{String})
        spl = getFields.(poslist)
        pos = [parseAsSVector(p[3:5]) for p in spl]
    end
    getInequivalentSites(filename::String) = filename |> readCifPositions |> getInequivalentSites

    """reading cif file returns the symmetries specified under _space_group_symop_operation_xyz as a string"""
    readCifSymmetries(filename::String) = readFileInfo(filename,"_space_group_symop_operation_xyz",CifFile())
    
    """given a vector of symmetric positions as a strings modify the strings such that they can be evaluated by julia"""
    function modifySymmetries!(symmetries::Vector{String})
        for (i,s) in enumerate(symmetries)
            s = replace(s,
            " " => "",
            "'" => "",
            "x" => "r[1]",
            "y" => "r[2]",
            "z" => "r[3]",
            )
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
