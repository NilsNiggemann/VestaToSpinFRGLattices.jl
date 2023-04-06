module CifToSpinFRGLattice
    
    using StaticArrays
    """reading cif file returns the symmetries specified under _space_group_symop_operation_xyz as a string"""
    function readCifSymmetries(filename::String)
        file = open(filename)
        symmetries = String[]
        for line in eachline(file)
            if occursin("_space_group_symop_operation_xyz",line)
                for line in eachline(file)
                    if occursin("_",line) || line == ""
                        break
                    end
                    push!(symmetries,line)
                end
                break
            end
        end
        close(file)
        return symmetries
    end
    
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

    function getSymmetries(symlist::Vector{String})
        [
            fncFromString(s) for s in symlist
        ]
    end

    getSymmetries(filename::String) = filename |> readCifSymmetries |> modifySymmetries! |> getSymmetries

    export readCifSymmetries,modifySymmetries!,modifySymmetries, fncFromString,getSymmetries
end # module SpinLatticeFromCif
