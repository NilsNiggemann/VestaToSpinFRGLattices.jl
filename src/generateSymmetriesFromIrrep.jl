function generateSymms(irreps::AbstractVector{<:FractionalTransformation};digits = 14,maxiter = 1000)
    syms = Set(irreps)
    roundsyms = Set(copy(irreps)) # rounded version to avoid duplicates

    iteration = 0
    for s1 in syms
        for s2 in syms
            if iteration > maxiter
                @warn "Too many iterations"
                return collect(syms)
            end
            Snew = s1*s2
            
            isInUnitCell(getTranslation(Snew)) || continue
            roundSnew = round(Snew; digits)
            if roundSnew ∉ roundsyms
                iteration += 1
                push!(syms,Snew)
                push!(roundsyms,roundSnew)
            end
        end
    end
    return collect(syms)
end

import Base.round
function round(s::FractionalTransformation;kwargs...)
    T = round.(s.WMatrix;kwargs...)
    FractionalTransformation(T)
end

function translate_back_to_UC(s::FractionalTransformation) 
    org = getTranslation(s)
    rot = getRotation(s)
    org = org - floor.(org)
    return FractionalTransformation(org,rot)
end
