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

function getUnitCellTranslation(R::Rvec,Basis)
    (;n1,n2,n3,b) = R

    T = FractionalTransformation(SA[1 0 0 -n1;
           0 1 0 -n2;
           0 0 1 -n3;
           0 0 0 1])

    return siteTransformation(T,Basis)
end

function generateSymms(irreps::AbstractVector{<:SiteTransformation};digits = 14,maxiter = 1000)
    syms = Set(irreps)
    Basis = first(irreps).Basis
    lower = minimum(minimum.(Basis.b))-1e-4
    upper = maximum(maximum.(Basis.b))+1e-4
    iteration = 0
    for s1 in syms
        for s2 in syms
            if iteration > maxiter
                @warn "Too many iterations"
                return collect(syms)
            end
            Snew = s1*s2
            
            
            isInUnitCell(getTranslation(Snew),lower,upper) || continue

            push!(syms,Snew)
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
