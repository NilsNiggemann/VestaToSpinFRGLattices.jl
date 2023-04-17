function generateSymms(irreps;digits = 14)
    syms = Set(irreps)
    roundsyms = Set(irreps) # rounded version to avoid duplicates

    iteration = 0
    for s1 in syms
        for s2 in syms
            iteration += 1
            if iteration > 1000
                @warn "Too many iterations"
                return collect(syms)
            end
            # iteration >1000 && error("Too many iterations")
            Snew = s1*s2

            isInUnitCell(getTranslation(Snew)) || continue

            roundSnew = round(Snew; digits)
            if roundSnew ∉ roundsyms
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

function testIrrs()
    transl = SA[0,0,0]
    xmirr = FractionalTransformation(transl,SA[-1 0 0; 0 1 0; 0 0 1] )
    ymirr = FractionalTransformation(transl,SA[1 0 0; 0 -1 0; 0 0 1] )
    zmirr = FractionalTransformation(transl,SA[1 0 0; 0 1 0; 0 0 -1] )
    xyRot = FractionalTransformation(transl,SA[0 -1 0; 1 0 0; 0 0 1] )
    xzRot = FractionalTransformation(transl,SA[0 0 -1; 0 1 0; 1 0 0] )
    yzRot = FractionalTransformation(transl,SA[1 0 0; 0 0 -1; 0 1 0] )
    return Set([xmirr,ymirr,zmirr,xyRot,xzRot,yzRot])
end