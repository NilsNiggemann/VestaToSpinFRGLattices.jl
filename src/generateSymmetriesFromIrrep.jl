function generateSymms(irreps)
    syms = Set(irreps)
    for iteration in 1:100
        for s1 in syms
            for s2 in syms
                iteration >1000 && error("Too many iterations")
                if s1 != s2 #probably need to get rid of translations
                    push!(syms,s1*s2)
                end
            end
        end
    end
    return collect(syms)
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