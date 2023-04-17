import CifToSpinFRGLattice as Cif
using FRGLatticePlotting,Plotly
using SpinFRGLattices
using FRGLatticePlotting.Plots
using SpinFRGLattices.StaticArrays
Plots.plotly()
##
function isInBox(s::Rvec,Basis::Basis_Struct,L)
    x,y,z = getCartesian(s,Basis)
    a1,a2,a3 = SpinFRGLattices.norm(Basis.a1),SpinFRGLattices.norm(Basis.a2),SpinFRGLattices.norm(Basis.a3)
    return 0 <= x <= L*a1 && 0 <= y <= L*a2 && 0 <= z <= L*a3
end
##
Basis = Cif.getBasis("../test/na6cu7bio4po44cl3_onlyCu.vesta")
# sites = generatePairSites(4,Basis)
sites = generateLUnitCells(1,Basis)
filter!(s-> isInBox(s,Basis,2),sites)
pairsPlot(sites, Basis)

Bonds = Cif.readBondsVesta("../test/na6cu7bio4po44cl3_onlyCu.vesta")

# sites1 = filter(s-> getSiteType(s,Basis) == 1,sites)
# sites2 = filter(s-> getSiteType(s,Basis) == 2,sites)
for b in Bonds
    plotDistBonds!(sites,Basis,minDist = b.minDist, maxDist = b.maxDist,lw = 10,color = Plots.Colors.RGB((b.colorRGB./255)...))
end
# plotDistBonds!(sites,Basis, minDist = Bonds[1].minDist, maxDist = Bonds[1].maxDist,lw = 10,color = :darkred)
# plotDistBonds!(sites,Basis,minDist = Bonds[2].minDist, maxDist = Bonds[2].maxDist,lw = 10,color = :black)
current()
##
"""Plot all sites and inequivalent pairs"""
function plotSystem2(System,Basis;
    plotAll = true,
    refSite = nothing,
    markersize = 5,
    inequivColor = "green",
    inequivalpha = 0.5,
    plotBonds=true,
    plotCouplings=true,
    CouplingColors = nothing,
    bondlw = 7,
    Bonds = [(minDist = Basis.NNdist-1e-3,maxDist = Basis.NNdist+1e-3,colorRGB = [0,0,0])],
    allpairs = unique!(SpinFRGLattices.sortedPairList(System.NLen,Basis)[1]),
    kwargs...)
    (;PairList,OnsitePairs )= System
    
    indices = copy(OnsitePairs)
    push!(indices,length(PairList)) # get final index
    if refSite === nothing 
        plotpairs = unique(PairList)
    else
        # allpairs = unique!(generatePairSites(System.NLen,Basis,Basis.refSites[refSite]))
        plotpairs = PairList[indices[refSite]:indices[refSite+1]]
    end
    filter!(x-> x in allpairs,plotpairs)

    plotAll || (allpairs = plotpairs)
    pl = pairsPlot(allpairs,Basis,markersize = markersize,aspect_ratio=:equal;kwargs...)
    if plotBonds
        for b in Bonds
            plotDistBonds!(allpairs,Basis,minDist = b.minDist, maxDist = b.maxDist,lw = bondlw,color = Plots.Colors.RGB((b.colorRGB./255)...))
        end
    end
    # plotBonds && plotDistBonds!(allpairs,Basis;color = Bondcolor,lw = bondlw, minDist = bondDist-1e-3, maxDist = bondDist+1e-3)

    plotAll && pairsPlot(plotpairs,Basis,pl,color = inequivColor,alpha = inequivalpha,markersize = 2*markersize)
    
    plotCouplings && plotCouplings!(System,Basis;refSite = refSite,colors = CouplingColors)
    return pl
end
##
Basis = Cif.getBasis("../test/na6cu7bio4po44cl3_onlyCu.vesta")
Bonds = Cif.readBondsVesta("../test/na6cu7bio4po44cl3_onlyCu.vesta")


function generateLayer(L,Basis,refsite)
    allpairs = generateLUnitCells(L,Basis,refsite)
    filter!(R->R.n3 == 0,allpairs)
    return allpairs
end
##
addSyms = let 
    refsitesFrac = [Basis.T*getCartesian(r,Basis) for r in Basis.refSites]
    traf(org,rot) = Cif.fractionaltransformation_origin(org,rot)
    xmirror = float(SA[-1 0 0;0 1 0;0 0 1])
    ymirror = float(SA[1 0 0;0 -1 0;0 0 1])
    xyrotation = float(SA[0 -1 0;1 0 0;0 0 1])
    inversion = float(SA[-1 0 0;0 -1 0;0 0 -1])

    Cif.FractionalTransformation(Basis.T * Basis.b[Basis.refSites[2].b],SA[-1 0. 0;0 -1 0;0 0 -1])

    Site2inversion = traf(refsitesFrac[2],inversion)

    println(refsitesFrac[2])
    println(traf(refsitesFrac[2],xmirror)(refsitesFrac[2]))

    [Site2inversion,traf(refsitesFrac[1],xmirror), traf(refsitesFrac[1],ymirror),traf(refsitesFrac[1],xyrotation)]
    # [traf(refsitesFrac[1],xmirror), traf(refsitesFrac[1],ymirror)]
    
end
CPyro = Cif.generateSystem(2,"../test/na6cu7bio4po44cl3_onlyCu.vesta",method = generateLayer;addSyms)
println(CPyro.PairList|> length)
##
allpairs = generateLayer(2,Basis,Basis.refSites[1])
plotSystem2(CPyro,Basis;refSite = 2,allpairs,bondDist = Basis.NNdist,Bonds,markersize = 3,plotAll = true,bondlw = 1)
zlims!(-20,30)
# pairsPlot(CPyro.PairList,Basis)
##