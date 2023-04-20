import CifToSpinFRGLattice as Cif
using FRGLatticePlotting,Plotly
using SpinFRGLattices
import SpinFRGLattices as SL
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
    xymirror = float(SA[0 1 0; 1 0 0;0 0 1])

    inversion = float(SA[-1 0 0;0 -1 0;0 0 -1])

    site1trafos = (xmirror,ymirror,xyrotation,xymirror,inversion)
    site1trafos = []

    site2trafos = (inversion)
    site2trafos = []
    # site3trafos = (ymirror,)
    site3trafos = []

    trafs = (site1trafos,site2trafos,site3trafos)
    
    syms = [traf(refsitesFrac[x],y) for x in eachindex(refsitesFrac) for y in trafs[x] ]
end
##
Squag = Cif.generateSystem(1,"../test/na6cu7bio4po44cl3_onlyCu.vesta",method = generateLayer;addSyms)
println(Squag.PairList|> length)
##
Basis = Cif.getBasis("../test/na6cu7bio4po44cl3_onlyCu.vesta")
Bonds = Cif.readBondsVesta("../test/na6cu7bio4po44cl3_onlyCu.vesta")

allpairs = generateLayer(1,Basis,Basis.refSites[1])
plotSystem(Squag,Basis;refSite = 2,allpairs,bondDist = Basis.NNdist,Bonds,markersize = 3,plotAll = true,bondlw = 1)
zlims!(-20,30)
# pairsPlot(CPyro.PairList,Basis)
##
Basis = Cif.getBasis("../test/CentredPyrochlore.vesta")
CPyro = Cif.generateSystem(5,"../test/CentredPyrochlore.vesta")
##
Bonds = Cif.readBondsVesta("../test/CentredPyrochlore.vesta")
plotSystem(CPyro,Basis;refSite = 1,bondDist = Basis.NNdist,Bonds,markersize = 3,plotAll = true,bondlw = (3,1))