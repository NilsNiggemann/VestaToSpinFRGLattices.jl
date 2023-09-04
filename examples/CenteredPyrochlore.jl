import VestaToSpinFRGLattices as Vesta
using VestaToSpinFRGLattices
using FRGLatticePlotting,Plotly
using SpinFRGLattices
import SpinFRGLattices as SL
using SpinFRGLattices.StaticArrays
using FRGLatticePlotting.Makie.Colors
using FRGLatticePlotting.Makie
##
function isInBox(s::Rvec,Basis::Basis_Struct,L)
    x,y,z = getCartesian(s,Basis)
    a1,a2,a3 = SpinFRGLattices.norm(Basis.a1),SpinFRGLattices.norm(Basis.a2),SpinFRGLattices.norm(Basis.a3)
    return 0 <= x <= L*a1 && 0 <= y <= L*a2 && 0 <= z <= L*a3
end
##
Basis = Vesta.getBasis("../test/na6cu7bio4po44cl3_onlyCu.vesta")
# sites = generatePairSites(4,Basis)
sites = generateLUnitCells(1,Basis)
filter!(s-> isInBox(s,Basis,2),sites)
pairsPlot(sites, Basis)

Bonds = Vesta.readBonds("../test/na6cu7bio4po44cl3_onlyCu.vesta")

for b in Bonds
    plotDistBonds!(sites,Basis,minDist = b.minDist, maxDist = b.maxDist,lw = 10,color = Colors.RGB((b.colorRGB./255)...))
end
current_figure()
##


function generateLayer(L,Basis,refsite)
    allpairs = generateLUnitCells(L,Basis,refsite)
    filter!(R->R.n3 == 0,allpairs)
    return allpairs
end
##
addSyms = let 
    refsitesFrac = [Basis.T*getCartesian(r,Basis) for r in Basis.refSites]
    traf(org,rot) = Vesta.fractionaltransformation_origin(org,rot)

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
Squag = Vesta.generateSystem(1,"../test/na6cu7bio4po44cl3_onlyCu.vesta",method = generateLayer;addSyms)
println(Squag.PairList|> length)
##
Basis = Vesta.getBasis("../test/na6cu7bio4po44cl3_onlyCu.vesta")
Bonds = Vesta.readBonds("../test/na6cu7bio4po44cl3_onlyCu.vesta")

allpairs = generateLayer(1,Basis,Basis.refSites[1])
fig,ax = FRGLatticePlotting.getStandardFigure(Rvec_3D;
xlabelvisible = false,xticklabelsvisible = false,xticksvisible = false,xspinesvisible = false,xgridvisible = false,
ylabelvisible = false,yticklabelsvisible = false,yticksvisible = false,yspinesvisible = false,ygridvisible = false,
zlabelvisible = false,zticklabelsvisible = false,zticksvisible = false,zspinesvisible = false,zgridvisible = false,
)
FRGLatticePlotting.plotSystem!(ax,Squag,Basis;refSite = 1,allpairs,bondDist = Basis.NNdist,Bonds,markersize = 20,plotAll = true,bondlw = 2,inequivScale = 2.0)
zlims!(ax,-10,20)
fig
# pairsPlot(CPyro.PairList,Basis)
##
Basis = getBasis("../test/CentredPyrochlore.vesta")
CPyro = generateSystem(14,"../test/CentredPyrochlore.vesta")
setNeighborCouplings!(CPyro,[5,2.],Basis)
##
Bonds = Vesta.readBonds("../test/CentredPyrochlore.vesta")
plotSystem(CPyro,Basis;refSite = 1,bondDist = Basis.NNdist,markersize = 3,Bonds,plotAll = true,plotCouplings = true,bondlw = (3,1))
