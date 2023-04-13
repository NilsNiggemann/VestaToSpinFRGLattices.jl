import CifToSpinFRGLattice as Cif
using FRGLatticePlotting,Plotly
using SpinFRGLattices
using FRGLatticePlotting.Plots
Plots.plotly()
##
Basis = Cif.getBasis("../test/CentredPyrochlore.vesta")
function isInBox(s::Rvec,Basis::Basis_Struct,L)
    x,y,z = getCartesian(s,Basis)
    return 0 <= x <= L && 0 <= y <= L && 0 <= z <= L
end
# sites = generatePairSites(4,Basis)
sites = generateLUnitCells(1,Basis)
filter!(s-> isInBox(s,Basis,1),sites)
pairsPlot(sites, Basis)

Bonds = Cif.readBondsVesta("../test/CentredPyrochlore.vesta")

# sites1 = filter(s-> getSiteType(s,Basis) == 1,sites)
# sites2 = filter(s-> getSiteType(s,Basis) == 2,sites)

plotDistBonds!(sites,Basis, minDist = Bonds[1].minDist, maxDist = Bonds[1].maxDist,lw = 10,color = :darkred)
plotDistBonds!(sites,Basis,minDist = Bonds[2].minDist, maxDist = Bonds[2].maxDist,lw = 10,color = :black)
##
Basis = Cif.getBasis("../test/CentredPyrochlore.vesta")
CPyro = Cif.generateSystem(6,"../test/CentredPyrochlore.vesta")
plotSystem(CPyro,Basis,refSite=1,bondDist = 0.35355,markersize = 3,plotAll = true)
plotSystem(CPyro,Basis,refSite=1,bondDist = 0.35355,markersize = 3,plotAll = true)