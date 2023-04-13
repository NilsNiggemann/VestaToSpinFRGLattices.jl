import CifToSpinFRGLattice as Cif
using FRGLatticePlotting,Plotly
using SpinFRGLattices
using FRGLatticePlotting.Plots
plotly()
##
Basis = Cif.getBasis("../test/CentredPyrochlore.vesta")
sites = generatePairSites(4,Basis)
# sites = generateLUnitCells(1,a)
pairsPlot(sites, Basis)
plotDistBonds!(sites,Basis, minDist = 0.35355 -1e-3, maxDist = 0.35355+1e-3,lw = 10)
plotDistBonds!(sites,Basis,lw = 4,color = :red)
##
CPyro = Cif.generateSystem(6,"../test/CentredPyrochlore.vesta")
plotSystem(CPyro,Basis,refSite=2,bondDist = 0.35355,markersize = 3)


##
plotDistBonds(CPyro,Basis;lw = 5)