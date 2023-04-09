using FRGLatticePlotting,Plotly
using SpinFRGLattices
using FRGLatticePlotting.Plots
plotly()
##
a = Cif.getBasis("../test/CentredPyrochlore.vesta")
sites = generatePairSites(4,a)
# sites = generateLUnitCells(1,a)
pairsPlot(sites, a)
plotDistBonds!(sites,a, minDist = 0.35355 -1e-3, maxDist = 0.35355+1e-3,lw = 10)
plotDistBonds!(sites,a,lw = 4,color = :red)