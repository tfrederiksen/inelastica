from Inelastica import MiscMath as mm
import numpy as N
from Inelastica import WriteXMGR as XMGR
import scipy.special as ss

# We try a simple Gaussian function:
x0 = 50
x = N.linspace(-x0,x0,1e5)
gauss = -1j*N.pi**-0.5*N.exp(-x**2)

# Output from our Hilbert function:
Hg,ker = mm.Hilbert(gauss)

# We can compare with this:
# Hilbert transform of a Gaussian function is related to Faddeva/w(z) functions:
# https://en.wikipedia.org/wiki/Dawson_function
ex = -1j*N.pi**0.5*ss.wofz(x)

# Collect data in a plot
data1 = XMGR.XYset(x,gauss.imag,legend='Gauss',Lwidth=2)
data2 = XMGR.XYset(x,ex.real,legend='Faddeva',Lwidth=2)
data3 = XMGR.XYset(x,N.pi*Hg.imag,legend='Inelastica',Lwidth=1)
graph1 = XMGR.Graph(data1,data2,data3)
graph1.SetXaxis(autoscale=True)
graph1.SetYaxis(autoscale=True)
graph1.ShowLegend()

# Compute error/difference between the two methods:
err = ex.real-N.pi*Hg.imag
data4 = XMGR.XYset(x,err,legend='Difference',Lwidth=2)
graph2 = XMGR.Graph(data4)
graph2.SetXaxis(autoscale=True)
graph2.SetYaxis(autoscale=True)
graph2.ShowLegend()

plot = XMGR.Plot('test.agr',graph1)
plot.AddGraphs(graph2)
plot.ArrangeGraphs(nx=2,ny=1)
plot.WriteFile()

# Print measure of the difference
print 'Max deviation between the two evaulations:',max(abs(err))

