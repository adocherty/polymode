from pylab import *
from numpy import *
from Polymode import *

#Load previously calculated modes
modes = load_data('ex2_modes.dat')

#Plot modes with vector field
figure(1)
Plotter.plot_modes_in_grid(modes, 'sz', cmap=Plotter.cm.gray, Nx=(50,60), rmax=3)
Plotter.plot_modes_in_grid(modes, 'vectore')

figure(2)
subplot(121)
modes[0].plot('hx', Nx=(50,60))

subplot(122)
modes[0].plot('hy', Nx=(50,60))

#This command might be needed to show a graph
show()

#To update a plot use
draw()
