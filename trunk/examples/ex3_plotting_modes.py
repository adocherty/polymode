from pylab import *
from numpy import *
from Polymode import *

#The modes
modes, wg = load_data("ex3_modes.dat")

#Sort the modes by reversed order of effective index
modes = sorted(modes, reverse=1)

#Print information about the modes
for mode in modes:
    mode.info()

#Create or select the first figure
figure(1)

#Plot the mode power and the H-polarization
Plotter.plot_modes_in_grid(modes, plottype='sz')
Plotter.plot_modes_in_grid(modes, plottype='ph', title='number')

#The second figure
figure(2)

#Clear it
clf()

#Select the second mode
m = modes[0]

#Plot Ez in a cartesian grid. Note that the waveguide is required to calculate Ez
subplot(121)
m.plot('ez', cartesian=1, rmax=10, wg=wg)
colorbar()                                      #A color bar gives the valueof the colors
wg.plot(fill=0)                         #Plot an outline of the waveguide on top

#Now plot Hz in the same way. We don't need the waveguide for this
subplot(122)
m.plot('hz', cartesian=1, rmax=10)
colorbar()                                      #A color bar gives the valueof the colors
wg.plot(fill=0)                         #Plot an outline of the waveguide on top

#To save a figure to a file use this command
savefig("fig_compressmodes.png")
