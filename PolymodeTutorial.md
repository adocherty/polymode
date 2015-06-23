![http://polymode.googlecode.com/svn/misc/images/polymode2.png](http://polymode.googlecode.com/svn/misc/images/polymode2.png)



## Running the tutorials ##
The code for all the tutorials can be downloaded from the "Downloads" section of the
[Polymode website](http://polymode.googlecode.com) or are also available in the "examples" directory
of the svn repository if you downloaded the code from svn.

The examples which contains a number of test programs that use the Polymode code. After installing Polymode
you can run the examples or your own code from any directory.

The **IPython** interactive python shell is assumed throughout this tutorial, although they
can be run directly using python non-interactively as well.
To run these examples open the IPython shell and cd to this directory.

```
cd _examples directory_
```

Run the  program (in this case named “ex1\_materials.py”) typing:

```
run _ex1_materials.py_
```

You can include the “.py” extension but it is not necessary.

These examples don't assume any familiarity with the Python programming language, however
for more advanced usage it is useful to have some knowledge of programming in it. A good
introduction to Python programming is available at the [Python tutorial](http://docs.python.org/tutorial/index.htm)
or the [Dive Into Python](http://www.diveintopython.org/) websites.

The numerical capabilities of Python are provided by the [Scipy and Numpy libraries](http://www.scipy.org), and more information on how to use Python for numerical calulations is given on their website. If you are unfamiliar with them a good place to get information is on their [Getting Started](http://scipy.org/Getting_Started) webpage.

The plotting is provided by [Matplotlib](http://matplotlib.sf.net) and a good introduction to commands available (which are all
similar to Matlab plotting commands) is available on their website.


---


## Tutorial 1 ##

### Materials ###
At the **IPython** interactive prompt run the  program (in this case named “ex1\_materials.py”) typing:

```
run  ex1_materials.py
```

You should get an output similar to this:

```
PMMA
Wavelength limits (um): 0.25->1.5
Temperature limits (K): 0->400
Cite: 
Refractive index at 1.45um is 1.48074+0i

Fused Silicon Oxide (SiO2)
Wavelength limits (um): 0.2->3.8
Temperature limits (K): 0->400
Cite: Malitson 1965, http://www.opticsinfobase.org/abstract.cfm?URI=josa-55-10-1205
Refractive index at 1.45um is 1.4452+0i

Silica doped with Germania
Wavelength limits (um): 0.6->1.8
Temperature limits (K): 0->400
Cite: Sunak, H.R.D.; Bastien, S.P., Photonics Technology Letters V1 N6 142-145, 1989
Refractive index at 1.45um is 1.45313+0i
```

The code itself is shown and explained following:

```
from numpy import *
from Polymode import *

# Basic materials 
air = Material.Air() 
m145 = Material.Fixed(1.45) 
```

The first two lines are needed in all programs using Polymode, they import the Polymode code for us to use and the numerical python code to enable mathematical operations on arrays.
The next two lines create two of the most common materials. `Material.Fixed` allows you to create a material
with a fixed refractive index, here 1.45, that doesn't change with wavelength.
`Material.Air` is the same as  `Material.Fixed(1)` - it gives a material with a fixed refractive index of 1.

```
silica = Material.Silica() 
germania = Material.Germania() 
pmma = Material.Polymer() 
```

These are some common materials for waveguides to be constructed from.
They have different refractive indices at different wavelengths which is automatically calculated.

```
sige = Material.SiO2GeO2(0.05) 
```

The inbuilt material `SiO2GeO2`, Silica doped with Germania , has a higher refractive index than silica and is used for the core of optical fibers. This approximation is accurate for a dopant concentration up to about 20%.

```
wl = 1.45

for mat in [pmma, silica, sige]: 
	ri = mat.index(wl) 
	mat.info() 
	mat.plot()
	print "Refractive index at %gμm is %s\n" % (wl, ri)
```

Here we iterate over the materials `pmma`, `silica` and `sige` and print their refractive index at a wavelength `wl` in micrometers. You can get some information about the material, including a reference for where the data is derived using the `.info()` method for any material. You are now at a command prompt for python and can enter commands and explore the materials that have been defined. The refractive index of the material versus the wavelength is plotted with the command `mat.plot()`

_Note that you may need the `draw()` and `show()` commands to refresh the plot window_

![http://polymode.googlecode.com/svn/wiki/figures/ex1_materials.png](http://polymode.googlecode.com/svn/wiki/figures/ex1_materials.png)

### The Waveguide ###
To construct an optical fiber we now need to create a waveguide and put different materials in it.

```
run ex1_the_waveguide
```

Try running the example, it should plot the waveguide and look like this:
![http://polymode.googlecode.com/svn/wiki/figures/ex2_waveguide.png](http://polymode.googlecode.com/svn/wiki/figures/ex2_waveguide.png)

The different colors represent different materials, here the green is the doped core and the surrounding grey the undoped silica.

Let’s see how this has been done, by looking at the code of `ex1_the_waveguide.py`.
There are three parts to creating the waveguide: choosing the waveguide parameters and materials, and then assembling the waveguide structure by placing particular shapes into it.

```
wg = Waveguide.Waveguide(material=cladding, symmetry=1) 
```

**Material**: The material that comprises the waveguide substrate.

**Symmetry**: To make the program more efficient we can exploit the rotational symmetry of the fiber. For example, if the waveguide has six identical sectors we only need to solve in a 60 degree segment and the symmetry is 6. A completely asymmetric structure would have a symmetry of 1. Lower symmetries are much slower to calculate.

```
s1 = Waveguide.Circle(core, center=(0,0), radius=2.0) 
wg.add_shape(s1) 
```

A waveguide needs to have some shapes added to it. There are a number of shaped available in the Waveguide module. The first parameter in any shape is the material it is made from, then there are different parameters for each shape type.

In this case we add the core as a circle at the origin with a radius of 2 μm. To plot the waveguide we now use the command:

```
wg.plot()
```

Tip: You can plot the refractive index profile that is used to numerically represent the structure by using the command `wg.plot_bitmap((N,,r,,,N,,ϕ,,))`
where `N,,r,,` is the radial resolution and `N,,ϕ,,` is the azimuthal resolution.


---


## Tutorial 2 ##

### Finding Modes of a Step Index Fiber ###

This exercise takes the waveguide constructed in the preceding problem and finds some supported modes.

Running the code by the command `run ex2_step_index_modes`.
Will construct the waveguide and solve for three modes giving the following output from the mode solver,
as the modes are solved:

```
I: Starting solve for 2 modes with neffs near 1.48183018032
I: NL Residual solver. Finding 'bound' modes, m0=0, wl=1
I: Mode #0 [15:0.77s] neff:1.46319 + 5.24e-18i res: 1.464e-12
I: Mode #1 [13:0.738s] neff:1.46303 + 6.18e-19i res: 1.11e-13

I: Starting solve for 1 modes with neffs near 1.48183018032
I: NL Residual solver. Finding 'bound' modes, m0=1, wl=1
I: Mode #0 [12:0.61s] neff:1.47428 - 2.33e-18i res: 2.1e-13
```

The information printed during the solving process gives information on the current mode found:

Mode _#number_ `[`_iterations_:_time_`]`, neff=_effective index_, res: _approx error_

The iterations information is how many iterations is needed to refine the leakage loss of the mode to the stated accuracy; highly leaky modes require more iterations. The value of the residue gives an indication of the accuracy of the modes - clearly neff should not be quoted to more significant figures than indicated the residue.

The code of `ex2_step_index_modes` is the same as the previous example where we created a step index fiber,
however we now create a solver with the specified waveguide and resolution.

```
Nx = 100,1	#number of radial & azimuthal points
wl = 1.0    #Wavelength

solver = NLSolver.DefaultSolver(wg, Nx)
modes = solver(wl, 0, number=2)
modes += solver(wl, 1, number=1)
```

To solve for modes in our created waveguide we need to give the resolution of the numerical solution Nx. This is a tuple of two numbers, (N<sub>r</sub>,N<sub>ϕ</sub>)) where N<sub>r</sub> is the radial resolution and N<sub>ϕ</sub> is the azimuthal resolution.

Note that the number of radial and azimuthal points defines the resolution of the calculation. Increasing the number of radial points linearly increases the time taken for the calculation, but the increase in calculation time as you increase the azimuthal resolution is much more dramatic. That is unfortunate as features like thin bridges with require a lot of points to resolve them properly. You should always check for proper convergence of the solution by running at different resolutions (for example, doubling the number of points). You need to be even more careful about convergence when considering the imaginary part of n<sub>eff</sub>.

The solver object is created with the wanted resolution and the waveguide as arguments. This solver can then be used to calculate any number of modes calling it with the wavenumber and  the symmetry class. By default the solver looks for modes with effective index near the highest refractive index in your structure; this may not always be the modes you are interested in, and different ways to specify which modes to look for will be detailed later.

The modes are returned as a list, if we just type `modes` and hit enter the modes in the list are printed.
Any mode in the list can be inspected interactively, let's investigate the first mode. It has some properties, the effective index, the loss, the symmetry class, the wavelength. Try these commands:
```
print modes[0].wl
print modes[0].neff
print modes[0].loss
print modes[0].m0
```

They will print the wavelength, effective index, loss and symmetry mode class of the first mode in the list respectively.

To go through each mode and print information on it we can use this command:

```
for mode in modes: 
	mode.info()
```

The output is:

```
VectorMode size: (100, 1), symmetry=1, m0=0, λ=1, rmax=2.1, vecs: RL {} 
 | neff=1.46319 + 6.10e-19i 	loss=3.327e-11dB/m 	res=2.58e-13 
 | Effective Area: 10.111 μm² 
 | Spot size: 0.4035 μm 
 | Numerical Aperture: 0.1747 

VectorMode size: (100, 1), symmetry=1, m0=0, λ=1, rmax=2.1, vecs: RL {} 
 | neff=1.46303 - 4.34e-19i 	loss=-2.368e-11dB/m 	res=1.47e-13 
 | Effective Area: 10.076 μm² 
 | Spot size: 0.3979 μm 
 | Numerical Aperture: 0.175 

VectorMode size: (100, 1), symmetry=1, m0=1, λ=1, rmax=2.1, vecs: RL {} 
 | neff=1.47428 - 7.99e-19i 	loss=-4.358e-11dB/m 	res=4.78e-12 
 | Effective Area: 5.7825 μm² 
 | Spot size: 0.7165 μm 
 | Numerical Aperture: 0.22842 
```

### Leaky Modes of a Step Index Fiber ###

Run this tutorial as usual by `run ex2_leaky_modes`.
This should print the list of found modes like the previous tutorial and then plot them.

```
I: Starting solve for 2 modes with neffs near 1.48183018032
I: NL Residual solver. Finding 'outward' modes, m0=0, wl=1
I: Mode #0 [24:1.31s] neff:1.46314 + 7.89e-04i res: 1.698e-11
I: Mode #1 [20:1.15s] neff:1.46301 + 8.07e-04i res: 5.775e-14

I: Starting solve for 1 modes with neffs near 1.48183018032
I: NL Residual solver. Finding 'outward' modes, m0=1, wl=1
I: Mode #0 [13:0.703s] neff:1.47439 + 1.69e-04i res: 1.101e-13
```

Look at this program file, the program is the same as the previous apart from the waveguide now has a substrate of polymer and an annlus of air that forms the cladding. Plot the refractive index profile of the waveguide:

```
wg.plot_bitmap(Nx)
```

This structure gives rise to modes that are guided by the core and cladding, but energy from them can tunnel through the this air gap and this causes leakage loss. This loss is found by the numerical solver as the imaginary part of the effective index.

```
solver = NLSolver.DefaultSolver(wg, Nx) 
modes = solver(wl, 0, number=2) 
modes += solver(wl, 0, number=1) 
```

We find modes in symmetry classes 0 and 1 just as previously. However the modes will now be lossy. The loss of each mode can be shown using the info() function:

```
modes[0].info() 

VectorMode size: (100, 1), symmetry=1, m0=0, λ=1, rmax=3.15, vecs: RL {} 
 | neff=1.46314 + 7.89e-04i 	loss=4.305e+04dB/m 	res=4.62e-12 
 | Effective Area: 10.618 μm² 
 | Spot size: 0.73925 μm 
 | Numerical Aperture: 0.17061
```

The entire mode list is now plotted with the command:

`Plotter.plot_modes_in_grid(modes, rmax=5)`

To plot a single mode we can use the command:

`modes[mode number].plot()`

Where the mode number is a number from 0 to _number of modes_ – 1.

The plot will default to plotting the power, (the z-component of the time-averaged Poynting vector). However other parts of the electromagnetic field can also be plotted. For example the command
`modes[0].plot('Hr')`
plots the radial component of the magnetic fieldof the first mode.

Other parts of the electromagnetic field that can be plotted are:

  * 'Sz' - z-component of the Poynting vector (power flow)
  * 'Hr' - radial-component of the magnetic field
  * 'Ha' - azimuthal-component of the magnetic field
  * 'Er' - radial-component of the electric field
  * 'Ea' - azimuthal-component of the electric field
  * 'Hz' - z-component of the magnetic field
  * 'Ez' - z-component of the electric field
  * 'Hx' - cartesian x-component of the magnetic field
  * 'Hy' - cartesian y-component of the magnetic field
  * 'Ex' - cartesian x-component of the electric field
  * 'Ey' - cartesian y-component of the electric field
  * 'pe' - polarization of the electric field
  * 'ph' - polarization of the magnetic field

The default is to plot the z-component of the Poynting vector.

This example used the save\_data and load\_data commands to save calculated modes and load them back again to plot them. They can be used to save any objects you like including waveguides, materials, shapes and modes.

```
save_data(data, filename)
data = load_data(filename)
```

Try running the next file `ex2_plotting_modes.py` to see some commands in action.


---


## Tutorial 3 ##

### Microstructured fibers ###
This tutorial file is `ex3_microstructured_fiber.py`.

A microstructured fiber can be created in the same way as the fibers we have created so far, however we will now have much more flexibility of what shapes to create in the waveguide.
Initially let's create a waveguide from PMMA polymer with six circular airholes. The structure we want looks like this:

![http://polymode.googlecode.com/svn/wiki/figures/ex3_microstructured_fiber.png](http://polymode.googlecode.com/svn/wiki/figures/ex3_microstructured_fiber.png)

The waveguide can be created in the same way, but now we use the fact that it is six-fold symmetric and set the symmetry to be six. We then only need to add one circule and the rest will be created automatically by the symmetry. The code to do this is:

```
air = Material.Air() 
polymer = Material.Polymer() 

wg = Waveguide.Waveguide(material=polymer, symmetry=6) 
s1 = Waveguide.Circle(air, center=(3,pi/6), radius=1.25) 
wg.add_shape(s1) 

wg.plot()
```

Try plotting the waveguide with the resolution N<sub>r</sub>=100, N<sub>ɸ</sub>=20 using the `plot_bitmap` method.

There are many other shapes that can be used to create different geometries. Some of them with their options are:

  * Waveguide.Circle( material, center=(r<sub>0</sub>, ɸ<sub>0</sub>), radius=radius, zorder=0, xy=False)
  * Waveguide.Ellipse( material, center=(r<sub>0</sub>, ɸ<sub>0</sub>), axes=(a,b), rot=0, zorder=0, xy=False)
  * Waveguide.Rectangle( material, center=(r<sub>0</sub>, ɸ<sub>0</sub>), axes=(a,b), rot=0, zorder=0, xy=False)
  * Waveguide.Polygon( material, nodes=`[` _list of node coords_ `]`,  zorder=0, xy=False)
  * Waveguide.Annulus( material,  r=(r<sub>min</sub>, r<sub>max</sub>, phi=(ɸ<sub>min</sub>, ɸ<sub>max</sub>), zorder=0)
  * Waveguide.Image( material, image=None, border=0)

The axes parameters for the ellipse and rectangle set the width and height of the objects relative to the coordinate system. They will be aligned along the radial and polar directions by default, however giving an argument to `rot` will change this.

It is important  to note that the coordinate system defaults to polar, however setting the option xy to True enables the center to be given in cartesian coordinates, and the major and minor axes or ellipses and rectangles to be in the cartesian directions.


### Solving Modes in Microstructured Fibers ###
This tutorial file is `ex3_finding_modes.py`.

We now want to find some modes in the previous fiber geometry. As before we create a solver object and run the solver with different symmetry indices to get different classes of modes. Note that the symmetry index with a fiber of symmetry N will range from -N<sub>sym</sub>/2 … N<sub>sym</sub>/2 if N<sub>sym</sub> is even and -(N<sub>sym</sub>-1)/2 … (N<sub>sym</sub>+1)/2 if N<sub>sym</sub> is odd. The modes for negative and positive indices are degenerate. We will come back to this point.

The code to solve for modes in the waveguide in `ex3_finding_modes.py`is:

```
Nx = 100,11
wl = 1.0

# Create the solver
solver = NLSolver.DefaultSolver(wg, Nx)

#Solve for modes in different symmetry classes.
#Notice the "modes +=" adds to the modes
modes = solver(wl, 0, 1.45, number=2)
modes += solver(wl, 1, 1.45, number=2)

# Save the modes and the waveguide
save_data((modes, wg), "ex3_modes.dat")

#Create a new plot window
figure()

Plotter.plot_modes_in_grid(modes)
```

We now specify the number of radial and azimuthal points. We solver for 2 modes with a symmetry index of 0
and three for a symmetry index of 1. The modes as well as the waveguide are saved to the file 'ex3\_modes.dat',
note the use of brackets to group them together. The modes are plotted in a grid as before.

![http://polymode.googlecode.com/svn/wiki/figures/ex3_finding_modes.png](http://polymode.googlecode.com/svn/wiki/figures/ex3_finding_modes.png)

### Hexagonal lattice ###
This tutorial file is `ex3_hexagonal_lattice.py`.

The `Waveguide` module contains a number of layout functions that can put shapes
in common configurations such as hexagonal and square lattices. To create a hexagonal
lattice we choose the symmetry of the waveguide, the number of layers and the layer
spacing. You can also specify the offset of the layers, which is the number of layers
to remove at the center to create the 'defect'.

```
sym = 6
layers = 2
D = 6.0
d = 4.0

circle = Waveguide.Circle(air, radius=d/2)
shapes = Waveguide.create_hexagonal_tiling(circle, layers=layers, D=D, symmetry=sym)
wg.add_shapes(shapes)

solver = NLSolver.DefaultSolver(wg, Nx)
modes = solver(wl, 0, nefflist=[1.48, 1.475,1.46])

Plotter.figure()
Plotter.plot_modes_in_grid(modes, wg=wg)
Plotter.show()
```

Note here we use the `nefflist` option to the solver to solve for modes near to these
approximate effective indices.


### Square lattice ###
This tutorial file is `ex3_square_lattice.py`.

Waveguide can be created programmatically which allows great ease
and flexibility in created different waveguides. As an example we
create a square lattice waveguide with a central lattice element removed.

We create this using a `for` loop to allow a variable number of rings of holes
in the lattice. First we set up the standard parameters and materials. Notice
the summetry of this waveguide is 4.
```
sym = 4
layers = 2
D = 5.0
d = 0.5*D

c = Waveguide.Circle(air, radius=d/2)
shapes = Waveguide.create_square_tiling(c, layers=layers, Dx=D, Dy=D, symmetry=sym)
wg.add_shapes(shapes)

#Find modes
solver = NLSolver.DefaultSolver(wg, Nx)
modes = solver(wl, 0, neffrange=[1.4,1.45], number=2)
modes += solver(wl, 1, neffrange=[1.4,1.45], number=2)
modes += solver(wl, 2, neffrange=[1.4,1.45], number=2)
```

Here we use the `neffrange` argument to the solver. This finds `number` modes within the effective index
range specified by `neffrange`.

The waveguide is shown below. The sector we have created is in black and the
gray holes are created by the automatic rotational symmetry.

![http://polymode.googlecode.com/svn/wiki/figures/ex3_square_lattice_wg.png](http://polymode.googlecode.com/svn/wiki/figures/ex3_square_lattice_wg.png)

### Spurious modes ###
This tutorial file is `ex3_square_lattice.py`.

The solver will try and find all modes within the specified range
but may not be able to do so, especially if there are not many modes supported
in the structure or there are many spurious modes then the solver may
not converge. This is a limitation of the current solver and will
hopefully be improved soon, meanwhile don't be alarmed if
_Rejecting unconverged mode_ is printed while trying to find modes.

Even when a mode is found it is possible that it is not a real physical
mode of the structure but a numerical artifact. These _spurious_ modes
have the majority of their power trapped between the structure and the
computational boundary. They are easy to spot and the solver tries to identify them
automatically.

Typing `modes` to get a list of the modes solved we see some modes marked with an
`[S]` at the end. These are the ones the solver has picked are spurious.
```
 <VectorMode: m0=0 wl=1 neff=1.44663 + 2.66e-03i r:6.3e-11  [S]>
 <VectorMode: m0=1 wl=1 neff=1.44816 + 4.44e-09i r:1e-10  []>
 <VectorMode: m0=1 wl=1 neff=1.4472 + 2.36e-03i r:2.7e-13  [S]>
 <VectorMode: m0=2 wl=1 neff=1.45018 + 2.37e-03i r:3.2e-11  [S]>
```
In this list all but the second mode are spurious. We can also see that the other three
modes have very large imaginary parts of the effective indes, or large losses.
This is due to them being artificially guided between the outside of the structure
and the inexact boundary conditions. This structure has a large number of such modes as
it has large gaps between the structure and the computational boundary.

When we plot them we see they're highlighted with a red background.

![http://polymode.googlecode.com/svn/wiki/figures/ex3_square_lattice.png](http://polymode.googlecode.com/svn/wiki/figures/ex3_square_lattice.png)

After manually examining them we can automatically remove the spurious modes
using the `Modes.filter_suprious_modes()` function. This removes all
spurious modes marked as such manually or automatically.

```
modes = Modes.filter_suprious_modes(modes)
```

If you don't agree with the automatic labelling of a mode as spurious or non spurious
you can manually change it by setting the `.is_spurious` property of the mode to
`True` or `False`.
For example we can set the 10th mode to be spurious by the command:
```
modes[10].is_spurious = True
```


### Plotting modes continued ###
This tutorial file is `ex3_plotting_modes.py`.

Here we use some more commands to plot modes in different ways. First we load the modes
and waveguide we solved for in the previous tutorial:
```
modes, wg = load_data("ex3_modes.dat")
```

And sort the modes by reversed order of effective index
```
modes = sorted(modes, reverse=1)
```

We use some plotting commants provided by Matplotlib (the Python plotting library) which are very similar
to Matlab plotting commands.

  * `figure(n)` - creates or selects the _nth_ plot window.
  * `clf()` - clears the current figure.
  * `subplot(nmj)` - creates a subplot figure with _n_ rows, _m_ columns and selects the _jth_ one.
  * `colorbar()` - creates and shows a colour scale for the plot.

There are many ways to plot the modes. Using the `Plotter.plot_modes_in_grid` also takes the `plottype` argument which is the same as the individual mode plot command, namely `'Sz', 'Hr', 'Ha', 'Er', 'Ea', 'Hz', 'Ez', 'Hx', 'Hy', 'Ex', 'Ey','pe', 'ph'`.

Issuing two different plot commands succesively will plot one plot over another:
```
Plotter.plot_modes_in_grid(modes, plottype='sz')
Plotter.plot_modes_in_grid(modes, plottype='ph', title='number')
```

The z-component of the electric and magnetic field is calculated from the transverse components,
therefore to calculate the z-component of the electric field the waveguide needs to be specified.

```
#Select the second mode
m = modes[1]

#Plot Ez in a cartesian grid. Note that the waveguide is required to calculate Ez
subplot(121)
m.plot('ez', cartesian=1, rmax=10, wg=wg)
colorbar()					#A color bar gives the valueof the colors
wg.plot(fill=0)				#Plot an outline of the waveguide on top

#Now plot Hz in the same way. We don't need the waveguide for this
subplot(122)
m.plot('hz', cartesian=1, rmax=10)
colorbar()					#A color bar gives the valueof the colors
wg.plot(fill=0)				#Plot an outline of the waveguide on top
```

Here we also plot the modes in caresian coordinates using the `cartesian` flag.

You can save a figure using the plot GUI or the `savefig` command
```
savefig("fig_compressmodes.png")
```

![http://polymode.googlecode.com/svn/wiki/figures/ex3_plotting_modes.png](http://polymode.googlecode.com/svn/wiki/figures/ex3_plotting_modes.png)

Documentation for more general plotting commands can be found at the
[Matplotlib website](http://matplotlib.sourceforge.net).



---


## Tutorial 4 - More advanced waveguides ##

### More shapes ###
The example `ex4_more_shapes.py` creates a waveguide with three different shapes and finds the fundamental mode in it.

![http://polymode.googlecode.com/svn/wiki/figures/ex4_more_shapes.png](http://polymode.googlecode.com/svn/wiki/figures/ex4_more_shapes.png)


### Coated holes ###
This tutorial file is `ex4_coated_holes.py`.

Shapes can be coated with another material of a specified thickness using the
```
    Material.coated_shape(shape, coating material, thickness)
```
function. This takes a previously created shape, a coating material and a thickness and returns a shape to use as the coating. Note that any shape can be coated, not just circular inclusions.
This is seen in the `ex4_coated_holes.py` file. The relavant portion is shown here, a circle is created and coated with a material of higher refractive index and both shapes are added to the waveguide.

```
hole = Waveguide.Circle(air, center=(6.75,0), radius=1.5) 

#Create a coating of material m17 from hole with a coating thickness of 0.2um
chole = Waveguide.coated_shape(hole, m17, 0.2)

#Add hole and coating to waveguide
wg.add_shape(chole)
```

The fundamental mode is found:

![http://polymode.googlecode.com/svn/wiki/figures/ex4_coated_holes.png](http://polymode.googlecode.com/svn/wiki/figures/ex4_coated_holes.png)

### Waveguide from a bitmap image ###
This tutorial file is `ex4_waveguide_from_image.py`.

The second way to enter a structure is to use an image. The example folder gives images of waveguide structures this way **MOF\_NRL.png**, an example taken from literature of a photonic crysal filber and and **hlfcollapsedKagome.png** a schematic of a kagome structure. Pictures of these structures are found in the directory called images.

Images are used in a similar way to other objects. However first the bitmap must be imported using the `Image.ImportBitmap(filename, size=(sx,sy), negative=True/False)` function.
This imports an image from the named file of a specified x and y dimentions in microns. The negative flag is false if the bright parts of the image are to be considered the material and false if the opposite is the case.
Then the waveguide shape can be created using the command `Waveguide.Image(material, image)`

The image will take a while to be converted to work with the waveguide as the image is optimized for the waveguide symmetry and converted to polar coordinates.

Selected code from the `ex4_waveguide_from_image.py` example file is explained following:

```
image_filename1 = os.path.join("images","MOF_NRL.png")
image_filename2 = os.path.join("images","hlfcollapsedKagome.png")
```
Note that the `os.path.join` function joins path names in a platform independent way so it will work on Linux and Windows. On Windows you can equivalently just write "images\MOF\_NRL.png".

```
image = Image.ImportBitmap(image_filename1, size=(20,20), offset=pi/6, negative=True)
image.contrast(scale=15, transition=0.45)
wg.add_shape(Waveguide.Image(air, image))
```
To load an image we use the `ImportBitmap` object to load the bitmap. The required arguments are the size of the
entire image in micrometers, the offset (to adjust the angle of the loaded bitmap appropriately) and weather
the lightest colors should be interpreted as the substrate or as the shape material. If `negative` is False (the default)
the lightest colors are converted to the shape's material, the darkest are assigned the backgrounsd index, for
the case that `negative` is True the opposite holdes.

Other actions can be applied to the image object, it may be necessary to apply some contrast to the image, in this case we adjust the contrast by changing the amount of contrast (scale) and the relative light/dark bias (transition).

The `plot_bitmap` member function gives the actual numerical values of the waveguide on a polar grid of
200x60 pixels (in this case), the sectors=1 specified plotting just one of the symmetric sectors and the
cmap command chooses the colormap.

```
wg.plot_bitmap((200,60), sectors=1, cmap=Plotter.cm.gray)
```

![http://polymode.googlecode.com/svn/wiki/figures/ex4_waveguide_from_image.png](http://polymode.googlecode.com/svn/wiki/figures/ex4_waveguide_from_image.png)


### Merging Shapes ###
This tutorial file is `ex4_merging_shapes.py`.

If shapes created in the waveguide overlap then they are automatically merged and one will override the other, appearing "in front" of the other. This behaviour can be changed with the shape property `zorder` - this defaults to 0 and the higher zorders appear in front of the lower zorders.

The code of  `ex4_merging_shapes.py` specifies the `zorder` of an shape using the `zorder=1` argument when creating the shape, and then changes the orders directly accessing the `shape.zorder` property.

```
s0 = Waveguide.Ellipse(m2, center=(0,0), axes=(3.0,1.5))
s1 = Waveguide.Circle(air, center=(3,pi/6), radius=1.5, zorder=1)
wg.add_shapes(s0, s1)

Plotter.subplot(121)
wg.plot()

#Change orders around
s0.zorder=1
s1.zorder=0

Plotter.subplot(122)
wg.plot()
```

The waveguide on the left should have the blue circles in front of the red ellipses.

![http://polymode.googlecode.com/svn/wiki/figures/ex4_merging_shapes.png](http://polymode.googlecode.com/svn/wiki/figures/ex4_merging_shapes.png)

Note the oscillations near the sharp changes in refractive index. This is the [Gibb's effect](http://en.wikipedia.org/wiki/Gibbs_effect) caused by the truncation of the Fourier series that is used internally.

### Refractive index profiles ###
This tutorial file is `ex4_parabolic_index.py`.

A non-step dependence on the refractive index of any shape can easily be defined.
The file `ex4_parabolic_index.py` creates a circular core and assigns it a parabolic refractive index.
A refractive index profile can be specified as a function of the distance from the objects center.

```
fn_parabolic = lambda d: 1-(d/radius)**2
```

This uses a Python “lambda function” to specify a function of d, distance. The function should be normalized to 0 and 1, a value of 1 will be assigned the refractive index of the object and 0 will be assigned the refractive index of the background. The defined function is 1 at the center of the circle and 0 when d=radius. The distance d is given to the function as the distance from the center of the object.

Any shape can be set to have an index function using the `set_index_function` method.

```
silica = Material.Silica() 
dsilica = Material.SiO2GeO2(0.2)
```

Here we create the waveguide with a given maximum computational radius, specifying `rmax=8`. If this is not given it is automatically assigned as the largest needed to enclose all the shapes in the waveguide, but sometimes you are required to specify it.
```
wg = Waveguide.Waveguide(rmax=8, material=silica, symmetry=1) 
```

Now we create the doped core and assign the index function to it
```
fn_parabolic = lambda d: 1-(d/radius)**2 

## Doped core 
s1 = Waveguide.Circle(dsilica, center=(0,0), radius=radius)
s1.set_index_function(fn_parabolic)

## Finally add the shape to the waveguide
wg.add_shape(s1) 
```

We can plot the actual refractive index profile that is used internally in the program using the `plot_bitmap` command. Try just using the `wg.plot` command and see what it plots.

```
wg.plot_bitmap(Nx, wl=1.0) 
```

The output of this example is shown below, the refractive index of Silica (1.4504) and Silica doped with 20% Germania (1.482) at a wavelength of λ=1.0μm is also plotted.

![http://polymode.googlecode.com/svn/wiki/figures/ex4_parabolic_index.png](http://polymode.googlecode.com/svn/wiki/figures/ex4_parabolic_index.png)


### Solving modes in an image waveguide ###
This tutorial file is `ex4_modes_with_a_image_waveguide.py`.

Here we use another image load it and solver 4 modes in it. The image does not have rotational symmetry so we
specify `symmetry=1`. The image is given a specified center as the automatic center calculation doesn't work without rotational symmetry.

```
# Create waveguide and specify the problem domain to have radius 15
wg = Waveguide.Waveguide(rmax=15, material=polymer, interior=polymer, symmetry=1)

# Load a png file and specify the size of the original image
# The center of the image (pixel offset) must also be specified for symmetry=1
image = Image.ImportBitmap("images/celery.png", size=(60,120), center=(230,140), negative=True)

# Add the shape to the waveguide
wg.add_shape(Waveguide.Image(air, image, border=2))

# Waveguide Parameters:
Nx=100,81		#number of radial & azimuthal points
wl=1.55			#Wavelength

# Create the solver
solver = NLSolver.DefaultSolver(wg, Nx)
modes = solver(wl, 0, number=4)

Plotter.plot_modes_in_grid(modes)
```

Note the `interior=polymer` is only typically needed for an image waveguide as the code cannot determine what the material of the core is. If the core material or material outside the calculation range (beyond rmax) is different from the waveguide's material then it needs to be specified:

```
Waveguide.Waveguide(material=_material_, interior=_material_, exterior=_material_)
```

Loading an image only allows you to have two different materials in the waveguide - the waveguide itself and the material assigned to the image. However, it is possible to add mode shapes to the waveguide on top of the image with different materials, even other images.

Be aware that using waveguides without symmetry will be very slow to solve for modes as the number of azimuthal points needed is large.

http://polymode.googlecode.com/svn/wiki/figures/ex4_modes_with_a_image_waveguide


---


## Tutorial 5 - Advanced solving ##

### Mode Quantities ###
This tutorial file is `ex5_modal_quantities.py`.

The `mode` object has a lot of information about the mode. The basic properties are:
  * The wavelength
  * The loss
  * The symmetry class
  * The propagation constant
  * The effective index

Other quantities can be calculated from the
  * Group index
  * Proportion of power in the core
  * Spot size
  * Nonlinear effective area, A<sub>eff</sub>

Some of these are accessed in the tutorial file
```
for m in modes:
	m.normalize(by='ext')

	print
	print "Mode effective index:", m.neff
	print "Loss: %.3g dB" % m.loss
	print "Symmetry class:", m.m0
	print "Propagation constant:", m.beta
	print "Wavelength:", m.wavelength

	print "Group index:", m.group_index(wg)
	print "Propagation contant from integral:", m.integral_propagation_lossless(wg)
	print "Proportion of power in core:", real(m.mode_power(r=wg.guess_core_size())/m.mode_power())
```



### Convergence ###
This tutorial file is `ex5_convergence_test.py`

It can be useful to check the accuracy of the calculated effective index or other modal properties. This can be done by running the solver with different numbers of points. The file `ex5_convergence_test.py` does just that, and presents the results.

Here instead of creating a solver and the running it immediately we create the solver and initialize it with the same arguments we would use to solve for the modes we want and then add the solver to a list.

This list of solvers is given to the `batch_solve` function which solves them sequentially, and writes the intermediate results to the filename specified. This means that if the calculation is interrupted then the results to that point can still be recovered.

```
wl = 1.45
Naz = 41
solvers = []
for Nr in [50,100,200,400,800]:
	solver = NLSolver.DefaultSolver((Nr,Naz), wg)
	solver.initialize(wl, 1, 1.45, number=1)
	solvers.append(solver)

modes = Solver.batch_solve(solvers, filename="convergencerun.solve")
```


```
#Extract neff and Nr
Nrs = [m.shape[0] for m in modes]
neff = [m.neff for m in modes]
```

To plot the convergence we use the most accurate run as an approximate to the true solution and plot the convergence versus the number of radial points. Python list syntax is used here, the basics are that the members of a list can be accessed using square brackets. If `x` is the list created by the command

```
x = [1,2,3,4]
```

Then `x[n]` will access the _n\_th element of the list, and if `n=-1` the_last_element is returned. So `x[0]` will return 1 and `x[-1]` will return 4. Sublists can be returned using a_slice

```
x[start:end:step]
```

which returns a new list composed of all elements from _start_ to _end_ in increments of _step_ - not including the element at _end_. So `x[1:-1]` is all elements in the list except the first and last, `[2,3]` in this case.
More information is available here  http://diveintopython.org/native_data_types/lists.html

The approximate error in the effective indices (except the final one used as an estimate to the true solution) are now plotted against the number of radial points with red circles 'ro'.

```
plot(Nrs[:-1], abs(neff [:-1] - neff[-1]), 'ro')
loglog()
xlabel(r'$N_r$')
ylabel(r'Error in $n_{eff}$')
```

We then change the axes to be logarithmic and put some labels on the axes. Note the use of latex in the axis labels.

![http://polymode.googlecode.com/svn/wiki/figures/ex5_convergence_test.png](http://polymode.googlecode.com/svn/wiki/figures/ex5_convergence_test.png)


### Eliptical polarization and degenerate modes ###
This tutorial file is `ex5_uranus_fiber.py`.

This tutorial calculates the modes in a three-fold symmetric fiber from the paper by
[Uranus and Hoekstra](http://www.opticsinfobase.org/oe/abstract.cfm?URI=oe-12-12-2795).
The waveguide is constructed of a single annulus with a span of 108 degrees. Note we convert
this to radians for the `Waveguide.Annulus` shape.

```
r1 = 1
r2 = 2
Daz = 108*pi/180
ring = Waveguide.Annulus(air, r=(r1,r2), phi=(-Daz/2,Daz/2))
wg = Waveguide.Waveguide(material=silica, symmetry=3)
wg.add_shape(ring)
```

We then find modes near the effective index 1.4 in mode classes 0 and 1. The modes and the electric polarization
are then plotted in a grid. Note we have sorted the modes in reverse order. Compare the effective indices to the modes in the Uranus paper. Increasing the resolution should improve the fit to the published effective indices.

```
vwe = NLSolver.DefaultSolver(wg, Nx)
modes = vwe(wl, 0, 1.4, number=2)
modes += vwe(wl, 1, 1.4, number=2)

modes.sort(reverse=True)

Plotter.plot_modes_in_grid(modes, 'sz,pe')
```

![http://polymode.googlecode.com/svn/wiki/figures/ex5_uranus_fiber.png](http://polymode.googlecode.com/svn/wiki/figures/ex5_uranus_fiber.png)

The modes in mode class 0 are linearly polarized however those in mode class 1 have a circular polarization.
The modes in class 1 are degenerate modes, however the solver will find a single circularly polarized mode that is a combination of both of them.
Compare this to figure 5 in Uranus and notice that the modes with circular polarization in the above figure correspond to the degenerate mode pairs `HE11(a,b)` and `HE21(a,b)`.

We'll look at extracting two degenerate linearly polarized modes from a circularly polarized mode in the next tutorial.


### Solving over a wavelength range ###
This tutorial file is `ex5_wavelength_scan.py`

The wavelength properties of a mode or modes can be evaluated by calculating the mode over a range of wavelengths.
There are two basic solvers for solving at a range of wavelengths, they are
`Solver.WavelengthScan` and `Solver.WavelengthTrack`.

The first one finds all modes with given search criteria at each wavelength. This can be quite slow
and there is no guarantee that they same modes found at one wavelength will be found at the next wavelength
step. If you are interested in the behaviour of many modes in the fiber and do not need the properties
of specific modes this may be useful.

The `Solver.WavelengthTrack` tracks specific modes across a range of wavelengths. This uses automatically selected
wavelength steps to speed the process. Here we show how to use this solver.

```
#Choose solver
solver = NLSolver.DefaultSolver(wg, Nx)

#Create wavelength track solver
wlsolver = Solver.WavelengthTrack(solver)

wlrange=[1.0,1.5]
modes = wlsolver(wlrange, 0, nefflist=[1.3])
```

Here the `Solver.WavelengthTrack` solver is created with the numerical mode solver as an argument.
It is then called with the wavelength range specified as the minimum and maximum wavelengths in the
range of interest.

```
subplot(121)
Plotter.plot_mode_properties(modes, 'neff', 'wl')
subplot(122)
Plotter.plot_mode_properties(modes, 'loss', 'wl')
```

Finally the effective index and loss of the modes versus the wavelength, in this case for the HE11-like mode. Try plotting the wavelength parameters for other types of modes too.

![http://polymode.googlecode.com/svn/wiki/figures/ex5_wavelength_scan.png](http://polymode.googlecode.com/svn/wiki/figures/ex5_wavelength_scan.png)



---


## Tutorial 6 - Calculating with modes ##

### Coordinates ###
This tutorial file is `ex6_coordinates.py`

To get the modal fields and perform calculations on them we need to evaluate the fields
and perform the calculations with a coordinate object.

There are two types that come with **Polymode**, they are
  * `coordinate.PolarCoord`
  * `coordinate.CartesianCoord`

They are created with the parameters
  * `coordinates.PolarCoord`(rrange=(r<sub>min</sub>, r<sub>min</sub>), arange=(a<sub>min</sub>, a<sub>max</sub>), N=_number of points_, border=_border_)
  * `coordinates.CartesianCoord`(X=±x<sub>max</sub>, Y=±y<sub>max</sub>,, N=_number of points_)

The number of points can also be specified for both axes as N=(N<sub>r</sub>, N<sub>a</sub>) for polar and  N=(N<sub>x</sub>, N<sub>y</sub>) for cartesian objects. The _border_ option forces the nodes to not go to the full minimum and maximum r. This should usually be set to `border=1` to avoid a zero minimum r.

The most important member functions of the coordinate objects are `.polar2d()` and `.cartesian2d()`
which return a 2d array of radial and azimuthal and x and y locations respectively of the node points in
cartesian object.

Other useful member functions are:
  * `.grad_t(f)` calculates the 2d vector gradiant the scalar data f
  * `.div_t(A)` calculates the 2d divergence of the vector data A
  * `.curl_t(A)`  calculates the z-component of the transverse vector data A
  * `.int_dA(f)`  calculates the area integral of the scalar data f

In the tutorial file we create two coordinate objects one for polar and one cartesian. The border is set to 1 for polar coordinates so we avoid problems with division by zero at the origin.

```
res = 200

#Coordinate objects
c1 = coordinates.PolarCoord(rrange=(0,2), border=1, arange=(-pi,pi), N=(res/2,res))
c2 = coordinates.CartesianCoord(X=2, Y=2, N=res)
```

Now the cartesian node locations are obtained for each grid:
```
x1,y1 = c1.cartesian2d()
x2,y2 = c2.cartesian2d()
```

And they are plotted using the _scatter_ function.
```
#Plot the point distribution
subplot(221)
scatter(x1.ravel(),y1.ravel(), 2)
axis('tight')

subplot(222)
scatter(x2.ravel(),y2.ravel(), 2)
axis('tight')
```

We have a function that is defined in cartesian coordinates. This could be also defined in radial coordinates as well.
```
f = lambda x,y: y*exp(-2*(x**2 + y**2))
```

The div(grad(f)) scalar function is calculated in both coordinates, this should be the same! We then plot them both in different subplots using the `Plotter.plot_v(coord, f)` function
```
g1=c1.div_t(c1.grad_t(f(x1,y1)))
g2=c2.div_t(c2.grad_t(f(x2,y2)))

subplot(223)
Plotter.plot_v(c1,real(g1))
axis('tight')

subplot(224)
Plotter.plot_v(c2,real(g2))
axis('tight')
```

We evaluate the integral of the absolute value of _f_ again with both coordinates. They again should be equivalent.
```
print "Function integral on polar grid:", c1.int_dA(abs(f(x1,y1)))
print "Function integral on cartesian grid:", c2.int_dA(abs(f(x2,y2)))
```

![http://polymode.googlecode.com/svn/wiki/figures/ex6_coordinates.png](http://polymode.googlecode.com/svn/wiki/figures/ex6_coordinates.png)


### Calculations with coordinates ###
This tutorial file is `ex6_calculating_with_coordinates.py`

Now most actions with modes can use coordinates. In particular we can plot what
the mode fields look like sampled on the points of the particular coordinate object
we have chosen simply by giving the `coord=` paramter to the plot command.
```
subplot(121)
m.plot(coord=c1)

subplot(122)
m.plot(coord=c2)
```

This plots the same mode on the polar and cartesian grids that we created in
the previous tutorial. All the same options for the `.plot()` method that
have been available can be used. Notice the small white "hole" in the center
of the lefthand plot. This is because we removed the origin of the polar coordinate
with the `border=1` parameter, try it again with `border=0`.
Try plotting the _x_ and _y_ components of the electric and magnetic fields
using both coordinates.

![http://polymode.googlecode.com/svn/wiki/figures/ex6_calculating_with_coordinates.png](http://polymode.googlecode.com/svn/wiki/figures/ex6_calculating_with_coordinates.png)

Also we can use the coord objects as options when we calculate properties of
the mode, such as the group index.
```
print "The group index using polar coordinates:", m.group_index(wg, coord=c1).real
print "The group index using cartesian coordinates:", m.group_index(wg, coord=c2).real
```

Now let's calculate something using a coordinate object. Here we calculate the
integral approximation to the propagation constant (Synder & Love, 1983, pg 222)

![http://mathurl.com/?img=ceem4w&x.png](http://mathurl.com/?img=ceem4w&x.png)
<a href='Hidden comment: 
\beta = Z_0 k \frac{ \int_A n2 \mathbf{e} \times h* \cdot \hat{z} dA}{ \int_A n2 |\mathbf{e}|2 dA}
'></a>

Where _Z₀_ is the impedance of free space. Note that througout the code the electric field is scaled by the impedence of free space, so we can take _Z₀_=1,

![http://mathurl.com/?img=bjhgob&.png](http://mathurl.com/?img=bjhgob&.png)
<a href='Hidden comment: 
E = Z_0 E_\textrm{calculated}
'></a>

```
e = m.electric_field(wg, coord=c1)
exhz = m.ExHc_t(coord=c1)
e2 = sum(abs(e)**2,axis=0)

n2 = m145.index(wl)**2

beta = m.k0*c1.int_dA(n2*exhz)/c1.int_dA(n2*e2)
```


### Retreiving mode fields ###
This tutorial file is `ex6_modal_fields.py`

Now we can extract the electric and magnetic fields from the modes we can do a lot of calculation with them.

The fields are extracted using the methods
  * _mode_.magnetic\_field(cartesian=True/False, coord=_coordinate object_)
  * _mode_.electric\_field(wg=_waveguide_, cartesian=True/False, coord=_coordinate object_)

These return the magnetic field vector and the electric field vector respectively. The magnetic field vector
will be h=(h\_x, h\_y, h\_z) if the coordinate object is cartesian or the cartesian option is given as true, and h=(h\_r, h\_ϕ, h\_z) if a polar coordinate object is given, or if the cartesian option is given as false. This is similar for the electric field.

Note that the waveguide can be given to the electric field to aid the calculation of the logitudinal component of the electric field, e\_z. If it is not given a constant refractive index is assumed throughout, which may or may not be a good approximation to the field depending on the location of the field intensity.

Similarly the transverse field components only can be extracted using the following matching methods

  * _mode_.magnetic\_transverse\_field(cartesian=True/False, coord=_coordinate object_)
  * _mode_.electric\_transverse\_field(cartesian=True/False, coord=_coordinate object_)

Following we use these commands to estimate the relative electric and magnetic fields in the _z_ direction,
relative to the transverse direction. We can use this to estimate if a mode is TM-like or TE-like.
```
hr,ha,hz = m.magnetic_field(coord=c1)
er,ea,ez = m.electric_field(wg, coord=c1)

print "TM like mode factor", sqrt(abs(hr)**2+abs(ha)**2).sum()/abs(hz).sum()
print "TE like mode factor", sqrt(abs(er)**2+abs(ea)**2).sum()/abs(ez).sum()
```

The same can be done for the cartesian coordinate, notice that in this case we are automatically returned the cartesian vector field.
```
h2 = m.magnetic_field(coord=c2)
e2 = m.electric_field(wg, coord=c2)
hx,hy,hz = h2
ex,ey,ez = e2
```

The Poynting vector is the cross product of the electric and conjugate magnetic fields, ![http://mathurl.com/?img=ao4vzb%.png](http://mathurl.com/?img=ao4vzb%.png)
<a href='Hidden comment: 
S = \frac{1}{2} e \times h^*
'></a>

```
poynting = cross(e2,conj(h2), axis=0)

#The x-component of the poynting vector
subplot(131)
Plotter.plot_v(c2, poynting[0].real)
wg.plot(fill=0, boundary=0)

subplot(132)
Plotter.plot_v(c2, poynting[1].real)
wg.plot(fill=0, boundary=0)

subplot(133)
Plotter.plot_v(c2, poynting[2].real)
wg.plot(fill=0, boundary=0)
```

We plot the real part of the three components of the Poynting vector as usual.

![http://polymode.googlecode.com/svn/wiki/figures/ex6_modal_fields.png](http://polymode.googlecode.com/svn/wiki/figures/ex6_modal_fields.png)


### Degenerate modes ###
This tutorial file is `ex6_degenerate_modes.py`

As was mentioned earlier modes in fibers with a symmetry 3 or greater in mode classes m₀ ≠ 0 are degenerate.
The mode calculated by **Polymode** will be elliptic and the two degenerate modes can be calculated from one elliptically
polarized mode.

Of course in fibers that have rotational symmetry that is not exploited by the code (the symmetry=1 option is used in the waveguide) will still have degenerate modes however together they will not be able to be generated from a single calculated mode.

To construct a pair of orthogonal degenerate modes we use the two functions

  * Modes.construct\_degenerate\_pair(_mode_)
  * Modes.construct\_lp\_degenerate\_pair(_mode_)

> Where `construct_degenerate_pair` constructs two circularly polarized modes and `Modes.construct_lp_degenerate_pair` constructs two linearly polarized modes.

```
mc1,mc2 = Modes.construct_degenerate_pair(modes[0])

m1,m2 = Modes.construct_lp_degenerate_pair(modes[0])

Plotter.plot_modes_in_grid([mc1,mc2,m1,m2], 'pe')
```

Run this example to plot the electric polarizations of the two mode pairs, and note that the circularly polarized modes are orthogonal - the arrows on the circles point in opposite directions.

![http://polymode.googlecode.com/svn/wiki/figures/ex6_degenerate_modes.png](http://polymode.googlecode.com/svn/wiki/figures/ex6_degenerate_modes.png)


<a href='Hidden comment: 



----

%== Tutorial 7 ==

%===  Compressing modes ===

Often when calculating large numbers of modes we don"t need to store the full mode
vectors. We can use the option to compress the mode vectors to a smaller size for
storage - then they will take up considerably less space when saved to a file.

The compressed modes can be plotted and used for calculation as normal, although
of course the accuracy of any computations will not be the same as when using the
full mode vectors.

To demonstrate compressing the modes we use the _pickle_ Python module directly.
In actual code you won"t need this - just use the load_data() and save_data()
functions of Polymode.
```
import pickle
```

First we create the regular solver and find 2 modes
```
solver = NLSolver.DefaultSolver(Nx, wg)
modes = solver(0, wl, 1.45, number=2)
```

Now to compress the modes calculated we use compress_to_size parameter
when creating the solver. It takes a _tuple_ containing two values, the number of
radial and azimuthal nodes to compress to as (N_r, N_ϕ).
```
solver = NLSolver.DefaultSolver(Nx, wg, compress_to_size = (10,5))
modes_comp = solver(0, wl, 1.45, number=2)
```

We can also choose to throw away all field information completely, this
will mean we can"t plot the mode at all but is useful if all that is needed is the
effective index of the mode.
```
solver = NLSolver.DefaultSolver(Nx, wg, discard_vectors = True)
modes_novec = solver(0, wl, 1.45, number=2)
```

Using the pickle module"s _dumps_ function we can estimate the file size of the
saved data for all three cases.
```
print "Size of uncompressed modes: %dkb" % (len(pickle.dumps(modes))/1024)
print "Size of compressed modes: %dkb" % (len(pickle.dumps(modes_comp))/1024)
print "Size of modes without vectors: %dkb" % (len(pickle.dumps(modes_novec))/1024)
```

We plot both modes. Notice there"s very little difference even with only 10x5 node points.
This is because the mode is very smooth - for modes with fine detail more will be lost.

http://polymode.googlecode.com/svn/wiki/figures/fig_compressmodes.png

%=== Batch solving over a range of neffs ===

%=== Batch solving over a range of wavelengths ===

%== Tutorial 8 ==

%=== Solving parallel jobs with MPI ===

'></a>
