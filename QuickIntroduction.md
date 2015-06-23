# Introduction to Polymode #

Polymode is object oriented, open source and is written in Python with extentions in C++.
It aims to allow the specification of a fiber geometry with the minimum of commands in the most intuitive and flexible was

The computation implements the following features:
  * Exploitation of azimuthal symmetry
  * Efficient use of boundary conditions at the internal core interface and at the
> > external cladding interface to reduce the computational domain
  * Multiple numerical algorithms to solve different structures.

> (Currently only the finite difference method is fully implemented)

The main objects that the user sees are:
  * Waveguide: The MOF structure is created as a Waveguide object
> > containing
  * Material: Each material contains the dispersion profile of the
> > material over a specified wavelength range. Polymode contains an internal
> > database of materials and can import definitions from different formats.
  * Mode: Modes calculated are stored with the field information and
> > effective index allowing easy manipulation and storage.
  * Solver: Different numerical solvers can be used for each problem to
> > find the modes of the structure

In the following sections we take a quick tour of each of the modules of Polymode.

## Materials ##

Polymode contains an inbuilt set of commonly used materials.

```
air = Material.Air()
m1 = Material.Silica()
m2 = Material.PMMA()
```

Materials have interactive properties:
```
m3.index(1.55)
>> 1.4761908338569294

m3.info()
>> Silica doped with Germania
>> Wavelength limits (um): 0.6->1.8
>> Temperature limits (K): 0->400
>> Cite: Sunak, H.R.D.; Bastien, S.P., Photonics Technology Letters V1 N6 142-145, 1989
```

Most objects including materials can be visualised with a the plot method, in this case plotting
the material dispersion.

```
m2.plot()
m3.plot()
```

![http://polymode.googlecode.com/svn/wiki/figures/qi_material_plot.png](http://polymode.googlecode.com/svn/wiki/figures/qi_material_plot.png)

More materials can be loaded from [Freesnell](http://swiss.csail.mit.edu/~jaffer/FreeSnell/) [Luxpop](http://www.luxpop.com/) or [Sopra](http://www.sopra-sa.com/) formatted files. It also contains data from the freely available Freesnell materials.

## Waveguide Structure ##
Once the materials have been defined, we can use them to form an arbitrary structure.
We can place different shapes including circules, elipses, spline curves and polygons.

The waveguide is created defining a subtrate material (which also extends past the computational boundary) and a symmetry.
```
wg = Waveguide.Waveguide(material=m2, symmetry=6)
```

For example to create a six hole waveguide we place a single hole, as we have specified the rotational symmetry of the waveguide to be 6 the hole implicitly occurs 6 times.
```
s = Waveguide.Circle(air, center=(6.75,0), radius=2.0)
wg.add_shape(s)
```

We can see this plotting the waveguide:
```
wg.plot()
```

![http://polymode.googlecode.com/svn/wiki/figures/qi_hex_lattice.png](http://polymode.googlecode.com/svn/wiki/figures/qi_hex_lattice.png)

## Modes ##
Once a structure is solved it will return a list of modes. Each mode object contains the effective index,
field information and other properties of that mode. The info method can be used to see some important
information abou the mode.

```
mode.info()
	>> Mode, size: (200, 30), symmetry: C6, m0: 0
	>> wl: 1, r: 0 -> 14.7
	>> neff=1.47711 + 2.38e-10i, loss=0.013dB/m, 
	>>  | Effective Area: 32.441 um
	>>  | Spot size: 1.2365 um
	>>  | Numerical Aperture: 0.098573
```

And the plot method will plot the mode power distribution; the default is to plot the z component of the
time averages Poynting vector, but any other component of the vector E and H fields can be visualized.

```
mode.plot()
```

![http://polymode.googlecode.com/svn/wiki/figures/qi_hex_modes.png](http://polymode.googlecode.com/svn/wiki/figures/qi_hex_modes.png)

## Example applications of Polymode ##

> Three quite different structures will be used to illustrate the program:
  1. A large solid core structure surrounded by a microstructure of holes connected by thin bridges
  1. A thin bridge kagome fibre
  1. A highly birefringent fibre
  1. A graded index microstructured polymer optical fibre (designed to have high bandwidth and low loss)

### Air-core kagome fiber ###

Complex strutures can be built from simple shapes in a variety of ways. There are also several helper functions that enable complicated patterns to be constructed quickly.

To create a Kagome fiber we create a central hexagonal shape and "coat" it with polymer.
This has the effect of creating six thin struts of thickness 

&lt;width&gt;

.

```
	p = Waveguide.RegularPolygon(air, Nsides=6, length=s1)
	cp = Waveguide.create_coated_shape(p, m1, width)
```

In the same way we create dodecagons for the first layer of fiber.
```
	p1 = Waveguide.RegularPolygon(air, Nsides=12, length=s2)
	cp1 = Waveguide.create_coated_shape(p1, m1, width)
```

We can add another layer in the same way, positioning them using the `Waveguide.create_hexagonal_tiling`
function.

![http://polymode.googlecode.com/svn/wiki/figures/qi_kagome_construction.png](http://polymode.googlecode.com/svn/wiki/figures/qi_kagome_construction.png)

After solving for 3 modes, we can plot them using the `Plotter.plot_modes_in_grid` function:

![http://polymode.googlecode.com/svn/wiki/figures/qi_kagome_modes.png](http://polymode.googlecode.com/svn/wiki/figures/qi_kagome_modes.png)

### A microstructured birefringent fiber ###

Silica based fibers doped with Flourine and Germania
```
	silica = Material.Silica()
	high = Material.SiO2GeO2(0.2)
	clad = Material.SiO2Fl(0.011)
```

We can create functional refractive index profiles, by specifying the relative index depending on the distance from the center of the object. This uses Python's `lambda` statement.
```
	fn_parabolic = lambda d: 1-(d/3)**2
```

Now create the shapes, a pure silica core and doped silica rods placed in a line
```
	core = Waveguide.Circle(silica, center=(0,0), radius=r1)
	wg.add_shape(inclusion)
	for i in range(Nrings):
	  inclusion = Waveguide.Circle(high, center=((i+1)*D,0), radius=r2, zorder=1)
	  inclusion.set_index_function(fn_parabolic, background=silica)
	  wg.add_shape(inclusion)
```

The resulting refractive index profile can be plotted:
![http://polymode.googlecode.com/svn/wiki/figures/qi_bifi_wg.png](http://polymode.googlecode.com/svn/wiki/figures/qi_bifi_wg.png)

And the two linearly polarized modes are found:
![http://polymode.googlecode.com/svn/wiki/figures/qi_bifi_modes.png](http://polymode.googlecode.com/svn/wiki/figures/qi_bifi_modes.png)

And finally we can plot the birefringence over the wavelength
\![http://polymode.googlecode.com/svn/wiki/figures/qi_bifi_wl.png](http://polymode.googlecode.com/svn/wiki/figures/qi_bifi_wl.png)

### GIMPOF Fiber ###

All of Polymode is object orientation, so custom classes can easily be created
```
class GIMPPOF(Waveguide):
  def add_holes(self, circles):
    air = Materials.Air()
    for x,y,r in circles:
      self.add_shape(Circle(air,center=(x,y),radius=r))
```

Use this custom waveguide
```
gimp90 = GIMPPOF(material=pmma,symmetry=9)
gimp90.add_holes([(19.158,6.973,1.261), ...])
```

![http://polymode.googlecode.com/svn/wiki/figures/qi_gimp_wg.png](http://polymode.googlecode.com/svn/wiki/figures/qi_gimp_wg.png)

The results of calculations in the structure for many modes of the highly multimodes
fiber. The group index is calcuated from the mode field.

![http://polymode.googlecode.com/svn/wiki/figures/qi_gimp_logloss.png](http://polymode.googlecode.com/svn/wiki/figures/qi_gimp_logloss.png)

This is just a brief introduction to what Polymode can do. For more detailed information see the [PolymodeTutorial](PolymodeTutorial.md).
