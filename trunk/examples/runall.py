#!/bin/env python
import os, sys
import matplotlib
matplotlib.use('agg')

from pylab import *

output_path='output'

def runtut(tut):
	#Test for existance
	figfilename = os.path.join(output_path,tut+'.png')
	if os.access(figfilename, os.F_OK):
		print "Skipping tutorial '%s'" % tut
	else:
		clf()
		print "\nRunning tutorial '%s'" % tut
		execfile(tut+'.py')
		savefig(figfilename)

runtut('ex1_materials')
runtut('ex1_the_waveguide')

runtut('ex2_step_index_modes')
runtut('ex2_leaky_modes')
runtut('ex2_plotting_modes')

runtut('ex3_microstructured_fiber')
runtut('ex3_finding_modes')
runtut('ex3_plotting_modes')
runtut('ex3_square_lattice')
runtut('ex3_more_shapes')

runtut('ex4_coated_holes')
runtut('ex4_merging_shapes')
runtut('ex4_modes_with_a_image_waveguide')
runtut('ex4_parabolic_index')
runtut('ex4_waveguide_from_image')
runtut('ex4_extra_materials')

runtut('ex5_modal_quantities')
runtut('ex5_mode_dispersion')
runtut('ex5_convergence_test')
runtut('ex5_uranus_fiber')
runtut('ex5_wavelength_scan')

runtut('ex6_coordinates')
runtut('ex6_calculating_with_coordinates')
runtut('ex6_modal_fields')
runtut('ex6_degenerate_modes')

runtut('ex7_compressing_modes')
runtut('ex7_form_birefringence')
runtut('ex7_birefringent_fiber')
runtut('ex7_birefringent_fiber_2')

runtut('ex8_kagome_fiber')

