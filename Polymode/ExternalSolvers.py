# _*_ coding=utf-8 _*_

"""
Solver class to interface with JCMsuite
"""
import os,sys,logging
import Waveguide, Modes

from numpy import *
from .Solver import Solve

class JCMWaveSolver(Solve):
    '''
    External solver class for the JCMSuite solvers
    supply the project_path where the files for the solver will be created
    and the main_path where the JCMsuite is installed
    '''
    jcmwave = None
    numthreads = 1
    shape_refinement = 3
    wg_domain_refinement = 3
    grid_refinements = 1
    fem_degree = 3

    path = "."
    results_path = "project_results"
    layout_filename = "layout.jcm"
    material_filename = "materials.jcm"
    triangulator_filename = "triangulator.jcm"
    project_filename = "project.jcmp"
    
    def __init__(self, wg, project_path=".", main_path=None, numthreads=1):
        self.path = project_path
        self.numthreads = numthreads
        
        #Check path
        if not (os.path.isdir(project_path) and os.access(project_path, os.W_OK)):
            logging.error("Project path does not exist!")
        
        #Append python path if specified
        if main_path is not None:
            sys.path.append(os.path.join(main_path,"ThirdPartySupport/Python"))

        #Try to load the JCM python interface
        try:
            self.jcmwave = __import__("jcmwave")
        except ImportError:
            import traceback
            result = traceback.format_exc()

            logging.error("Couldn't import jcmwave interface!")
            logging.info("%s" % result)

        #Run the rest of the configuration
        Solve.__init__(self,wg)

    def setup(self, geo_accuracy=1, refinement_levels=1, tol=1e-8, xy_domain=False):
        self.wg_domain_xy = xy_domain
        self.shape_refinement = geo_accuracy
        self.wg_domain_refinement = 1
        self.grid_refinements = refinement_levels

        #Create wg file
        self.create_materials_list()
        self.export_wg(domainxy=self.wg_domain_xy)

        #Create triangulation
        self.export_triangulation(view=0)
        
        #Clear mode list
        self.modes = []

    def create_materials_list(self):
        #Start material list:
        material_id = 1

        #Add waveguide material, and internal and external materials
        self.materials = {}
        self.materials[self.wg.material] = material_id

        #Also update zorder
        self.zorder_min = inf

        for shape in self.wg.shapes:
            if shape.material not in self.materials:
                material_id += 1
                self.materials[shape.material]= material_id

            self.zorder_min = min(shape.zorder, self.zorder_min)

        return self.materials

    def export_wg(self, domainxy=False, sectors=None):
        '''
        Export the internal PolyMode waveguide shapes to a JCMSuite file
        '''
        if sectors is None:
            sectors = self.wg.symmetry
        
        pml_sf = 2.0
        
        wglims = self.wg.extents(xy=domainxy)
        if domainxy:
            Cx = (wglims[1]+wglims[0])/2; Cy = (wglims[3]+wglims[2])/2
            Dx = wglims[1]-wglims[0]; Dy = wglims[3]-wglims[2]

            wg_domain_shape = "Parallelogram"
            wg_domain_segment_labels = ["South", "West", "East", "North"]
            wg_domain_size = "SideLengthA = %g\n SideLengthB = %g\n" % (Dx,Dy)
            wg_domain_center = "GlobalPositionX = %g\n GlobalPositionY = %g\n" % (Cx,Cy)
            pml_type = "NormalConstruction"
        else:
            wg_domain_shape = "Circle"
            wg_domain_segment_labels = ["South", "West", "East", "North"]
            wg_domain_size = "Radius = %g\n" % wglims[1]
            wg_domain_center = "GlobalPositionX = %g\n GlobalPositionY = %g\n" % (0,0)
            pml_type = "RadialConstruction"
        
        wg_material_id = self.materials[self.wg.material]
        quad_construction = \
 """   QuadConstruction {
     RadialScalingFactor = %g
     QuadConstructionType = %s
     MaterialId = %d
  }""" % (pml_sf, pml_type, wg_material_id)

        with open(os.path.join(self.path,self.layout_filename), "w") as jcmfile:
            jcmfile.write("Layout {\n")
            jcmfile.write(' Name = "Polymode Layout"\n')
            jcmfile.write(" UnitOfLength = %e\n" % 1e-6)
    
            #Computational domain
            jcmfile.write(" %s {\n" % wg_domain_shape)
            jcmfile.write(' Name = "Waveguide"\n')
            jcmfile.write(" %s" % wg_domain_size)
            jcmfile.write(" %s" % wg_domain_center)
            jcmfile.write(" MaterialId = %d\n" % wg_material_id)
            jcmfile.write(" Priority = -1\n")
            jcmfile.write(" RefineAll = %d\n" % self.wg_domain_refinement)

            for ii in range(4):
                jcmfile.write("  BoundarySegment {\n")
                jcmfile.write("   Segment = %s\n" % wg_domain_segment_labels[ii])
                jcmfile.write("   BoundaryClass = TransparentBoundary\n")
                jcmfile.write("   %s\n" % quad_construction)
                jcmfile.write("  }\n")
            jcmfile.write("  }\n")
    
            #Save object for each shape
            for ii,shape in enumerate(self.wg.shapes):
                for kk in range(sectors):
                    krotate = kk*2*pi/self.wg.symmetry

                    #If shape overlaps the origin then don't copy it .. JCMgeo can complain
                    rmin,rmax,_,_ = shape.extents()
                    if rmin==0 and kk>0:
                        continue
            
                    #Rotated center
                    r0,phi0 = shape.get_center()
                    centerxy = array([r0*cos(phi0+krotate), r0*sin(phi0+krotate)])

                    if isinstance(shape, Waveguide.Circle):
                        jcmfile.write(" Circle {\n")
                        jcmfile.write('   Name = "Circle %d.%d"\n' % (ii,kk))
                        jcmfile.write("   Radius = %g\n" % shape.radius)
                        jcmfile.write("   RefineAll = %s\n" % self.shape_refinement)
                
                    elif isinstance(shape, Waveguide.Rectangle):
                        jcmfile.write(" Parallelogram {\n")
                        jcmfile.write('   Name = "Rectangle %d.%d"\n' % (ii,kk))
                        jcmfile.write("   SideLengthA = %g\n   SideLengthB = %g\n" % shape.axes)
                        jcmfile.write("   RefineAll = %s\n" % self.shape_refinement)

                    elif isinstance(shape, Waveguide.Polygon):
                        jcmfile.write(" Polygon {\n")
                        jcmfile.write('   Name = "Polygon %d.%d"\n' % (ii,kk))
                        jcmfile.write("   Points = [\n")

                        xypoints = shape.to_nodelist(rotate=krotate)

                        #Ensure shapes are "open"
                        if xypoints[-1]==xypoints[0]:
                            xypoints=xypoints[:-1]

                        #Transform points - we should fix this to work in PolyMode itself
                        for xy in xypoints:
                            xy = array(xy)-centerxy
                            jcmfile.write("  %g %g\n" % tuple(xy))
                        jcmfile.write(" ]\n")

                    else:
                        print "Shape %s couldn't be converted" % shape
                        break

                    jcmfile.write("   GlobalPositionX = %g\n" % centerxy[0])
                    jcmfile.write("   GlobalPositionY = %g\n" % centerxy[1])
                    jcmfile.write("   MaterialId = %d\n" % self.materials[shape.material])
                    jcmfile.write("   Priority = %g\n" % (shape.zorder-self.zorder_min+1))
                    jcmfile.write(" }\n")   #Close shape

            jcmfile.write("}\n")    #Close Layout

    def export_materials(self, wl):
            #Write materials:
            with open(os.path.join(self.path,self.material_filename), "w") as jcmfile:
                for mat in self.materials:
                    #Calulate permittivity at wavelength
                    eps = mat.index(wl)**2
                    jcmfile.write("Material {\n")
                    jcmfile.write("  Id = %d\n" % self.materials[mat])
                    jcmfile.write('  Name = "%s"\n' % mat.__class__.__name__)
                    jcmfile.write("  RelPermeability = %g\n" % 1.0)
                    if abs(eps.imag)>abs(eps.real*1e-12):
                        jcmfile.write("  RelPermittivity = (%g,%g)\n" % (eps.real, eps.imag))
                    else:
                                            jcmfile.write("  RelPermittivity = %g\n" % eps.real)
                    jcmfile.write("}\n")    #Close material

    def export_triangulation(self, view=False):
        with open(os.path.join(self.path,self.triangulator_filename), "w") as jcmfile:

            jcmfile.write("""GeoControl {
  Output {
    TriangulationFormatJCM = yes
  }""")
    
            if view: jcmfile.write("""  Graphic {
    Triangulation = yes
    DependencyTree = no
    ExactGeometry = no
    Delay = -1
  }""")
            else: jcmfile.write("""  Graphic {
    Triangulation = no
    DependencyTree = no
    ExactGeometry = no
    Delay = 0
  }""")
    
            jcmfile.write("""  Triangulator {
    Triangulation = yes
    MinimumAngle = 20.0
  }
}""")

        #Run triangulator
        if self.jcmwave:
            self.jcmwave.set_num_threads(self.numthreads)
            self.jcmwave.geo(self.path)

    def external_solve(self, neffguess, number=1, run=True):
        keys = {}
        keys['wavelength'] = float(self.wl*1e-6)
        keys['n_eigenvalues'] = int(number)
        keys['eigenvalue_guess_re'] = float(real(neffguess))
        keys['eigenvalue_guess_im'] = float(imag(neffguess))
        keys['fem_degree'] = self.fem_degree
        keys['max_n_refinement_steps'] = int(self.grid_refinements)
        keys['leaky_mode'] = 1

        project_filename = os.path.join(self.path,self.project_filename)

        if self.jcmwave:
            self.jcmwave.runembeddedscripts(project_filename+'t',  keys)
            self.jcmwave.set_num_threads(self.numthreads)
            if run: self.jcmwave.solve('--solve', project_filename)

    def load_external_solutions(self, filename=None):
        if filename is None:
            filename = os.path.join(self.path, self.results_path, 'eigenvalues.jcm')
    
        jcmeigenvalues = self.jcmwave.loadtable(filename)
        solutions = []
        for column in jcmeigenvalues['columns']:
            if column['name'] == 'effective_refractive_index':
                solutions = column['data']

        #Construct modes
        modes = []
        k0 = 2*pi/self.wl
        for neff in solutions:
            ev = neff**2*k0**2
            modes.append(Modes.Mode(wl=self.wl, m0=self.m0, evalue=ev, symmetry=1))

        return modes
    
    def calculate(self, number=None):
        #Update materials
        self.export_materials(self.wl)

        finished = 0
        kk = 0
        while not finished and self.numbercalculated<self.totalnumber:
            #Choose next neffguess
            if self.nefflist is not None:
                neffguess = self.nefflist[kk]
            elif self.modelist is not None:
                neffguess = self.modelist[kk]
            else:
                neffguess = self.bracket[1]

            logging.info("Looking for modes near n_eff=%s" % neffguess)
            
            #Call JCMsolve to solve for the modes
            self.external_solve(neffguess, number=self.totalnumber, run=True)
            modes = self.load_external_solutions()
            
            self.numbercalculated += len(modes)
            self.modes += modes

            #Are we finished?
            if any([m.neff<self.bracket[0] for m in modes]):
                finished = 1
        
        #Sort modes in finalization method
        self.modes.sort()
        return self.modes



def export_to_nader(solver, prefix=""):
    "Export the waveguide and paramters file for Nader's solver"
    #Write waveguide to file
    wgfilename = prefix+"waveguide.in"
    paramfilename = prefix+"input_param.in"
    print "Saving waveguide to '%s'" % wgfilename

    wg = solver.wg
    wl = solver.wl
    Nx = solver.base_shape
    solver.equation.setup(Nx,wg,solver.m0,wl)
    wg.amerge.write_to_file(wl, wgfilename)

    rmax = wg.rmax
    ni = complex(wg.material.index(wl))
    next = complex(wg.exterior_index(wl))

    neffstart = sqrt(solver.evapprox)*wl/(2*pi)
    
    #Write input param to file
    print "Saving parameters to %s" % paramfilename
    spacer = "-"*80 + "\n"

    paramfile = open(paramfilename, "w")
    paramfile.write(spacer)
    paramfile.write("Vector or scalar (V/S)\t\tPlot mode fields (Y/N)\t\tTrack modes during wavelength scan (Y/N)\n")
    paramfile.write("V\t\tY\t\tY\n")
    paramfile.write(spacer)
    paramfile.write("R (Domain radius)\tReal[n_host]\tImag[n_host]\tReal[n_external]\tImag[n_external]\n")
    paramfile.write("%g\t\t%g\t\t%g\t\t%g\t\t%g\n" % (rmax, ni.real, ni.imag, next.real, next.imag))
    paramfile.write(spacer)
    paramfile.write("Number of different Resolutions\tWaveguide Symmetry\n")
    paramfile.write("%d\t\t%d\n" % (1, wg.symmetry))
    paramfile.write( spacer )
    paramfile.write("Radial resolution (n_max)\tAngular resolution (m_max)\n")
    paramfile.write("%d\t\t%d\n" % (Nx[0], (Nx[1]//2)*wg.symmetry))
    paramfile.write( spacer )
    paramfile.write("Starting wavelength\tEnding wavelength\tNumber of sampled wavelengths in between\n")
    paramfile.write("%g\t\t%g\t\t%d\n" % (wl,wl,1))
    paramfile.write( spacer )
    paramfile.write("Number of different search effective indices\n%d\n" % 1)
    paramfile.write( spacer )
    paramfile.write("Real[Search n_eff]\tImag[Search n_eff]\n")

    paramfile.write("%g\t\t%g\n" % (neffstart,0))

    paramfile.write( spacer )
    paramfile.write("Number of mode classes\n%d\n" % 1)
    paramfile.write( spacer )
    paramfile.write("Mode class (m0)\tNumber of modes in this class\n")
    paramfile.write("%d\t\t%d\n" % (solver.m0, solver.totalnumber))
    paramfile.write( spacer )
    paramfile.write("Number of ABC refinements\n3\n")
    paramfile.write( spacer )

