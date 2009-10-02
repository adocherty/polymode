# _*_ coding=utf-8 _*_
from __future__ import division
from pylab import *
from numpy import *

from Polymode import *

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas

from PyQt4.QtCore import *
from PyQt4.QtGui import *

class Calculate(QThread):
    def __init__(self, parent = None):
        QThread.__init__(self, parent)
        self.exiting = False

        ## Solver parameters
        self.Nx=300,60
        self.m0 = 1

        ## Waveguide Parameters
        self.radius1 = 3.5
        self.radius2 = 3
        self.D = 6.8
        self.Nrings = 2

    def create_waveguide(self):
        ## Materials
        silica = Material.Silica()
        high = Material.SiO2GeO2(0.2)   #Silica with 20% Germania
        clad = Material.SiO2Fl(0.011)           #Silica with 1.1% Flourine

        ## Create waveguide
        self.wg = Waveguide.Waveguide(material=clad, symmetry=2)

        #Refractive index function
        fn_parabolic = lambda d: 1-(d/3)**2

        #Construct the waveguide
        core = Waveguide.Circle(silica, center=(0,0), radius=self.radius1)
        rings = []
        for i in range(self.Nrings):
            rings += [ Waveguide.Circle(silica, center=((i+1)*self.D,0), radius=self.radius1) ]
            inclusion = Waveguide.Circle(high, center=((i+1)*self.D,0), radius=self.radius2, zorder=1)
            inclusion.set_index_function(fn_parabolic, background=silica)
            rings += [ inclusion ]

            self.wg.add_shapes(core, rings)
        return self.wg

    def __del__(self):
        self.exiting = True
        self.wait()

    def run(self):
        wg = self.create_waveguide()

        #Create the solver
        solver = NLSolver.DefaultSolver(wg, self.Nx, compress_to_size=(50,10))

        #Solve at difference wavelengths
        wls=arange(1.,1.5,0.0025)

        self.allmodes=[]
        self.birefringence=[]
        neffapprox=wg.shapes[0].material.index(wls[0])-1e-3
        bifi=0
        for wl in wls:
            modes = solver(wl, self.m0, neffapprox, number=2)

            if self.exiting: break

            if len(modes)>1:
                #Find x polarized mode:
                pa = modes[0].polarization_angle()

                if abs(pa)<pi/4:        #x-polarized
                    bifi = modes[1].beta-modes[0].beta
                else:   #y-polarized
                    bifi = modes[0].beta-modes[1].beta
                neffapprox = modes[0].neff

            self.birefringence.append(bifi)
            self.allmodes+=[modes]
            self.emit(SIGNAL("output(float,float)"), wl, bifi.real)
        save_data((self.birefringence, self.allmodes), "ex7_bifi.dat")



class DisplayPlot(FigureCanvas):
    '''The Main QT application interface'''
    def __init__(self, parent=None):
        FigureCanvas.__init__(self, Figure())
        self.setParent(parent)

        self.ax = self.figure.add_subplot(111)
        #self.ax.hold(True)
        self.ax.grid()
        self.draw()

    def addpoint(self, x,y):
        print "Add point", x, y
        self.ax.scatter(x,y, marker='o')
        self.draw()


class ApplicationWindow(QMainWindow):
    about_message = "This program is a simple example of a Qt4 application embedding matplotlib canvases."
    def __init__(self):
        QMainWindow.__init__(self)
        self.setAttribute(Qt.WA_DeleteOnClose)
        self.setWindowTitle("Birefringence Example")
        self.setWindowIcon(QIcon('images/polymode_icon.png'))

        #Add menus
        self.file_menu = QMenu('&File', self)
        self.file_menu.addAction('&Quit', self.close,Qt.CTRL+Qt.Key_Q)
        self.menuBar().addMenu(self.file_menu)

        self.help_menu = QMenu('&Help', self)
        self.menuBar().addSeparator()
        self.menuBar().addMenu(self.help_menu)

        self.help_menu.addAction('&About', self.about)

        #container widget
        self.main_widget = QWidget(self)

        #Add plot window
        dp = DisplayPlot(self.main_widget)

        #Add buttons in horizontal layout
        quit = QPushButton('Close', self)
        start = QPushButton('Start', self)

        #Paramters
        quit = QLabel('Radius 1', self)
        start = QPushButton('Radius 2', self)
        quit = QLabel('Separation', self)
        start = QPushButton('Number', self)


        hbox = QHBoxLayout()
        hbox.addStretch(1)
        hbox.addWidget(quit)
        hbox.addWidget(start)

        grid = QGridLayout(self.main_widget)
        grid.addWidget(dp)

        vbox = QVBoxLayout(self.main_widget)
        vbox.addStretch(1)
        vbox.addWidget(dp)
        vbox.addLayout(hbox)
        grid.addLayout(hbox)

        self.setLayout(vbox)

        self.connect(quit, SIGNAL('clicked()'), self.close)

        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)

        self.statusBar().showMessage("All hail matplotlib!", 2000)

        #Setup calculation thread
        self.thread = Calculate()
        #self.connect(self.thread, SIGNAL("finished()"), self.updateUi)
        #self.connect(self.thread, SIGNAL("terminated()"), self.updateUi)
        #self.connect(self.thread, SIGNAL("output(float,float)"), self.addpoint)

        #Start calculations
        #self.thread.start()

    def closeEvent(self, event):
        "Called on window closure"
        event.accept()
        return
        reply = QMessageBox.question(self, 'Message', "Are you sure?", QMessageBox.Yes, QMessageBox.No)
        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

    def about(self):
        QMessageBox.about(self, about_message)

app = QApplication(sys.argv)

aw = ApplicationWindow()
aw.show()
sys.exit(app.exec_())
