# _*_ coding=utf-8 _*_
#
#---------------------------------------------------------------------------------
#Copyright © 2009 Andrew Docherty
#
#This program is part of Polymode.
#Polymode is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
#---------------------------------------------------------------------------------
"""
Functions for easy plotting of modes and modal properties

This is a light wrapper to matplotlib: http://matplotlib.sourceforge.net
"""

import logging

try:
    import pylab
    from pylab import *
except:
    logging.error("Plot library unavailable, plots cannot be made")
    pylab = None

from numpy import (pi,append,c_,cos,sin,newaxis,log,absolute,arange,exp,fft,real,imag,iterable,atleast_1d)


__spurious_warning_color__ = 'crimson'

def plot_v(coord, f, style='pcolor', cmap=None, color=None, alpha=1.0, aspect=None, plotstyle={}):
    """
    Plot the 2-d data f using a coordinate object.
        
    Parameters:
     coord: The coordinate object associated to the data f
     f: The 2-d array of data to be plotted
     style: plotting style, one of 'pcolor'*, 'contour', 'line', 'vector', 'circ'
     cmap: color map from pylab.cm
     color: line color for the line-based plotters
     alpha: the alpha transparancy of the plot
     plotstyle: specific matplotlib styles
     aspect: the plot aspect ratio
    """
    if cmap==None: cmap = pylab.cm.jet

    #Calculate the bases
    rm,phim = coord.polar2d(interval=0)
    xm,ym = coord.cartesian2d(interval=0)

    irm,iphim = coord.polar2d(interval=1)
    ixm,iym = coord.cartesian2d(interval=1)

    #Close plot if spans 2*pi and plotting with contour plot
    if (style.startswith('cont') or style.startswith('line')) and hasattr(coord,'arange'):
        if abs(abs(diff(coord.arange))-2*pi)<1e-6:
            xm = append(xm, xm[:,:1], axis=1)
            ym = append(ym, ym[:,:1], axis=1)
            f = append(f, f[:,:1], axis=1)

    #The actual plotting commands       
    autoaspect = 'equal'
    if style.startswith('cont'):
        V=10
        c=pylab.contourf(xm,ym,f.real,V,colors=color,cmap=cmap,alpha=alpha,linestyles=None)
        ax = gca()
        
    elif style.startswith('line'):
        c=pylab.contour(xm,ym,f.real,colors=color,cmap=cmap,alpha=alpha)
        ax = gca()

    elif style.startswith('pcol'):
        c=pylab.pcolor(ixm,iym,f.real,cmap=cmap,shading='flat',alpha=alpha)
        ax = gca()

    #Vector plot
    elif style.startswith('vector'):
        fx,fy = 0.1*f.real/abs(f).max()
        c=pylab.quiver(xm,ym,fx,fy, pivot='middle', scale=2, color=color)
        ax = gca()

    #Plot circular polarization
    elif style.startswith('circ'):
        from matplotlib import patches
        ax=gca()

        #Could make these a little more configurable!
        size = 0.4*min(coord.characteristic_length)
        dp=0.02
        arrowp = 0.2
        width=0.5
        phis = arange(dp,1+dp,dp)*2*pi

        fx,fy = size*f/abs(f).max()
        for ii in ndindex(fx.shape):
            cx,cy = fx[ii], fy[ii]
            xy = real([xm[ii] + cx*exp(1j*phis), ym[ii] + cy*exp(1j*phis)])
            e = patches.Polygon(xy.T, fill=0, ec=color, **plotstyle)
            ax.add_artist(e)
            
            #Only add arrow if ellipse is large enough
            if linalg.norm([cx,cy])>0.5*size:
                dx,dy = real([cx*exp(1j*arrowp)-cx, cy*exp(1j*arrowp)-cy])
                xyt = array([xm[ii]+real([cx-width*dy,cx+width*dy,cx+dx]),\
                            ym[ii]+real([cy+width*dx,cy-width*dx,cy+dy])])
                arrow = patches.Polygon(xyt.T, fill=1, ec=color, fc=color, **plotstyle)
                ax.add_artist(arrow)

        ax.axis([xm.min(), xm.max(), ym.min(), ym.max()])
                    
    elif style.startswith('3d'):
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt

        ax = Axes3D(gcf())
        ax.plot_surface(xm,ym,f.real, rstride=1, cstride=1, cmap=cmap)

    elif style.startswith('polar'):
        pylab.pcolor(irm,iphim,f.real,cmap=cmap,shading='flat',alpha=alpha)
        ax = gca()
        autoaspect = 'auto'
    
    else:
        raise NotImplementedError, "Plot type not implemented"
    
    #Set aspect ratio for plot
    ax.set_aspect(autoaspect if aspect is None else aspect)
    draw()

def plot_modes_in_grid(modes, plottype='Sz', wg=None, title=None, bar=False,
                    Nrows=None, axis_off=None, **plotargs):
    """
    Plot all given modes on one figure.
    
    Parameters:
     - plottype: Standard mode.plot() type
     - wg: if given the waveguide will be plotted
     - title: how to title the plots, either 'full', 'neff', 'number' or ''
     - bar: plot a colorbar on each plot
     - axis_off: hide the axes
     - Nx: resolution of plot
     - cartesian: True/False to plot cartesian or polar
    """

    #fig = pylab.figure()
    fig = gcf()
    plottype = plottype.replace(',',' ').split()
    modes = np.atleast_1d(modes)
    
    nmodes = len(modes)
    ntypes = len(plottype)

    #Set appropriate title
    if title is None and nmodes*ntypes<6: title='full'
    elif title is None: title=''
    
    #The dimensions of the grid
    if Nrows is None:
        Nrows = floor(sqrt(nmodes*ntypes))
    Ncols = ceil(nmodes*ntypes/Nrows)
    
    axis_ticks = nmodes<4
    ii = 0
    for ii in range(nmodes):
        for jj in range(ntypes):
            ax = fig.add_subplot(Nrows,Ncols,ii*ntypes+jj+1)
            
            #Only put ticks if there aren't too many modes
            if not axis_ticks:
                ax.set_xticks([])
                ax.set_yticks([])
            
            if axis_off:
                ax.set_axis_off()
            
            #if wg: wg.plot(fill=False)

            if title == 'number': title_str = "(%d)" % ii
            elif title == 'full': title_str = r"%(type)s, $n_{\mathrm{eff}}=%(tneff)s$"
            elif title == 'neff': title_str = r"$%(tneff)s$"
            else: title_str = title

            modes[ii].plot(plottype[jj], wg=wg, title=title_str, **plotargs)
            
            if wg:
                wg.plot(fill=0)
            
            #Colorize background if a spurious mode is detected
            if modes[ii].is_spurious:
                ax.set_axis_bgcolor(__spurious_warning_color__)
            
            if bar: np.colorbar()

    return fig

def extract_data_from_modes(modes=[], datatype='', return_label=False):
    """
    Return a list containing the requested modal paramters for each mode given
        modes: list of modes
        datatype: 'neff', 'loss', 'wl', 'ineff', 'Na', 'Nr'
        return_label: if true also return a latex formatted data label
    """
    #Deal with groups of modes
    modes = flatten(modes)

    if datatype.startswith('loss'):
        y = [ m.loss for m in modes ]
        lab = r'Loss, db/m'
    elif datatype.startswith('neff'):
        y = [ real(m.neff) for m in modes ]
        lab = r'$Re(n_{\rm eff})$'
    elif datatype.startswith('ineff'):
        y = [ imag(m.neff) for m in modes ]
        lab = r'$Im(n_{\rm eff})$'
    elif datatype.startswith('disp'):
        y = [ m.dispersion for m in modes ]
        lab = r'Dispersion'
    elif datatype.startswith('w'):
        y = [ m.wl for m in modes ]
        lab = r'Wavelength, $\mu$m'
    elif datatype.startswith('nr'):
        y = [ m.coord.Nr for m in modes ]
        lab = r'Radial resolution, $N_r$'
    elif datatype.startswith('na'):
        y = [ m.coord.Naz for m in modes ]
        lab = r'Azimuthal resolution, $N_\phi$'
    elif datatype in modes[0].label:
        y = [ float(m.label[datatype]) for m in modes ]
        lab = "%s" % datatype

    if return_label:
        return y,lab
    else:
        return y
    
def plot_mode_properties(modes=[], ydata='loss', xdata='wl', style='', sort2d=False):
    """
    Plot a graph of the specified modal properties
    modes: list of modes to extract property
    ydata: name of property on y axis
    xdata: name of property on x axis
    style: matplotlib linestyle to plot with,
            see help in pylab.plot
    
    xdata, ydata can be one of:
    'neff': the real part of the mode effective index
    'ineff': the imaginary part of the mode effective index
    'loss': the loss in db/km for the mode
    'dispersion': the estimated dispersion of the mode
    'wavelength': the mode wavelength
    'nr': the radial resolution of the calculation
    'naz': the azimuthal resolution of the calculation
    """
    ax = gca()
    
    x,xlab = extract_data_from_modes(modes, xdata.lower(),True)
    y,ylab = extract_data_from_modes(modes, ydata.lower(),True)

    if sort2d:
        xall = sort(unique(x))
        x2d = []; y2d = []
        done=False; jj=0
        while not done:
            done=True
            x_mm = []; y_mm = []
            for ii in range(len(xall)):
                yii = find(x==xall[ii])
                ys = sorted(array(y)[yii], reverse=True)
                if jj<len(ys):
                    done=False
                    x_mm += [ xall[ii] ]
                    y_mm += [ ys[jj] ]

            if not done:
                plot(x_mm, y_mm)
                x2d+=[x_mm]; y2d+=[y_mm]
                jj+=1
    else:
        ax.plot(x,y,style)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)

    return ax
    
