"""

:mod:`Inelastica.io.xmgrace`
============================

This module provides a `Python`_ interface to write `XMGR/GRACE <xmgrace_>`_ files.

It is a renamed version of the code `WriteXMGR.py` written by Thomas Frederiksen, July 2007.

Scripting
---------

To create a `GRACE <xmgrace_>`_ plot from a given data set, a minimal
script structure takes the following form:

.. code-block:: bash

    >>> from Inelastica.io.xmgrace import *
    >>> x,y = [1, 2], [3, 4]
    >>> data = XYset(x,y)
    >>> graph = Graph(data)
    >>> plot = Plot('test.agr',graph)
    >>> plot.WriteFile()

More graphs and data sets are easily appended. Continuing on the
example above:

.. code-block:: bash

    >>> x2,y2 = [0, 1], [2, 1]
    >>> data2 = XYset(x2,y2)
    >>> graph2 = Graph(data2)
    >>> plot.AddGraphs(graph2)
    >>> plot.ArrangeGraphs(nx=2,ny=1)
    >>> plot.WriteFile('test2.agr')

With these additional lines we end up with a plot containing two graphs.
The function call to ``ArrangeGraphs(...)`` arranges them nicely on the
canvas in a 2x1 layout, i.e., side by side.

A more detailed demonstration of this module can be generated with

.. code-block:: bash

    >>> import Inelastica.io.xmgrace as xmgr
    >>> xmgr.demo()

creating the files `demo.agr` and `demo.eps`.

.. currentmodule:: Inelastica.io.xmgrace

Classes
-------

.. autosummary::
   :toctree:

   Dataset
   XYset
   XYDXset
   XYDYset
   XYDXDYset
   XYSIZEset
   Graph
   Plot

"""
from __future__ import print_function

import time
import math
import os
import numpy as N

symbols = {'\alpha': '\\f{Symbol}a',
           '\beta': '\\f{Symbol}b',
           '\gamma': '\\f{Symbol}g',
           '\delta': '\\f{Symbol}d',
           '\lambda': '\\f{Symbol}l',
           '\pi': '\\f{Symbol}p',
           '\sigma': '\\f{Symbol}s',
           '\Gamma': '\\f{Symbol}G',
           '\Delta': '\\f{Symbol}D'
          }
flag = {True: 'on', False: 'off'}


def Array2XYsets(A, **keywords):
    """
    Converts a matrix  A = [X, Y1, Y2, ...] into a list of XYset objects.
    Formatting keywords (Lcolor, Stype, etc.) are applied to all the generated XYsets.
    """
    sets = []
    A = N.array(A)
    x = A[:, 0]
    for i in range(1, len(A[0])):
        y = A[:, i]
        sets.append(XYset(x, y, **keywords))
    return sets


def Datafile2XYsets(fn, Sort=False, **keywords):
    """
    Reads an ascii file and applies ``Array2XYsets(...)``
    """
    print('io.xmgrace.Datafile2XYsets: reading', fn)
    f = open(fn, 'r')
    A = []
    for line in f:
        if line[0] != '#' and len(line) > 0:
            l = line.split()
            for i, lval in enumerate(l):
                l[i] = float(lval)
            if len(l) >= 2:
                A.append(l)
    f.close()
    A = N.array(A)
    if Sort:
        print(' ... sorting along x-values')
        tmp = N.transpose(A)
        argi = tmp.argsort()[0] # sort x-values
        A = A[argi]
    print(' ... data shape =', N.shape(A))
    return Array2XYsets(A, **keywords)


class Dataset(object):

    """
    Generic class for any XMGR dataset
    """

    def __init__(self, **keywords):
        "Returning instance of class."
        self.legend = ''
        self.Ltype = ''
        self.Lstyle = ''
        self.Lwidth = ''
        self.Lcolor = ''
        self.Stype = ''
        self.Ssize = ''
        self.Scolor = ''
        self.Sfillcolor = ''
        self.Sfillpattern = ''
        self.Slinewidth = ''
        self.Slinestyle = ''
        self.EBsize = ''
        self.EBcolor = ''
        self.EBlinewidth = ''
        self.SetFormat(**keywords)

    def SetFormat(self, legend='', Ltype='', Lstyle='', Lwidth='', Lcolor='',
                  Stype='', Ssize='', Scolor='', Sfillcolor='', Sfillpattern='',
                  Slinewidth='', Slinestyle='',
                  EBsize='', EBcolor='', EBlinewidth=''):
        "Sets the formatting options (linestyle, symbol type, etc.)."
        if legend != '': self.SetLegend(legend)
        if Ltype != '': self.SetLtype(Ltype)
        if Lstyle != '': self.SetLstyle(Lstyle)
        if Lwidth != '': self.SetLwidth(Lwidth)
        if Lcolor != '': self.SetLcolor(Lcolor)
        if Stype != '': self.SetStype(Stype)
        if Ssize != '': self.SetSsize(Ssize)
        if Scolor != '': self.SetScolor(Scolor)
        if Sfillcolor != '': self.SetSfillcolor(Sfillcolor)
        if Sfillpattern != '': self.SetSfillpattern(Sfillpattern)
        if Slinewidth != '': self.SetSlinewidth(Slinewidth)
        if Slinestyle != '': self.SetSlinestyle(Slinestyle)
        if EBsize != '': self.SetEBsize(EBsize)
        if EBcolor != '': self.SetEBcolor(EBcolor)
        if EBlinewidth != '': self.SetEBlinewidth(EBlinewidth)

    def SetLegend(self, legend):
        "Sets the legend (string)."
        self.legend = legend

    def SetLtype(self, Ltype):
        "Sets the line type (integer)."
        self.Ltype = Ltype

    def SetLstyle(self, Lstyle):
        "Sets the line style (integer)."
        self.Lstyle = Lstyle

    def SetLwidth(self, Lwidth):
        "Sets the line width (float)."
        self.Lwidth = Lwidth

    def SetLcolor(self, Lcolor):
        "Sets the line color (integer)."
        self.Lcolor = Lcolor

    def SetStype(self, Stype):
        "Sets the symbol type (integer)."
        self.Stype = Stype

    def SetSsize(self, Ssize):
        "Sets the symbol size (float)."
        self.Ssize = Ssize

    def SetScolor(self, Scolor):
        "Sets the symbol color (integer)."
        self.Scolor = Scolor

    def SetSfillcolor(self, Sfillcolor):
        "Sets the symbol fill color (integer)."
        self.Sfillcolor = Sfillcolor

    def SetSfillpattern(self, Sfillpattern):
        "Sets the symbol fill pattern (integer)."
        self.Sfillpattern = Sfillpattern

    def SetSlinewidth(self, Slinewidth):
        "Sets the symbol linewidth (float)."
        self.Slinewidth = Slinewidth

    def SetSlinestyle(self, Slinestyle):
        "Sets the symbol line style (integer)."
        self.Slinestyle = Slinestyle

    def SetEBsize(self, EBsize):
        "Sets the error bar size (float)."
        self.EBsize = EBsize

    def SetEBcolor(self, EBcolor):
        "Sets the error bar color (integer)."
        self.EBcolor = EBcolor

    def SetEBlinewidth(self, EBlinewidth):
        "Sets the error bar linewidth (float)."
        self.EBlinewidth = EBlinewidth

    def __GetXMGRstring__(self, graphnr, setnr):
        string = ''
        if self.legend != '':
            string += '@ s%i legend \"%s\" \n'%(setnr, self.legend)
        if self.Ltype != '':
            string += '@ s%i line type %i \n'%(setnr, self.Ltype)
        if self.Lstyle != '':
            string += '@ s%i line linestyle %i \n'%(setnr, self.Lstyle)
        if self.Lwidth != '':
            string += '@ s%i line linewidth %f \n'%(setnr, self.Lwidth)
        if self.Lcolor != '':
            string += '@ s%i line color %i \n'%(setnr, self.Lcolor)
        if self.Stype != '':
            string += '@ s%i symbol %i \n'%(setnr, self.Stype)
        if self.Ssize != '':
            string += '@ s%i symbol size %f \n'%(setnr, self.Ssize)
        if self.Scolor != '':
            string += '@ s%i symbol color %i \n'%(setnr, self.Scolor)
        if self.Sfillcolor != '':
            string += '@ s%i symbol fill color %i \n'%(setnr, self.Sfillcolor)
        if self.Sfillpattern != '':
            string += '@ s%i symbol fill pattern %i \n'%(setnr, self.Sfillpattern)
        if self.Slinewidth != '':
            string += '@ s%i symbol linewidth %f \n'%(setnr, self.Slinewidth)
        if self.Slinestyle != '':
            string += '@ s%i symbol linestyle %i \n'%(setnr, self.Slinestyle)
        if self.EBsize != '':
            string += '@ s%i errorbar size %f \n'%(setnr, self.EBsize)
        if self.EBcolor != '':
            string += '@ s%i errorbar color %i \n'%(setnr, self.EBcolor)
        if self.EBlinewidth != '':
            string += '@ s%i errorbar linewidth %f \n'%(setnr, self.EBlinewidth)
            string += '@ s%i errorbar riser linewidth %f \n'%(setnr, self.EBlinewidth)
        return string


class XYset(Dataset):

    """
    Class for 2-dimensional (X,Y) datasets
    """

    def __init__(self, x, y, **keywords):
        "Returning instance of class."
        Dataset.__init__(self, **keywords)
        self.x = N.array(x)
        self.y = N.array(y)
        if len(x) != len(y):
            print('ERROR: x and y data lists are not of equal length!!!')
            error
        self.xmin = float(min(self.x))
        self.ymin = float(min(self.y))
        self.xmax = float(max(self.x))
        self.ymax = float(max(self.y))

    def GetWorld(self):
        "Returns tuple (xmin,ymin,xmax,ymax) of the data set."
        return self.xmin, self.ymin, self.xmax, self.ymax

    def GetXMGRstring(self, graphnr, setnr):
        "Returns a string containing the GRACE commands that describe the data set."
        string = '# ------------------ Dataset %i.%i ----------------- \n'%(graphnr, setnr)
        string += '@ target g%i.s%i \n@ type xy \n'%(graphnr, setnr)
        string += self.__GetXMGRstring__(graphnr, setnr)
        for i in range(len(self.x)):
            string += '   %.8f  %.8f \n'%(self.x[i], self.y[i])
        return string


class XYDXset(XYset):

    """
    Class for 3-dimensional (X,Y,dX) datasets
    """

    def __init__(self, x, y, dx, **keywords):
        "Returning instance of class."
        XYset.__init__(self, x, y, **keywords)
        self.dx = N.array(dx)

    def GetXMGRstring(self, graphnr, setnr):
        "Returns a string containing the GRACE commands that describe the data set."
        string = '# ------------------ Dataset %i.%i ----------------- \n'%(graphnr, setnr)
        string += '@ target g%i.s%i \n@ type xydx \n'%(graphnr, setnr)
        string += self.__GetXMGRstring__(graphnr, setnr)
        for i in range(len(self.x)):
            string += '   %.8f  %.8f  %.8f \n'%(self.x[i], self.y[i], self.dx[i])
        return string


class XYDYset(XYset):

    """
    Class for 3-dimensional (X,Y,dY) datasets
    """

    def __init__(self, x, y, dy, **keywords):
        "Returning instance of class."
        XYset.__init__(self, x, y, **keywords)
        self.dy = N.array(dy)

    def GetXMGRstring(self, graphnr, setnr):
        "Returns a string containing the GRACE commands that describe the data set."
        string = '# ------------------ Dataset %i.%i ----------------- \n'%(graphnr, setnr)
        string += '@ target g%i.s%i \n@ type xydy \n'%(graphnr, setnr)
        string += self.__GetXMGRstring__(graphnr, setnr)
        for i in range(len(self.x)):
            string += '   %.8f  %.8f  %.8f \n'%(self.x[i], self.y[i], self.dy[i])
        return string


class XYDXDYset(XYset):

    """
    Class for 4-dimensional (X,Y,dX,dY) datasets
    """

    def __init__(self, x, y, dx, dy, **keywords):
        "Returning instance of class."
        XYset.__init__(self, x, y, **keywords)
        self.dx = N.array(dx)
        self.dy = N.array(dy)

    def GetXMGRstring(self, graphnr, setnr):
        "Returns a string containing the GRACE commands that describe the data set."
        string = '# ------------------ Dataset %i.%i ----------------- \n'%(graphnr, setnr)
        string += '@ target g%i.s%i \n@ type xydxdy \n'%(graphnr, setnr)
        string += self.__GetXMGRstring__(graphnr, setnr)
        for i in range(len(self.x)):
            string += '   %.8f  %.8f  %.8f  %.8f \n'%(self.x[i], self.y[i], self.dx[i], self.dy[i])
        return string


class XYSIZEset(XYset):

    """
    Class for 3-dimensional (X,Y,size) datasets
    """

    def __init__(self, x, y, size, **keywords):
        "Returning instance of class."
        XYset.__init__(self, x, y, **keywords)
        self.size = N.array(size)

    def GetXMGRstring(self, graphnr, setnr):
        "Returns a string containing the GRACE commands that describe the data set."
        string = '# ------------------ Dataset %i.%i ----------------- \n'%(graphnr, setnr)
        string += '@ target g%i.s%i \n@ type xysize \n'%(graphnr, setnr)
        string += self.__GetXMGRstring__(graphnr, setnr)
        for i in range(len(self.x)):
            string += '   %.8f  %.8f  %.8f \n'%(self.x[i], self.y[i], self.size[i])
        return string


class Graph(object):

    """
    Class for a graph (containing datasets)
    """

    def __init__(self, *datasets):
        "Returning instance of class."
        self.datasets = []
        self.world = [0., 0., 1., 1.]
        self.view = [0.15, 0.15, 1.15, 0.85]
        self.string = ''
        self.legend = False
        self.SetWorld()
        self.SetView()
        if len(datasets) > 0:
            self.AddDatasets(*datasets)

    def SetTitle(self, title, font=0, size=1.5, color=1):
        "Sets the title of the graph."
        self.string += '@ title \"%s\"\n'%title
        self.string += '@ title font %i \n'%font
        self.string += '@ title size %.8f \n'%size
        self.string += '@ title color %i \n'%color

    def SetSubtitle(self, subtitle, font=0, size=1.0, color=1):
        "Sets the subtitle of the graph."
        self.string += '@ subtitle \"%s\"\n'%subtitle
        self.string += '@ subtitle font %i \n'%font
        self.string += '@ subtitle size %.8f \n'%size
        self.string += '@ subtitle color %i \n'%color

    def SetWorld(self, xmin=False, xmax=False, ymin=False, ymax=False):
        "Sets the graph world (defining the axis limits)."
        if xmin: self.world[0] = xmin
        if xmax: self.world[2] = xmax
        if ymin: self.world[1] = ymin
        if ymax: self.world[3] = ymax

    def GetWorldOfDatasets(self):
        "Returns tuple (xmin,ymin,xmax,ymax) of the data included in the graph."
        xmin, ymin, xmax, ymax = self.datasets[0].GetWorld()
        for dset in self.datasets:
            if dset.xmin < xmin: xmin = dset.xmin
            if dset.ymin < ymin: ymin = dset.ymin
            if dset.xmax > xmax: xmax = dset.xmax
            if dset.ymax > ymax: ymax = dset.ymax
        return xmin, ymin, xmax, ymax

    def SetView(self, xmin='', xmax='', ymin='', ymax=''):
        "Sets the graph view (defining the position of the graph on the canvas)."
        if xmin != '':
            self.view[0] = float(xmin)
        if xmax != '':
            self.view[2] = float(xmax)
        if ymin != '':
            self.view[1] = float(ymin)
        if ymax != '':
            self.view[3] = float(ymax)
        # Default position of legend box is upper left corner of the graph
        [xmin, ymin, xmax, ymax] = self.view
        self.SetLegend(showlegend=self.legend, xpos=xmin+0.01*(xmax-xmin), ypos=ymax-0.01*(ymax-ymin))

    def Add2View(self, xmin='', xmax='', ymin='', ymax=''):
        "Adds value(s) to the graph view."
        if xmin != '':
            self.view[0] += float(xmin)
        if xmax != '':
            self.view[2] += float(xmax)
        if ymin != '':
            self.view[1] += float(ymin)
        if ymax != '':
            self.view[3] += float(ymax)
        self.SetView()

    def ShowLegend(self):
        "Sets the graph legend flag to True."
        self.SetLegend(showlegend=True)

    def HideLegend(self):
        "Sets the graph legend flag to False."
        self.SetLegend(showlegend=False)

    def SetLegend(self, showlegend='', xpos='', ypos='', boxLstyle=0, fillpattern=0, Csize=1.0):
        "Sets the graph legend flag and its position/style."
        if showlegend != '':
            self.legend = bool(showlegend)
            self.string += '@ legend %s \n'%flag[showlegend]
        if xpos != '' and ypos != '':
            self.string += '@ legend %.8f, %.8f \n'%(xpos, ypos)
        if xpos != '' and ypos == '':
            self.string += '@ legend %.8f, %.8f \n'%(xpos, 0.80)
        if xpos == '' and ypos != '':
            self.string += '@ legend %.8f, %.8f \n'%(0.20, ypos)
        if boxLstyle != '':
            self.string += '@ legend box linestyle %i \n'%boxLstyle
        if fillpattern != '':
            self.string += '@ legend box fill pattern %i \n'%fillpattern
        if Csize != '':
            self.string += '@ legend char size %.8f \n'%Csize

    def SetXaxis(self, vmin='', vmax='', label='', labelsize='', labelautopos='', labelpospar='',
                 labelposper='', majorUnit='', minorUnit='', useticks='',
                 useticklabels='', ticklabelsize='', majorGridlines='', minorGridlines='',
                 autoscale='', scale=''):
        """
        Sets the following axis parameters (only those specified):

        Parameters
        ----------

        vmin/vmax : float, optional
        label : string, optional
        labelsize : float, optional
        labelautopos : bool, optional
        labelpospar : float, optional
        labelposper : float, optional
        majorUnit : float, optional
        minorUnit : float, optional
        useticks : bool, optional
        useticklabels : bool, optional
        ticklabelsize : float, optional
        majorGridlines : bool, optional
        minorGridlines : bool, optional
        autoscale : bool, optional
        scale : \"Normal\"/\"Logarithmic\", optional
        """
        if vmin != '':
            self.SetWorld(xmin=float(vmin))
        if vmax != '':
            self.SetWorld(xmax=float(vmax))
        if label != '':
            self.string += '@ xaxis label \"%s\" \n'%str(label)
        if labelsize != '':
            self.string += '@ xaxis label char size %.3f \n'%float(labelsize)
        if labelautopos != '':
            if labelautopos == True:
                self.string += '@ xaxis label place auto \n'
            else:
                self.string += '@ xaxis label place spec \n'
            if labelpospar != '' and labelposper != '':
                self.string += '@ xaxis label place %.8f, %.8f \n'%(float(labelpospar), float(labelposper))
        if majorUnit != '':
            self.string += '@ xaxis tick major %.8f \n'%float(majorUnit)
        if minorUnit != '':
            self.string += '@ xaxis tick minor %.8f \n'%float(minorUnit)
        if useticks != '':
            self.string += '@ xaxis tick %s \n'%flag[useticks]
        if useticklabels != '':
            self.string += '@ xaxis ticklabel %s \n'%flag[useticklabels]
        if ticklabelsize != '':
            self.string += '@ xaxis ticklabel char size %.3f \n'%float(ticklabelsize)
        if majorGridlines != '':
            self.string += '@ xaxis tick major grid %s \n'%flag[majorGridlines]
        if minorGridlines != '':
            self.string += '@ xaxis tick minor grid %s \n'%flag[minorGridlines]
        if autoscale != '':
            xmin, ymin, xmax, ymax = self.GetWorldOfDatasets()
            rng = xmax-xmin
            if rng > 1e-20:
                unit = 2*10**math.floor(math.log(rng, 10))
            else:
                unit = 1
            self.SetXaxis(vmin=xmin, vmax=xmax, majorUnit=unit, minorUnit=unit/2.)
        if scale == 'Logarithmic':
            self.string += '@ xaxes scale %s \n'%scale
            # Scale axis to positive numbers
            xmin, ymin, xmax, ymax = self.GetWorldOfDatasets()
            if xmin <= 0.: self.SetXaxis(vmin=1e-10)
            else: self.SetXaxis(vmin=xmin)
            if xmax <= 0.: self.SetXaxis(vmax=1e10)
            else: self.SetXaxis(vmax=xmax)

    def SetYaxis(self, vmin='', vmax='', label='', labelsize='', labelautopos='',
                 labelpospar='', labelposper='', majorUnit='', minorUnit='', useticks='',
                 useticklabels='', ticklabelsize='', majorGridlines='', minorGridlines='',
                 autoscale='', scale=''):
        """
        Sets the following axis parameters (only those specified):

        Parameters
        ---------------------

        vmin/vmax : float, optional
        label : string, optional
        labelsize : float, optional
        labelautopos : bool, optional
        labelpospar : float, optional
        labelposper : float, optional
        majorUnit : float, optional
        minorUnit : float, optional
        useticks : bool, optional
        useticklabels : bool, optional
        ticklabelsize : float, optional
        majorGridlines : bool, optional
        minorGridlines : bool, optional
        autoscale : bool, optional
        scale : \"Normal\"/\"Logarithmic\", optional
        """
        if vmin != '':
            self.SetWorld(ymin=vmin)
        if vmax != '':
            self.SetWorld(ymax=vmax)
        if label != '':
            self.string += '@ yaxis label \"%s\" \n'%str(label)
        if labelsize != '':
            self.string += '@ yaxis label char size %.3f \n'%float(labelsize)
        if labelautopos != '':
            if labelautopos == True:
                self.string += '@ yaxis label place auto \n'
            else:
                self.string += '@ yaxis label place spec \n'
            if labelpospar != '' and labelposper != '':
                self.string += '@ yaxis label place %.8f, %.8f \n'%(float(labelpospar), float(labelposper))
        if majorUnit != '':
            self.string += '@ yaxis tick major %.8f \n'%float(majorUnit)
        if minorUnit != '':
            self.string += '@ yaxis tick minor %.8f \n'%float(minorUnit)
        if useticks != '':
            self.string += '@ yaxis tick %s \n'%flag[useticks]
        if useticklabels != '':
            self.string += '@ yaxis ticklabel %s \n'%flag[useticklabels]
        if ticklabelsize != '':
            self.string += '@ yaxis ticklabel char size %.3f \n'%float(ticklabelsize)
        if majorGridlines != '':
            self.string += '@ yaxis tick major grid %s\n'%flag[majorGridlines]
        if minorGridlines != '':
            self.string += '@ yaxis tick minor grid %s \n'%flag[minorGridlines]
        if autoscale != '':
            xmin, ymin, xmax, ymax = self.GetWorldOfDatasets()
            rng = ymax-ymin
            if rng > 1e-20:
                unit = 2* 10**math.floor(math.log(rng, 10))
            else:
                unit = 1
            self.SetYaxis(vmin=ymin, vmax=ymax, majorUnit=unit, minorUnit=unit/2.)
        if scale == 'Logarithmic':
            self.string += '@ yaxes scale %s \n'%scale
            # Scale axis to positive numbers
            xmin, ymin, xmax, ymax = self.GetWorldOfDatasets()
            if ymin <= 0.: self.SetYaxis(vmin=1e-10)
            else: self.SetYaxis(vmin=ymin)
            if ymax <= 0.: self.SetYaxis(vmax=1e10)
            else: self.SetYaxis(vmax=ymax)

    def SetXaxisSpecialTicks(self, ticklist):
        "Sets the axis ticks according to the ticklist [[value0,label0],[value1,label1],...]."
        self.string += '@ xaxis tick spec type both\n@ xaxis tick spec 11\n'
        for i in range(len(ticklist)):
            x, lab = ticklist[i][0], ticklist[i][1]
            if lab in symbols:
                # replace with xmgr code for symbol
                lab = symbols[lab]
            self.string += '@ xaxis tick major %i, %.8f \n'%(i, x)
            self.string += '@ xaxis ticklabel %i, \"%s\" \n'%(i, lab)

    def SetYaxisSpecialTicks(self, ticklist):
        "Sets the axis ticks according to the ticklist [[value0,label0],[value1,label1],...]."
        self.string += '@ yaxis tick spec type both\n@ yaxis tick spec 11\n'
        for i in range(len(ticklist)):
            y, lab = ticklist[i][0], ticklist[i][1]
            if lab in symbols:
                # replace with xmgr code for symbol
                lab = symbols[lab]
            self.string += '@ yaxis tick major %i, %.8f \n'%(i, y)
            self.string += '@ yaxis ticklabel %i, \"%s\" \n'%(i, lab)

    def SetSpecial(self, string):
        "Inserts a special GRACE command string to the graph."
        if string[-1] != '\n': string += '\n'
        self.string += string

    def AddDatasets(self, *datasets):
        "Append data set object(s) to graph."
        for dset in datasets:
            if isinstance(dset, Dataset):
                self.datasets.append(dset)
            elif isinstance(dset, tuple):
                self.AddDatasets(*dset)
            elif isinstance(dset, list):
                self.AddDatasets(*tuple(dset))

    def GetXMGRstring(self, graphnr):
        "Returns a string containing the GRACE commands that describe the graph and its data sets."
        string = '# ==================== GRAPH %i ===================== \n'%graphnr
        string += '@ g%i on\n@ g%i type XY\n@ with g%i\n'%(graphnr, graphnr, graphnr)
        string += '@ world %.8f, %.8f, %.8f, %.8f\n'%tuple(self.world)
        string += '@ view %.8f, %.8f, %.8f, %.8f\n'%tuple(self.view)
        string += self.string
        for setnr in range(len(self.datasets)):
            dset = self.datasets[setnr]
            string += dset.GetXMGRstring(graphnr, setnr)
        return string


class Plot(object):

    """
    Class for a plot (containing graphs)
    """

    def __init__(self, filename='Default.agr', *graphs):
        "Returning instance of class."
        self.filename = filename
        self.graphs = []
        self.frame = [0.15, 0.15, 1.15, 0.85]
        self.string = '# Grace file created on %s\n'%time.ctime()
        self.SetDefaults()
        if len(graphs) > 0:
            self.AddGraphs(*graphs)

    def AddGraphs(self, *graphs):
        "Append graph object(s) to plot."
        for graph in graphs:
            print('Graph %i appended to plot \"%s\"'%(len(self.graphs), self.filename))
            self.graphs.append(graph)

    def ArrangeGraphs(self, nx=1, ny=1, vspace=0.1, hspace=0.1, order=0):
        "Arranges the graphs on the plot frame with nx rows and ny columns."
        print('Arranging graphs in format %ix%i using type %i ordering'%(nx, ny, order))
        x0, y0 = self.frame[0], self.frame[1]
        x1, y1 = self.frame[2], self.frame[3]
        dx = ((self.frame[2]-self.frame[0])-(nx-1)*hspace)/nx
        dy = ((self.frame[3]-self.frame[1])-(ny-1)*vspace)/ny
        for i in range(len(self.graphs)):
            graph = self.graphs[i]
            if order == 1:
                print('... %i: Filling top to down, starting from upper left'%i)
                x = x0+(i/ny)*(dx+hspace)
                y = y1-dy-(i%ny)*(dy+vspace)
            elif order == 2:
                print('... %i: Filling left to right, starting from lower left'%i)
                x = x0+(i%nx)*(dx+hspace)
                y = y0+(i/nx)*(dy+vspace)
            elif order == 3:
                print('... %i: Sorting down to top, starting from lower left'%i)
                x = x0+(i/ny)*(dx+hspace)
                y = y0+(i%ny)*(dy+vspace)
            elif order == 4:
                print('... %i: Sorting right to left, starting from upper right'%i)
                x = x1-dx-(i%nx)*(dx+hspace)
                y = y1-dy-(i/nx)*(dy+vspace)
            elif order == 5:
                print('... %i: Sorting top to down, starting from upper right'%i)
                x = x1-dx-(i/ny)*(dx+hspace)
                y = y1-dy-(i%ny)*(dy+vspace)
            elif order == 6:
                print('... %i: Sorting right to left, starting from lower right'%i)
                x = x1-dx-(i%nx)*(dx+hspace)
                y = y0+(i/nx)*(dy+vspace)
            elif order == 7:
                print('... %i: Sorting down to top, starting from lower right'%i)
                x = x1-dx-(i/ny)*(dx+hspace)
                y = y0+(i%ny)*(dy+vspace)
            else:
                print('...%i: Filling left to right, starting from upper left (default)'%i)
                x = x0+(i%nx)*(dx+hspace)
                y = y1-dy-(i/nx)*(dy+vspace)
            graph.view = [x, y, x+dx, y+dy]
            graph.SetView()

    def SetAxesFontSizes(self, labelsize=1.0, ticklabelsize=1.0):
        "Sets font sizes on all axis labels in the plot."
        for g in self.graphs:
            g.SetXaxis(labelsize=labelsize, ticklabelsize=ticklabelsize)
            g.SetYaxis(labelsize=labelsize, ticklabelsize=ticklabelsize)

    def DefineColor(self, colornr, RGB, description='UserDefined'):
        "Defines an arbitrary RGB color for use in the plot."
        R, G, B = tuple(RGB)
        self.string += '@map color %i to (%i, %i, %i), \"%s\"\n' %(colornr, R, G, B, description)

    def SetDefaults(self, backgroundcolor=0, backgroundfill=False,
                    linewidth=1.0, linestyle=1, color=1, pattern=1,
                    font=0, charsize=1.0, symbolsize=1.0):
        "Sets the defaults for the plot."
        self.string += '@map font 0 to "Times-Roman","Times-Roman"\n'
        self.string += '@autoscale onread none\n'
        self.string += '@background color %i\n' %backgroundcolor
        self.string += '@page background fill %s\n' %flag[backgroundfill]
        self.string += '@default linewidth %.3f\n' %linewidth
        self.string += '@default linestyle %i\n' %linestyle
        self.string += '@default color %i\n' %color
        self.string += '@default pattern %i\n' %pattern
        self.string += '@default font %i\n' %font
        self.string += '@default char size %.6f\n' %charsize
        self.string += '@default symbol size %.6f\n' %symbolsize

    def SetSpecial(self, string):
        "Inserts a special GRACE command string to the plot."
        if string[-1] != '\n': string += '\n'
        self.string += string

    def PutText(self, string, posx, posy, color=1, rot=0, font=0, just=0, charsize=1.0):
        "Inserts text in plot at specified position"
        self.string += '@with string\n'
        self.string += '@    string on\n'
        self.string += '@    string loctype view\n'
        self.string += '@    string %.6f, %.6f\n'%(posx, posy)
        self.string += '@    string color %i\n'%color
        self.string += '@    string rot %i\n'%rot
        self.string += '@    string font %i\n'%font
        self.string += '@    string just %i\n'%just
        self.string += '@    string char size %.6f\n'%charsize
        self.string += '@    string def "%s"\n'%string

    def ShowTimestamp(self, xpos=0.03, ypos=0.03, color=1, rot=0, font=0, charsize=0.75):
        "Inserts a timestamp into the plot."
        self.string += '@ timestamp on\n'
        self.string += '@ timestamp %.3f, %.3f\n'%(xpos, ypos)
        self.string += '@ timestamp color %i\n'%color
        self.string += '@ timestamp rot %i\n'%rot
        self.string += '@ timestamp font %i\n'%font
        self.string += '@ timestamp char size %.3f\n'%charsize

    def WriteFile(self, filename=''):
        "Writes the GRACE commands that describe the whole plot to a file."
        if filename == '': filename = self.filename
        print('Writing plot to file \"%s\"'%(filename))
        f = open(filename, 'w')
        f.write(self.string)
        for i in range(len(self.graphs)):
            graph = self.graphs[i]
            f.write(graph.GetXMGRstring(i))
        f.close()

    def Print2File(self, printfile, path2grace='xmgrace'):
        """
        Prints the GRACE plot to a file [PS/EPS/MIF/SVG/PDF/PNM/JPEG/PNG/Metafile]
        according to the file extension (provided that the necessary device library
        is installed).
        """
        dev = '-hdevice '
        if printfile.endswith('.ps'): dev = ''
        if printfile.endswith('.eps'): dev += 'EPS'
        if printfile.endswith('.mif'): dev += 'MIF'
        if printfile.endswith('.svg'): dev += 'SVG'
        if printfile.endswith('.pdf'): dev += 'PDF'
        if printfile.endswith('.pnm'): dev += 'PNM'
        if printfile.endswith('.jpg') or printfile.endswith('.jpeg'): dev += 'JPEG'
        if printfile.endswith('.png'): dev += 'PNG'
        if printfile.endswith('.gmf'): dev += 'Metafile'
        os.system('%s %s -hardcopy %s -printfile %s' %(path2grace, self.filename, dev, printfile))
        print('Printing plot to file \"%s\"'%printfile)

    def Print(self, printcmd='lp -c', path2grace='xmgrace'):
        "Sends plot to printer."
        tmpfile = 'tmp.xmgr.ps'
        self.Print2File(tmpfile, path2grace)
        os.system(printcmd+tmpfile)
        print('Plot send to printer (\"%s %s\")'%(printcmd, tmpfile))


def demo():
    """
    Generation of a demo plot
    """
    # Arbitrary data sets
    x = [0.5, 1.5, 3, 5]
    y1 = [2, 0.5, 5, 9.2]
    y2 = [0.75, 5, 4, 2]
    dy = [0.1, 0.2, 0.3, 0.4]

    # Create instances of data set
    d1 = XYDYset(x, y1, dy, legend='A', EBsize=0.5, EBlinewidth=2.0, EBcolor=4, Stype=1)
    d2 = XYset(x, y2, legend='B', Lwidth=3)
    d3 = XYset(x, y2, legend='C', Ltype=2, Lstyle=3, Lwidth=2, Lcolor=99,
                 Stype=1, Ssize=1.1, Scolor=3, Sfillcolor=4, Sfillpattern=1,
                 Slinewidth=3, Slinestyle=1)
    dsets = Array2XYsets(10.0*N.transpose(N.array([x, y1, y2])))

    # Create instance of a graph
    g = Graph(d1, d2, d3)
    g.SetXaxis(label='x-label', autoscale=True)
    g.SetYaxis(label='y-label', autoscale=True)
    g.SetTitle('Demo of Inelastica.io.xmgrace module', size=1.1)
    g.SetSubtitle('Graph subtitle')

    # Create instance of a plot
    p = Plot('demo.agr', g)
    p.DefineColor(99, [200, 10, 200])

    # Add more graphs to plot
    for i in range(1, 9):
        g = Graph(d1, d2)
        g.SetSubtitle('nr %i'%i)
        g.SetXaxis(autoscale=True)
        g.SetYaxis(autoscale=True)
        if i == 1:
            # Show legend
            g.ShowLegend()
        if i == 2:
            # Use special labels and ticks
            g.SetXaxisSpecialTicks([[2, '\alpha'], [3, '\beta'], [4, '\gamma']])
            g.SetYaxisSpecialTicks([[3, '\Gamma'], [5, '\Delta']])
        if i == 3:
            # No ticks and labels
            g.SetXaxis(useticks=False, useticklabels=False)
            g.SetYaxis(useticks=False, useticklabels=False)
        if i == 4:
            # MajorGridlines
            g.SetXaxis(majorGridlines=True, minorGridlines=True, scale='Logarithmic')
            g.SetYaxis(majorGridlines=True, minorGridlines=True, scale='Logarithmic')
        if i == 5:
            # Fixed coordinate ranges
            g.ShowLegend()
            g.SetXaxis(vmin=-4, vmax=4)
            g.SetYaxis(vmin=-2, vmax=8)
            g.AddDatasets(d3)
        if i == 6:
            g.AddDatasets(dsets)
            g.SetXaxis(autoscale=True)
            g.SetYaxis(autoscale=True)
        p.AddGraphs(g)

    # Arrange the plots in a 3x3 structure
    p.ArrangeGraphs(nx=3, ny=3, hspace=0.1, vspace=0.1)

    # Change font size on all graphs
    p.SetAxesFontSizes(labelsize=0.9, ticklabelsize=0.9)

    # Change size of last graph
    g.Add2View(xmin=-0.05, xmax=0.05, ymin=0.05, ymax=0.05)

    # Finally, write the plot file
    p.ShowTimestamp()
    p.WriteFile()
    p.Print2File('demo.eps')
    #p.Print()

if __name__ == '__main__':
    demo()
