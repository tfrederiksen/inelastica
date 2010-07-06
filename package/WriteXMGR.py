"""
###########################################################################
#   Module WriteXMGR.py                                                   #
#   Written by Thomas Frederiksen, July 2007.                             #
#   Copyright (c), All Rights Reserved                                    #
#                                                                         #
#   This program is free software; you can redistribute it and/or modify  #
#   it under the terms of the GNU General Public License as published by  #
#   the Free Software Foundation; either version 3 of the License, or     #
#   (at your option) any later version.                                   #
#                                                                         #
#   This program is distributed in the hope that it will be useful,       #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#   GNU General Public License for more details.                          #
#                                                                         #
#   You should have received a copy of the GNU General Public License     #
#   along with this program.  If not, see <http://www.gnu.org/licenses/>. #
###########################################################################


This module provides a Python interface to XMGR/GRACE (get the newest
information about GRACE and download the latest version at
http://plasma-gate.weizmann.ac.il/Grace/).

To create a GRACE file from a given data set, the minimal script
structure will take the following form:

    >>> data = XYset(x,y)
    >>> graph = Graph(data)
    >>> plot = Plot('test.agr',graph)
    >>> plot.WriteFile()

...Done!

More graphs and data sets are easily appended. Continuing on
the example above:

    >>> data2 = XYset(x2,y2)
    >>> graph2 = Graph(data2)
    >>> plot.AddGraphs(graph2)
    >>> plot.ArrangeGraphs(nx=2,ny=1)
    >>> plot.WriteFile('test2.agr')

With these additional lines we end up with a plot containing
two graphs. The function \"ArrangeGraphs(...)\" is called
to arrange the graphs nicely on the canvas (side by side).

A more detailed demonstration of this module can be generated
by typing \"python WriteXMGR.py\" in the command promt. This
will create a file \"test.agr\".
"""

import time, math, os
import numpy as N

symbols = {'\alpha':'\\f{Symbol}a',
           '\beta':'\\f{Symbol}b',
           '\gamma':'\\f{Symbol}g',
           '\delta':'\\f{Symbol}d',
           '\Gamma':'\\f{Symbol}G',
           '\Delta':'\\f{Symbol}D'
          }
flag = {True:'on',False:'off'}


def Array2XYsets(A,**keywords):
    """
    Converts a matrix A = [[x0,y0(x0),y1(x0),...],[x1,y0(x1),y1(x1),...],...]
    into a list of XYset objects. Keywords for the format (Lcolor, Stype, etc.)
    are applied to all the generated XYsets.
    """
    sets = []
    A = N.array(A)
    x = A[:,0]
    for i in range(1,len(A[0])):
        y = A[:,i]
        sets.append(XYset(x,y,**keywords))
    return sets



class Dataset:
    
    def __init__(self,**keywords):
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

    def SetFormat(self,legend='',Ltype='',Lstyle='',Lwidth='',Lcolor='',
                  Stype='',Ssize='',Scolor='',Sfillcolor='',Sfillpattern='',
                  Slinewidth='',Slinestyle='',
                  EBsize='',EBcolor='',EBlinewidth=''):
        "Sets the formatting options (linestyle, symbol type, etc.)."
        if legend!='': self.SetLegend(legend)
        if Ltype!='': self.SetLtype(Ltype)
        if Lstyle!='': self.SetLstyle(Lstyle)
        if Lwidth!='': self.SetLwidth(Lwidth)
        if Lcolor!='': self.SetLcolor(Lcolor)
        if Stype!='': self.SetStype(Stype)
        if Ssize!='': self.SetSsize(Ssize)
        if Scolor!='': self.SetScolor(Scolor)
        if Sfillcolor!='': self.SetSfillcolor(Sfillcolor)
        if Sfillpattern!='': self.SetSfillpattern(Sfillpattern)
        if Slinewidth!='': self.SetSlinewidth(Slinewidth)
        if Slinestyle!='': self.SetSlinestyle(Slinestyle)
        if EBsize!='': self.SetEBsize(EBsize)
        if EBcolor!='': self.SetEBcolor(EBcolor)
        if EBlinewidth!='': self.SetEBlinewidth(EBlinewidth)

    def SetLegend(self,legend):
        "Sets the legend (string)."
        self.legend = legend

    def SetLtype(self,Ltype):
        "Sets the line type (integer)."
        self.Ltype = Ltype

    def SetLstyle(self,Lstyle):
        "Sets the line style (integer)."
        self.Lstyle = Lstyle

    def SetLwidth(self,Lwidth):
        "Sets the line width (float)."
        self.Lwidth = Lwidth
        
    def SetLcolor(self,Lcolor):
        "Sets the line color (integer)."
        self.Lcolor = Lcolor
        
    def SetStype(self,Stype):
        "Sets the symbol type (integer)."
        self.Stype = Stype

    def SetSsize(self,Ssize):
        "Sets the symbol size (float)."
        self.Ssize = Ssize
        
    def SetScolor(self,Scolor):
        "Sets the symbol color (integer)."
        self.Scolor = Scolor
        
    def SetSfillcolor(self,Sfillcolor):
        "Sets the symbol fill color (integer)."
        self.Sfillcolor = Sfillcolor

    def SetSfillpattern(self,Sfillpattern):
        "Sets the symbol fill pattern (integer)."
        self.Sfillpattern = Sfillpattern
        
    def SetSlinewidth(self,Slinewidth):
        "Sets the symbol linewidth (float)."
        self.Slinewidth = Slinewidth
        
    def SetSlinestyle(self,Slinestyle):
        "Sets the symbol line style (integer)."
        self.Slinestyle = Slinestyle

    def SetEBsize(self,EBsize):
        "Sets the error bar size (float)."
        self.EBsize = EBsize

    def SetEBcolor(self,EBcolor):
        "Sets the error bar color (integer)."
        self.EBcolor = EBcolor

    def SetEBlinewidth(self,EBlinewidth):
        "Sets the error bar linewidth (float)."
        self.EBlinewidth = EBlinewidth



    def __GetXMGRstring__(self,graphnr,setnr):
        string = ''
        if self.legend!='':
            string += '@ s%i legend \"%s\" \n'%(setnr,self.legend)
        if self.Ltype!='':
            string += '@ s%i line type %i \n'%(setnr,self.Ltype)
        if self.Lstyle!='':
            string += '@ s%i line linestyle %i \n'%(setnr,self.Lstyle)
        if self.Lwidth!='':
            string += '@ s%i line linewidth %f \n'%(setnr,self.Lwidth)
        if self.Lcolor!='':
            string += '@ s%i line color %i \n'%(setnr,self.Lcolor)
        if self.Stype!='':
            string += '@ s%i symbol %i \n'%(setnr,self.Stype)
        if self.Ssize!='':
            string += '@ s%i symbol size %f \n'%(setnr,self.Ssize)
        if self.Scolor!='':
            string += '@ s%i symbol color %i \n'%(setnr,self.Scolor)
        if self.Sfillcolor!='':
            string += '@ s%i symbol fill color %i \n'%(setnr,self.Sfillcolor)
        if self.Sfillpattern!='':
            string += '@ s%i symbol fill pattern %i \n'%(setnr,self.Sfillpattern)
        if self.Slinewidth!='':
            string += '@ s%i symbol linewidth %f \n'%(setnr,self.Slinewidth)
        if self.Slinestyle!='':
            string += '@ s%i symbol linestyle %i \n'%(setnr,self.Slinestyle)
        if self.EBsize!='':
            string += '@ s%i errorbar size %f \n'%(setnr,self.EBsize)
        if self.EBcolor!='':
            string += '@ s%i errorbar color %i \n'%(setnr,self.EBcolor)
        if self.EBlinewidth!='':
            string += '@ s%i errorbar linewidth %f \n'%(setnr,self.EBlinewidth)
            string += '@ s%i errorbar riser linewidth %f \n'%(setnr,self.EBlinewidth)
        return string

    
class XYset(Dataset):

    def __init__(self,x,y,**keywords):
        "Returning instance of class."
        Dataset.__init__(self,**keywords)
        self.x = N.array(x)
        self.y = N.array(y)
        if len(x) != len(y):
            print 'ERROR: x and y data lists are not of equal length!!!'
            error
        self.xmin = float(min(self.x))
        self.ymin = float(min(self.y))
        self.xmax = float(max(self.x))
        self.ymax = float(max(self.y))

    def GetWorld(self):
        "Returns tuple (xmin,ymin,xmax,ymax) of the data set."
        return self.xmin,self.ymin,self.xmax,self.ymax
    

    def GetXMGRstring(self,graphnr,setnr):
        "Returns a string containing the GRACE commands that describe the data set."
        string = '# ------------------ Dataset %i.%i ----------------- \n'%(graphnr,setnr)
        string += '@ target g%i.s%i \n@ type xy \n'%(graphnr,setnr)
        string += self.__GetXMGRstring__(graphnr,setnr)
        for i in range(len(self.x)):
            string += '   %.8f  %.8f \n'%(self.x[i],self.y[i])
        return string
    

class XYDXset(XYset):
    
    def __init__(self,x,y,dx,**keywords):
        "Returning instance of class."
        XYset.__init__(self,x,y,**keywords)
        self.dx = N.array(dx)

    def GetXMGRstring(self,graphnr,setnr):
        "Returns a string containing the GRACE commands that describe the data set."
        string = '# ------------------ Dataset %i.%i ----------------- \n'%(graphnr,setnr)
        string += '@ target g%i.s%i \n@ type xydx \n'%(graphnr,setnr)
        string += self.__GetXMGRstring__(graphnr,setnr)
        for i in range(len(self.x)):
            string += '   %.8f  %.8f  %.8f \n'%(self.x[i],self.y[i],self.dx[i])
        return string


class XYDYset(XYset):
   
    def __init__(self,x,y,dy,**keywords):
        "Returning instance of class."
        XYset.__init__(self,x,y,**keywords)
        self.dy = N.array(dy)

    def GetXMGRstring(self,graphnr,setnr):
        "Returns a string containing the GRACE commands that describe the data set."
        string = '# ------------------ Dataset %i.%i ----------------- \n'%(graphnr,setnr)
        string += '@ target g%i.s%i \n@ type xydy \n'%(graphnr,setnr)
        string += self.__GetXMGRstring__(graphnr,setnr)
        for i in range(len(self.x)):
            string += '   %.8f  %.8f  %.8f \n'%(self.x[i],self.y[i],self.dy[i])
        return string


class XYDXDYset(XYset):
   
    def __init__(self,x,y,dx,dy,**keywords):
        "Returning instance of class."
        XYset.__init__(self,x,y,**keywords)
        self.dx = N.array(dx)
        self.dy = N.array(dy)

    def GetXMGRstring(self,graphnr,setnr):
        "Returns a string containing the GRACE commands that describe the data set."
        string = '# ------------------ Dataset %i.%i ----------------- \n'%(graphnr,setnr)
        string += '@ target g%i.s%i \n@ type xydxdy \n'%(graphnr,setnr)
        string += self.__GetXMGRstring__(graphnr,setnr)
        for i in range(len(self.x)):
            string += '   %.8f  %.8f  %.8f  %.8f \n'%(self.x[i],self.y[i],self.dx[i],self.dy[i])
        return string



class Graph:

    def __init__(self,*datasets):
        "Returning instance of class."
        self.datasets = []
        self.world = [0.,0.,1.,1.]
        self.view = [0.15,0.15,1.15,0.85]
        self.string = ''
        self.legend = False
        self.SetWorld()
        self.SetView()
        if len(datasets)>0:
            self.AddDatasets(*datasets)
            
    def SetTitle(self,title,font=0,size=1.5,color=1):
        "Sets the title of the graph."
        self.string += '@ title \"%s\"\n'%title
        self.string += '@ title font %i \n'%font
        self.string += '@ title size %.8f \n'%size
        self.string += '@ title color %i \n'%color
    
    def SetSubtitle(self,subtitle,font=0,size=1.0,color=1):
        "Sets the subtitle of the graph."
        self.string += '@ subtitle \"%s\"\n'%subtitle
        self.string += '@ subtitle font %i \n'%font
        self.string += '@ subtitle size %.8f \n'%size
        self.string += '@ subtitle color %i \n'%color
    
    def SetWorld(self,xmin=False,xmax=False,ymin=False,ymax=False):
        "Sets the graph world (defining the axis limits)."
        if not xmin: xmin = self.world[0]
        if not xmax: xmax = self.world[2]
        if not ymin: ymin = self.world[1]
        if not ymax: ymax = self.world[3]
        if xmin < xmax:
          self.world[0] = xmin
          self.world[2] = xmax
        if ymin < ymax:
          self.world[1] = ymin
          self.world[3] = ymax 

    def GetWorldOfDatasets(self):
        "Returns tuple (xmin,ymin,xmax,ymax) of the data included in the graph."
        xmin,ymin,xmax,ymax = self.datasets[0].GetWorld()
        for set in self.datasets:
            if set.xmin < xmin: xmin = set.xmin
            if set.ymin < ymin: ymin = set.ymin
            if set.xmax > xmax: xmax = set.xmax
            if set.ymax > ymax: ymax = set.ymax
        return xmin,ymin,xmax,ymax
    
    def SetView(self,xmin='',xmax='',ymin='',ymax=''):
        "Sets the graph view (defining the position of the graph on the canvas)."
        if xmin!='':
            self.view[0] = float(xmin)
        if xmax!='':
            self.view[2] = float(xmax)
        if ymin!='':
            self.view[1] = float(ymin)
        if ymax!='':
            self.view[3] = float(ymax)
        # Default position of legend box is upper left corner of the graph
        [xmin,ymin,xmax,ymax] = self.view
        self.SetLegend(showlegend=self.legend,xpos=xmin+0.01*(xmax-xmin),ypos=ymax-0.01*(ymax-ymin))

    def Add2View(self,xmin='',xmax='',ymin='',ymax=''):
        "Adds value(s) to the graph view."
        if xmin!='':
            self.view[0] += float(xmin)
        if xmax!='':
            self.view[2] += float(xmax)
        if ymin!='':
            self.view[1] += float(ymin)
        if ymax!='':
            self.view[3] += float(ymax)
        self.SetView()

    def ShowLegend(self):
        "Sets the graph legend flag to True."
        self.SetLegend(showlegend=True)

    def HideLegend(self):
        "Sets the graph legend flag to False."
        self.SetLegend(showlegend=False)

    def SetLegend(self,showlegend='',xpos='',ypos='',boxLstyle=0,fillpattern=0,Csize=1.0):
        "Sets the graph legend flag and its position/style."
        if showlegend!='':
            self.legend = bool(showlegend)
            self.string += '@ legend %s \n'%flag[showlegend]
        if xpos!='' and ypos!='':
            self.string += '@ legend %.8f, %.8f \n'%(xpos,ypos)
        if boxLstyle!='':
            self.string += '@ legend box linestyle %i \n'%boxLstyle
        if fillpattern!='':
            self.string += '@ legend box fill pattern %i \n'%fillpattern
        if Csize!='':
            self.string += '@ legend char size %.8f \n'%Csize
    
    def SetXaxis(self,min='',max='',label='',labelsize='',labelautopos='',labelpospar='',
                 labelposper='',majorUnit='',minorUnit='',useticks='',
                 useticklabels='',ticklabelsize='',majorGridlines='',minorGridlines='',
                 autoscale='',scale=''):
        """
        Sets the axis parameters (if specified):
        - min/max        = [float]
        - label          = [string]
        - labelsize      = [float]
        - labelautopos   = [bool]
        - labelpospar    = [float]
        - labelposper    = [float]
        - majorUnit      = [float]
        - minorUnit      = [float]
        - useticks       = [bool]
        - useticklabels  = [bool]
        - ticklabelsize  = [float]
        - majorGridlines = [bool]
        - minorGridlines = [bool]
        - autoscale      = [bool]
        - scale          = [\"Normal\"/\"Logarithmic\"]
        """
        if min!='':
            self.SetWorld(xmin=float(min))
        if max!='':
            self.SetWorld(xmax=float(max))
        if label!='':
            self.string += '@ xaxis label \"%s\" \n'%str(label)
        if labelsize!='':
            self.string += '@ xaxis label char size %.3f \n'%float(labelsize)
        if labelautopos!='':
           if labelautopos==True:  
              self.string += '@ xaxis label place auto \n'
           else:
              self.string += '@ xaxis label place spec \n'
           if labelpospar!='' and labelposper!='':
              self.string += '@ xaxis label place %.8f, %.8f \n'%(float(labelpospar),float(labelposper))
        if majorUnit!='':
            self.string += '@ xaxis tick major %.8f \n'%float(majorUnit)
        if minorUnit!='':
            self.string += '@ xaxis tick minor %.8f \n'%float(minorUnit)
        if useticks!='':
            self.string += '@ xaxis tick %s \n'%flag[useticks]
        if useticklabels!='':
            self.string += '@ xaxis ticklabel %s \n'%flag[useticklabels]
        if ticklabelsize!='':
            self.string += '@ xaxis ticklabel char size %.3f \n'%float(ticklabelsize)
        if majorGridlines!='':
            self.string += '@ xaxis tick major grid %s \n'%flag[majorGridlines]
        if minorGridlines!='':
            self.string += '@ xaxis tick minor grid %s \n'%flag[minorGridlines]
        if autoscale!='':
            xmin,ymin,xmax,ymax = self.GetWorldOfDatasets()
            rng = xmax-xmin
            if rng>1e-20:
                unit = 0.2*10**math.floor(math.log(rng,10))
            else:
                unit = 1 
            self.SetXaxis(min=xmin,max=xmax,majorUnit=unit,minorUnit=unit/2.)
        if scale =='Logarithmic':
            self.string += '@ xaxes scale %s \n'%scale
            # Scale axis to positive numbers
            xmin,ymin,xmax,ymax = self.GetWorldOfDatasets()
            if xmin<=0.: self.SetXaxis(min=1e-10)
            else: self.SetXaxis(min=xmin)
            if xmax<=0.: self.SetXaxis(max=1e10)
            else: self.SetXaxis(max=xmax)
            
    def SetYaxis(self,min='',max='',label='',labelsize='',labelautopos='',
                 labelpospar='',labelposper='',majorUnit='',minorUnit='',useticks='',
                 useticklabels='',ticklabelsize='',majorGridlines='',minorGridlines='',
                 autoscale='',scale=''):
        """
        Sets the axis parameters (if specified):
        - min/max        = [float]
        - label          = [string]
        - labelsize      = [float]
        - labelautopos   = [bool]
        - labelpospar    = [float]
        - labelposper    = [float]
        - majorUnit      = [float]
        - minorUnit      = [float]
        - useticks       = [bool]
        - useticklabels  = [bool]
        - ticklabelsize  = [float]
        - majorGridlines = [bool]
        - minorGridlines = [bool]
        - autoscale      = [bool]
        - scale          = [\"Normal\"/\"Logarithmic\"]
        """
        if min!='':
            self.SetWorld(ymin=min)
        if max!='':
            self.SetWorld(ymax=max)
        if label!='':
            self.string += '@ yaxis label \"%s\" \n'%str(label)
        if labelsize!='':
            self.string += '@ yaxis label char size %.3f \n'%float(labelsize)
        if labelautopos!='':
           if labelautopos==True:  
              self.string += '@ yaxis label place auto \n'
           else:
              self.string += '@ yaxis label place spec \n'
           if labelpospar!='' and labelposper!='':
              self.string += '@ yaxis label place %.8f, %.8f \n'%(float(labelpospar),float(labelposper))
        if majorUnit!='':
            self.string += '@ yaxis tick major %.8f \n'%float(majorUnit)
        if minorUnit!='':
            self.string += '@ yaxis tick minor %.8f \n'%float(minorUnit)
        if useticks!='':
            self.string += '@ yaxis tick %s \n'%flag[useticks]
        if useticklabels!='':
            self.string += '@ yaxis ticklabel %s \n'%flag[useticklabels]
        if ticklabelsize!='':
            self.string += '@ yaxis ticklabel char size %.3f \n'%float(ticklabelsize)
        if majorGridlines!='':
            self.string += '@ yaxis tick major grid %s\n'%flag[majorGridlines]
        if minorGridlines!='':
            self.string += '@ yaxis tick minor grid %s \n'%flag[minorGridlines]
        if autoscale!='':
            xmin,ymin,xmax,ymax = self.GetWorldOfDatasets()
            rng = ymax-ymin
            if rng>1e-20:
                unit = 0.2*10**math.floor(math.log(rng,10))
            else:
                unit = 1
            self.SetYaxis(min=ymin,max=ymax,majorUnit=unit,minorUnit=unit/2.)
        if scale =='Logarithmic':
            self.string += '@ yaxes scale %s \n'%scale
            # Scale axis to positive numbers
            xmin,ymin,xmax,ymax = self.GetWorldOfDatasets()
            if ymin<=0.: self.SetYaxis(min=1e-10)
            else: self.SetYaxis(min=ymin)
            if ymax<=0.: self.SetYaxis(max=1e10)
            else: self.SetYaxis(max=ymax)

                    
    def SetXaxisSpecialTicks(self,ticklist):
        "Sets the axis ticks according to the ticklist [[value0,label0],[value1,label1],...]."
        self.string += '@ xaxis tick spec type both\n@ xaxis tick spec 11\n'
        for i in range(len(ticklist)):
            x,lab = ticklist[i][0],ticklist[i][1]
            try: lab = symbols[lab]
            except: pass
            self.string += '@ xaxis tick major %i, %.8f \n'%(i,x)
            self.string += '@ xaxis ticklabel %i, \"%s\" \n'%(i,lab)

    def SetYaxisSpecialTicks(self,ticklist):
        "Sets the axis ticks according to the ticklist [[value0,label0],[value1,label1],...]."
        self.string += '@ yaxis tick spec type both\n@ yaxis tick spec 11\n'
        for i in range(len(ticklist)):
            y,lab = ticklist[i][0],ticklist[i][1]
            try: lab = symbols[lab]
            except: pass
            self.string += '@ yaxis tick major %i, %.8f \n'%(i,y)
            self.string += '@ yaxis ticklabel %i, \"%s\" \n'%(i,lab)

    def SetSpecial(self,string):
        "Inserts a special GRACE command string to the graph."
        if string[-1]!='\n': string += '\n'
        self.string += string

    def AddDatasets(self,*datasets):
        "Append data set object(s) to graph."
        for set in datasets:
            self.datasets.append(set)

    def GetXMGRstring(self,graphnr):
        "Returns a string containing the GRACE commands that describe the graph and its data sets."
        string = '# ==================== GRAPH %i ===================== \n'%graphnr
        string += '@ g%i on\n@ g%i type XY\n@ with g%i\n'%(graphnr,graphnr,graphnr)
        string += '@ world %.8f, %.8f, %.8f, %.8f\n'%tuple(self.world)
        string += '@ view %.8f, %.8f, %.8f, %.8f\n'%tuple(self.view)
        string += self.string
        for setnr in range(len(self.datasets)):
            set = self.datasets[setnr]
            string += set.GetXMGRstring(graphnr,setnr)
        return string


class Plot:

    def __init__(self,filename='Default.agr',*graphs):
        "Returning instance of class."
        self.filename = filename
        self.graphs = []
        self.frame = [0.15,0.15,1.15,0.85]
        self.string = '# Grace file created on %s\n'%time.ctime()
        self.SetDefaults()
        if len(graphs)>0:
            self.AddGraphs(*graphs)

    def AddGraphs(self,*graphs):
        "Append graph object(s) to plot."
        for graph in graphs:
            print 'Graph %i appended to plot \"%s\"'%(len(self.graphs),self.filename)
            self.graphs.append(graph)

    def ArrangeGraphs(self,nx=1,ny=1,vspace=0.1,hspace=0.1,order=0):
        "Arranges the graphs on the plot frame with nx rows and ny columns."
        print 'Arranging graphs in format %ix%i using type %i ordering'%(nx,ny,order)
        x0,y0 = self.frame[0],self.frame[1]
        x1,y1 = self.frame[2],self.frame[3]
        dx = ((self.frame[2]-self.frame[0])-(nx-1)*hspace)/nx
        dy = ((self.frame[3]-self.frame[1])-(ny-1)*vspace)/ny
        for i in range(len(self.graphs)):
            graph = self.graphs[i]
            if order==1:
               print '... %i: Filling top to down, starting from upper left'%i
               x = x0+(i/ny)*(dx+hspace)
               y = y1-dy-(i%ny)*(dy+vspace)
            elif order==2:
               print '... %i: Filling left to right, starting from lower left'%i
               x = x0+(i%nx)*(dx+hspace)
               y = y0+(i/nx)*(dy+vspace)
            elif order==3:
               print '... %i: Sorting down to top, starting from lower left'%i
               x = x0+(i/ny)*(dx+hspace)
               y = y0+(i%ny)*(dy+vspace)
            elif order==4:
               print '... %i: Sorting right to left, starting from upper right'%i
               x = x1-dx-(i%nx)*(dx+hspace)
               y = y1-dy-(i/nx)*(dy+vspace)
            elif order==5:
               print '... %i: Sorting top to down, starting from upper right'%i
               x = x1-dx-(i/ny)*(dx+hspace)
               y = y1-dy-(i%ny)*(dy+vspace)
            elif order==6:
               print '... %i: Sorting right to left, starting from lower right'%i
               x = x1-dx-(i%nx)*(dx+hspace)
               y = y0+(i/nx)*(dy+vspace)
            elif order==7:
               print '... %i: Sorting down to top, starting from lower right'%i
               x = x1-dx-(i/ny)*(dx+hspace)
               y = y0+(i%ny)*(dy+vspace)        
            else:
               print '...%i: Filling left to right, starting from upper left (default)'%i
               x = x0+(i%nx)*(dx+hspace)
               y = y1-dy-(i/nx)*(dy+vspace)
            graph.view = [x,y,x+dx,y+dy]
            graph.SetView()

    def SetAxesFontSizes(self,labelsize=1.0,ticklabelsize=1.0):
        "Sets font sizes on all axis labels in the plot."
        for g in self.graphs:
            g.SetXaxis(labelsize=labelsize,ticklabelsize=ticklabelsize)
            g.SetYaxis(labelsize=labelsize,ticklabelsize=ticklabelsize)

    def DefineColor(self,colornr,RGB,description='UserDefined'):
        "Defines an arbitrary RGB color for use in the plot."
        R,G,B = tuple(RGB)
        self.string += '@map color %i to (%i, %i, %i), \"%s\"\n' %(colornr,R,G,B,description)

    def SetDefaults(self,backgroundcolor=0,backgroundfill=False,
                    linewidth=1.0,linestyle=1,color=1,pattern=1,
                    font=0,charsize=1.0,symbolsize=1.0):
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

    def SetSpecial(self,string):
        "Inserts a special GRACE command string to the plot."
        if string[-1]!='\n': string += '\n'
        self.string += string

    def ShowTimestamp(self,xpos=0.03,ypos=0.03,color=1,rot=0,font=0,charsize=0.75):
        "Inserts a timestamp into the plot."
        self.string += '@ timestamp on\n'
        self.string += '@ timestamp %.3f, %.3f\n'%(xpos,ypos)
        self.string += '@ timestamp color %i\n'%color
        self.string += '@ timestamp rot %i\n'%rot
        self.string += '@ timestamp font %i\n'%font
        self.string += '@ timestamp char size %.3f\n'%charsize

    def WriteFile(self,filename=''):
        "Writes the GRACE commands that describe the whole plot to a file."
        if filename=='': filename = self.filename
        print 'Writing plot to file \"%s\"'%(filename)
        f = open(filename,'w')
        f.write(self.string)
        for i in range(len(self.graphs)):
            graph = self.graphs[i]
            f.write(graph.GetXMGRstring(i))
        f.close()

    def Print2File(self,printfile,path2grace='xmgrace'):
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
        os.system('%s %s -hardcopy %s -printfile %s' %(path2grace,self.filename,dev,printfile))
        print 'Printing plot to file \"%s\"'%printfile

    def Print(self,printcmd='lp -c',path2grace='xmgrace'):
        "Sends plot to printer."
        tmpfile = 'tmp.xmgr.ps'
        self.Print2File(tmpfile,path2grace)
        os.system(printcmd+tmpfile)
        print 'Plot send to printer (\"%s %s\")'%(printcmd,tmpfile)

            
if __name__ == '__main__':
    print
    print '----------------------------------------------------'
    print 'Generating an example of the \"WriteXMGR.py\" module'
    print '----------------------------------------------------'
    print

    # Arbitrary data sets
    x = [0.5,1.5,3,5]
    y1 = [2,0.5,5,9.2]
    y2 = [0.75,5,4,2]
    dy = [0.1,0.2,0.3,0.4]

    # Create instances of data set
    d1 = XYDYset(x,y1,dy,legend='A',EBsize=0.5,EBlinewidth=2.0,EBcolor=4,Stype=1)
    d2 = XYset(x,y2,legend='B',Lwidth=3)
    d3 = XYset(x,y2,legend='C',Ltype=2,Lstyle=3,Lwidth=2,Lcolor=99,
                 Stype=1,Ssize=1.1,Scolor=3,Sfillcolor=4,Sfillpattern=1,
                 Slinewidth=3,Slinestyle=1)

    # Create instance of a graph
    g = Graph(d1,d2,d3)
    g.SetXaxis(label='x-label',autoscale=True)     
    g.SetYaxis(label='y-label',autoscale=True)
    g.SetTitle('Graph title',size=1.2)      
    g.SetSubtitle('Graph subtitle')

    # Create instance of a plot
    p = Plot('test.agr',g)
    p.DefineColor(99,[200,10,200])

    # Add more graphs to plot
    for i in range(1,9):
        g = Graph(d1,d2)
        g.SetSubtitle('nr %i'%i)
        g.SetXaxis(autoscale=True)
        g.SetYaxis(autoscale=True)
        if i==1:
            # Show legend
            g.ShowLegend()
        if i==2:
            # Use special labels and ticks
            g.SetXaxisSpecialTicks([[2,'\alpha'],[3,'\beta'],[4,'\gamma']])
            g.SetYaxisSpecialTicks([[3,'\Gamma'],[5,'\Delta']])
        if i==3:
            # No ticks and labels
            g.SetXaxis(useticks=False,useticklabels=False)
            g.SetYaxis(useticks=False,useticklabels=False)
        if i==4:
            # MajorGridlines
            g.SetXaxis(majorGridlines=True,minorGridlines=True,scale='Logarithmic')
            g.SetYaxis(majorGridlines=True,minorGridlines=True,scale='Logarithmic')
        if i==5:
            # Fixed coordinate ranges
            g.ShowLegend()
            g.SetXaxis(min=-4,max=4)
            g.SetYaxis(min=-2,max=8)
            g.AddDatasets(d3)
        p.AddGraphs(g)

    # Arrange the plots in a 3x3 structure
    p.ArrangeGraphs(nx=3,ny=3,hspace=0.1,vspace=0.1)

    # Change font size on all graphs
    p.SetAxesFontSizes(labelsize=0.9,ticklabelsize=0.9)
    
    # Change size of last graph
    g.Add2View(xmin=-0.05,xmax=0.05,ymin=0.05,ymax=0.05)

    # Finally, write the plot file
    p.ShowTimestamp()
    p.WriteFile()
    p.Print2File('test.eps')
    #p.Print()

