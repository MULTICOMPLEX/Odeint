import numpy as np
# Interpolation Functions
from scipy.interpolate import RectBivariateSpline, interp1d, interp2d
# Plotting
import matplotlib.pyplot as plt
# Debugging libraries
import sys
import os
# Import time to time the simulation
import time

from bokeh.io import push_notebook
from bokeh.layouts import gridplot
from bokeh.layouts import row, widgetbox
from bokeh.plotting import figure, show, output_notebook
from bokeh.models import Range1d,CustomJS,Slider,Span,ColumnDataSource
from scipy.interpolate import interp1d
from bokeh.models.widgets import Slider
from ipywidgets import interact, widgets, FloatSlider

def generateMu(n,a0):
    N       = int(n*(n+1)/2)
    mu      = np.zeros((3,N))
    mu[2,:] = -1 # Set all values to negative 
    d       = 2*a0/(n+1) # vertical/horizontal distance
    mu[0,0] = -a0+1.5*d # y - value
    mu[1,0] = -a0+0.5*d # x - value
    h=mu[0,0]
    count=0
    for i in range(1,n):
        h=h+d
        for j in range(i+1):
            count=count+1
            mu[0,count]=h
            if (np.abs(mu[0,count])<10**-10):
                mu[0,count]=0
            if j==0:
                mu[1,count]=mu[1,0]
            else:
                mu[1,count]=mu[1,count-1]+d    
            if (np.abs(mu[1,count])<10**-10):
                mu[1,count]=0
    return mu

a0=.2
n=50
mu=generateMu(n,a0)


p1 = figure(title="Hysterons Distribution")
p1.grid.grid_line_alpha=0.3
p1.xaxis.axis_label = '\u0251'
p1.yaxis.axis_label = '\u03B2'
p1.x_range = Range1d(-a0, a0)
p1.y_range = Range1d(-a0, a0)
colors_va=[]
for i in range(len(mu[1,:])):
    colors_va = np.append(colors_va,["navy"])
    
pdata = {'x_val': mu[1,:],
        'y_val': mu[0,:],
        'colors': colors_va}

psource = ColumnDataSource(data=pdata)
r0=p1.circle('x_val', 'y_val', size=1, color='colors', alpha=0.5,source=psource)
p1.legend.location = "bottom_right"




s2 = ColumnDataSource(data=dict(x=[-a0, a0], ym=[-a0, a0]))
r1=p1.line(x='x', y='ym', color="black", line_width=1, alpha=0.6, source=s2)
s3 = ColumnDataSource(data=dict(x=[-a0, a0], ym=[-a0, -a0]))
r2=p1.line(x='x', y='ym', color="orange", line_width=5, alpha=0.6, source=s3)

p2 = figure(title="Cycle")
p2.grid.grid_line_alpha = 0.3
p2.xaxis.axis_label = 'Input'
p2.yaxis.axis_label = 'Output'
p2.x_range = Range1d(-a0, a0)
p2.y_range = Range1d(-(n*(n+1)/2), (n*(n+1)/2))
r3=p2.line([-a0], [-a0*n],line_width=2)
p2.legend.location = "top_left"

# output notebook
output_notebook()

def update(value):
    if     ((r2.data_source.data['x'][1]==a0) & (r2.data_source.data['x'][0]==-a0)):
        vert = 0
        hor  = 1
    elif   ((r2.data_source.data['ym'][1]==a0) & (r2.data_source.data['ym'][0]==-a0)):
        vert = 1
        hor  = 0
    if ((r2.data_source.data['ym'][0]<=value) & hor):   # change postive and hor
        r2.data_source.data['ym'] = [value,value]      # set location of horizontal line
    elif ((r2.data_source.data['ym'][0]>value) & hor): # change negative and hor
        r2.data_source.data['ym'] = [-a0,a0]           # change to vertical line
        r2.data_source.data['x']  = [value,value]      # set location
    elif ((r2.data_source.data['x'][0]>=value) & vert):# change negative and vertical 
        r2.data_source.data['x']  = [value,value]      # update location 
    elif ((r2.data_source.data['x'][0]<value) & vert):# change negative and vertical 
        r2.data_source.data['ym'] = [value,value]      # set location of horizontal line
        r2.data_source.data['x']  = [-a0,a0]           # change to horizontal line
    if hor:
        for i,val in enumerate(mu[0,:]):
            if (value>val):
                mu[2,i]=1
                r0.data_source.data['colors'][i] = "pink"
    if vert:
        for i,val in enumerate(mu[1,:]):
            if (value<val):
                mu[2,i]=-1
                r0.data_source.data['colors'][i] = "navy"
    newx = np.append(r3.data_source.data['x'],value)
    newy = np.append(r3.data_source.data['y'],np.sum(mu[2,:]))
    r3.data_source.data = {'x' : newx, 'y' : newy}
    push_notebook()
    
show(gridplot([[p1,p2]], plot_width=300, plot_height=300))

val_widget = FloatSlider(min=-a0, max=a0, step=a0/20, value=-a0)


interact(update,value=val_widget)