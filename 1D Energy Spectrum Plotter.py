from plotly.offline import download_plotlyjs, init_notebook_mode, iplot, plot
from plotly.graph_objs import *
init_notebook_mode()

import numpy as np
import csv

filename = 'background(all-48hr-CW20k-CD7573)ch-1-2-time-slices-align_to_ch1.csv' #csv file name
with open(filename, 'r') as f:
    spectra = np.loadtxt(f, delimiter = ',')

x = spectra[:,0] # Energy range
y = spectra[:,4] #first time slice

#for i in x:
#	x[i]*=3.779 #calib factor

axisrange = np.amax(x) + 100

trace0 = Scattergl( #points
    x=x,
    y=y,
    mode = 'lines')



data= [trace0] #, trace1]

layout = Layout(
    title= '1D energy spectrum (2D coin - time slice)',
    autosize='false',
    width='750',
    height='750',
    hovermode='closest',

    xaxis=dict(title='HPGE - Energy (KeV)',range=[0,axisrange],nticks= 20),
    yaxis=dict(title='Counts',nticks= 20,) #type='log')
)

fig = dict(data=data, layout = layout)

name=filename.replace(".csv", "")
newname=name+"_Spectrum(300 to 500 ns).html"

plot(fig, filename=newname)
