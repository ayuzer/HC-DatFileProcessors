from plotly.offline import download_plotlyjs, init_notebook_mode, iplot, plot
from plotly.graph_objs import *
init_notebook_mode()

import numpy as np
import csv

filename = 'background(all-48hr-CW20k-CD7573)ch-1-2-time-ribbon-calibrated.csv' #csv file name
with open(filename, 'r') as f:
    spectra = np.loadtxt(f, delimiter = ',')

y_raw = spectra[:,0] # Energy range
sample_size = spectra.shape[1]-1 #time region sample slices (# of columns except 1st col)
traces =[]

fill_colors = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854','#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854','#66c2a5', '#fc8d62']
tinterval= ['below -266 ns', '-266 : -133 ns', '-133 : 0 ns', '0 : 133 ns', '133 : 266 ns', '266 : 399 ns' , '399 : 533 ns' , '533 : 666 ns', '666 : 799 ns', '799 : 933 ns', '933 : 1066 ns', 'above 1066 ns' ]
for i in range(1,sample_size+1):#iterate through columns
    z_raw = spectra[:,i] #grab all rows energy histrogram counts for each time frame
    x = []
    y = []
    z = []

    ci = int(255/sample_size*i) #color index

    for j in range(0, len(z_raw)): #iterate through energy range row
        z.append(z_raw[j]) #counts
        y.append( y_raw[j] ) #loop of energy
        x.append( (i*10) ) #width of ribbon

    traces.append(dict(
        type='scatter3d',
        mode='lines',
        z=z,
        x=x,
        y=y,
        surfaceaxis=-1,
        #showlegend= False,
        name= tinterval[i-1],
        line=dict(
            color=z,
            colorscale='Viridis',
            width=4
        ),
    ))


axisrange = max(y_raw)+100

layout = Layout(
     title= 'HPGE - Top Veto Coincidence (HPGE energy spectrum through various time window)',
     scene = dict(
         aspectmode='cube',
         xaxis=dict(title='Delta Time Window',showgrid=False,showticklabels=False),
         yaxis=dict(title='HPGE energy (KeV)',range=[0,axisrange]),
         zaxis=dict(title='Counts')
                 )
         )

fig = dict(data=traces, layout = layout)

name = filename.replace(".csv", "")
newname = name + "-plot.html"
plot(fig, filename=newname)

