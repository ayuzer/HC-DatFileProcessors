from plotly.offline import download_plotlyjs, init_notebook_mode, iplot, plot
from plotly.graph_objs import *
init_notebook_mode()

import numpy as np
import csv

filename = '2Dcoin-yellowknife.csv' #csv file name
with open(filename, 'r') as f:
    data_iter= csv.reader(f, delimiter=',')
    
    filedata = list(data_iter)
    fCoinArray= np.array(filedata) #contains 2d Coin array
    CoinArray= fCoinArray.astype(np.float)

x=np.array(CoinArray[:,0])
y=np.array(CoinArray[:,1])
z=np.array(CoinArray[:,2])

axisrange = max(np.amax(x),np.amax(y))

trace0 = Scatter3d( #points
    x=x,
    y=y, 
    z=z,
    mode = 'markers',
    
    marker = dict(
        size=4,
        color=z,
        colorscale='Viridis',
        showscale=True,
        opacity=0.8
        )
    )
               

data= [trace0] #, trace1]

layout = Layout(
    title= '2D coincidence plot',
    scene = dict(
        aspectmode='cube',
        bgcolor='rgb(227,252,220)',
        xaxis=dict(title='CH0 energy', range=[0,axisrange]),
        yaxis=dict(title='CH1 energy',range=[0,axisrange]),
        zaxis=dict(title='Counts')
                )
        )

fig = dict(data=data, layout = layout)

name=filename.replace(".csv", "")
newname=name+"_3DCoinPlot_colorscale.html"

plot(fig, filename=newname)

#####################################################




