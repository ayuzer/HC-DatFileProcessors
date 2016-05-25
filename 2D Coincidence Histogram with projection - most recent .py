from plotly.offline import download_plotlyjs, init_notebook_mode, iplot, plot
from plotly.graph_objs import *
init_notebook_mode()

import numpy as np
import csv
import math

filename = 'UraniumStd(U-002)ch-0-1-doublecoin-data.csv' #csv file name
with open(filename, 'r') as f:
    data_iter= csv.reader(f, delimiter=',')
    
    filedata = list(data_iter)
    fCoinArray= np.array(filedata) #contains 2d Coin array
    CoinArray= fCoinArray.astype(np.float)


x=np.array(CoinArray[:,0]) #x coord
y=np.array(CoinArray[:,1]) #y coord


axisrange = max(np.amax(x),np.amax(y)) + 100

#number of bins
k= int(axisrange/2)

name = filename.replace(".csv", "")
plotname = name+" 2D Coincidence"

#removing trace1 (individual data points due to lag)
#trace1 = Scatter(
#    showlegend=False,x=x, y=y, mode='markers', name='points',
#    marker=dict(color='rgb(102,0,0)', size=3, opacity=0.2),
#    hoverinfo='none'
#)

trace2 = Histogram2d(
    x=x, y=y, name='density', nbinsx=k, nbinsy=k,
    colorscale=[[0, 'rgb(255,255,255)'],[0.01, 'rgb(0,0,255)'], [0.1, 'rgb(50,255,255)'], [0.25, 'rgb(0,255,0)'], [0.5, 'rgb(255,255,0)'], [0.75, 'rgb(255,0,0)'], [1, 'rgb(0,0,0)']], reversescale=False,zmin=3, zauto=False, zsmooth= 'best'
)
trace3 = Histogram(
    showlegend=False,x=x, name='HPGE 1 Energy (KeV)',
    marker=dict(color='rgb(102,0,0)'),
    nbinsx=500,
    nbinsy=500,
    yaxis='y2'
)
trace4 = Histogram(
    showlegend=False,y=y, name='HPGE 2 Energy (KeV)', marker=dict(color='rgb(102,0,0)'),
    nbinsx=500,
    nbinsy=500,
    xaxis='x2'
)
data = [trace2, trace3, trace4]

layout = Layout(
    title = plotname,
    autosize=False,
    width=800,
    height=750,
    xaxis=dict(
        #type = 'log', #set x acis to log
        title='HPGE 1 Energy (KeV)',
        domain=[0, 0.85],
        showgrid=False,
        range=[0,axisrange],nticks= 10
    ),
    yaxis=dict(
        #type = 'log', #set y axis to log
        title='HPGE 2 Energy (KeV)',
        domain=[0, 0.85],
        showgrid=False,
        range=[0,axisrange],nticks= 10
    ),
    margin=dict(
        t=50
    ),
    hovermode='closest',
    bargap=0,
    xaxis2=dict(
        domain=[0.85, 1],
        showgrid=False,
        zeroline=False
    ),
    yaxis2=dict(
        domain=[0.85, 1],
        showgrid=False,
        zeroline=False
    )
)



fig = dict(data=data, layout = layout)


newname = name + "-2D histogram Plot.html"
plot(fig, filename=newname)
