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
    CoinArray = CoinArray[np.argsort(CoinArray[:,2])] #sort by 3rd column (z value)

x=np.array(CoinArray[:,0]) #x coord
y=np.array(CoinArray[:,1]) #y coord
z=np.array(CoinArray[:,2]) # number of counts


axisrange = max(np.amax(x),np.amax(y)) + 100

trace0 = Scattergl( #points
    x=x,
    y=y, 
    mode = 'markers',
    text= z,
    
    marker = dict(
        size=10,
        color=z,
        colorscale='Viridis',
        opacity=0.8,
        showscale=True
        )
    )
               

data= [trace0] #, trace1]

layout = Layout(
    title= '2D coincidence plot',
    autosize='false',
    width='750',
    height='750',
    hovermode='closest',
    
    xaxis=dict(title='CH0 energy',range=[0,axisrange],nticks= 20),
    yaxis=dict(title='CH1 energy',range=[0,axisrange],nticks= 20)
)    
#    scene = dict(
#        
#        xaxis=dict(title='CH0 energy',range=[0,axisrange]),
#        yaxis=dict(title='CH1 energy',range=[0,axisrange])
#                )
#        )

fig = dict(data=data, layout = layout)

name=filename.replace(".csv", "")
newname=name+"_2DCoinPlot_colorscale.html"

plot(fig, filename=newname)