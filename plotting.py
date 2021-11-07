import pandas as pd
import numpy as np
import plotly.io as pio
import plotly.graph_objs as go
import plotly.express as px
import os
import glob
import re

pio.templates["colour"] = go.layout.Template(
    layout=go.Layout(
        colorway=['#3366CC', '#DC3912', '#FF9900', '#109618', '#990099', '#0099C6', '#DD4477', '#66AA00', '#B82E2E', '#316395']
    )
)
pio.templates.default = "plotly_white+presentation+colour"
os.chdir('Dropbox/Honours/Results/G9_flare/scaled/averaged_each/single_channel')

#G9
#6.7
G9_neg_027 = glob.glob('*G9*6.7*1[89]?.csv')

G9_epoch_6GHz = pd.read_csv(G9_neg_027[0], sep=',').loc[:,'Epoch'].to_numpy()
G9_neg_027_list = []
for i in range(len(G9_neg_027)):
    G9_neg_027_list.append(pd.read_csv(G9_neg_027[i], sep=',').loc[:,'Intensity'].to_numpy())

G9_neg_027_mean = np.mean(np.array(G9_neg_027_list), axis = 0)


G9_neg_085 = glob.glob('*G9*6.7*17[78].csv')

G9_neg_085_list = []
for i in range(len(G9_neg_085)):
    G9_neg_085_list.append(pd.read_csv(G9_neg_085[i], sep=',').loc[:,'Intensity'].to_numpy())

G9_neg_085_mean = np.mean(np.array(G9_neg_085_list), axis = 0)


G9_peak6 = glob.glob('*G9*6.7*22[56].csv')

G9_peak6_list = []
for i in range(len(G9_peak6)):
    G9_peak6_list.append(pd.read_csv(G9_peak6[i], sep=',').loc[:,'Intensity'].to_numpy())

G9_peak6_mean = np.mean(np.array(G9_peak6_list), axis = 0)

scaling_6 = 1313 / G9_neg_027_mean

#12.1
G9_peak12 = glob.glob('*G9*12.1*20[67].csv')
G9_epoch_12GHz = pd.read_csv(G9_peak12[0], sep=',').loc[:,'Epoch'].to_numpy()

G9_peak12_list = []
for i in range(len(G9_peak12)):
    G9_peak12_list.append(pd.read_csv(G9_peak12[i], sep=',').loc[:,'Intensity'].to_numpy())

G9_peak12_mean = np.mean(np.array(G9_peak12_list), axis = 0)


G9_neg_09 = glob.glob('*G9*12.1*16[012].csv')

G9_neg_09_list = []
for i in range(len(G9_neg_09)):
    G9_neg_09_list.append(pd.read_csv(G9_neg_09[i], sep=',').loc[:,'Intensity'].to_numpy())

G9_neg_09_mean = np.mean(np.array(G9_neg_09_list), axis = 0)

scaling_12 = 79.41 / G9_neg_09_mean


fig = go.Figure()
fig.add_trace(go.Scatter(x = G9_epoch_6GHz, y = G9_peak6_mean * scaling_6,
        mode='lines+markers',
        legendgroup="6",
        legendgrouptitle_text="Velocity (kms<sup>-1</sup>)",
        name='1.3 (Peak)',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array= (5 / 1313) * G9_peak6_mean * scaling_6,
            visible=True)
        ))
fig.add_hline(
    y=6550,
    line=dict(dash='dash', color='#3366CC'),
    name='1.3 (Peak) Pre-Flare Average'
)
fig.update_layout(
    title=str('G9.621+0.196 6.7 GHz 1.3 kms<sup>-1</sup> During Flare'),
    title_x=0.5,
    xaxis_title="Epoch (Day of 2021)",
    xaxis=dict(tickformat="%j"),
    yaxis_title="Flux Density (Jy)",
    yaxis_range=[0, 8000],
)
fig.show()


fig = go.Figure()
fig.add_trace(go.Scatter(x = G9_epoch_12GHz, y = G9_peak12_mean * scaling_12,
        mode='lines+markers',
        legendgroup="12",
        legendgrouptitle_text="12.1 GHz<br>Velocity (kms<sup>-1</sup>)",
        name='1.3 (Peak)',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=(5 / 83) * G9_peak12_mean * scaling_12,
            visible=True)
        ))
"""
fig.add_trace(go.Scatter(x = G9_epoch_12GHz, y = G9_neg_09_mean * scaling_12,
        mode='lines+markers',
        legendgroup="12",
        name='-0.9',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.001, 0],
            visible=True)
        ))
"""
fig.add_hline(
    y=284,
    line=dict(dash='dash', color='#3366CC'),
    name='1.3 (Peak) Pre-Flare Average'
)
fig.update_layout(
    title=str('G9.621+0.196 12.1 GHz 1.3 kms<sup>-1</sup> During Flare'),
    title_x=0.5,
    xaxis_title="Epoch (Day of 2021)",
    xaxis=dict(tickformat="%j"),
    yaxis_title="Flux Density (Jy)",
    yaxis_range=[0, 1000],
)
fig.show()



