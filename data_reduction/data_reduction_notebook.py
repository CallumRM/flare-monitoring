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
os.chdir('Dropbox/Honours/Results/realtive_intensity/')


#First investigations:
#G9
#~~~~~~~~~~~~~~~~~~~~~~
G9_list = glob.glob('*G9*')

G9_df_list = []
for i in range(len(G9_list)):
    G9_df_list.append(pd.read_csv(G9_list[i], sep=','))
G9_name_list = np.array(G9_list, dtype=str)

fig = go.Figure()
for i in range(len(G9_df_list)):
    name_split = re.split('[.]', G9_name_list[i])[1]
    name = name_split[-10:]

    fig.add_trace(go.Scatter(x = G9_df_list[i].loc[:,'Epoch'].to_numpy(), y = G9_df_list[i].loc[:,'Intensity'].to_numpy(),
        mode='lines+markers',
        name=name
    ))
fig.update_layout(
    title=str('G9_6GHz'),
    yaxis_range=[0, 1]
    )

fig.show()

#~~~~~~~~~~~~~~~~~~~~~~

#G318
#~~~~~~~~~~~~~~~~~~~~~~
G318_list = glob.glob('*318*')

G318_df_list = []
for i in range(len(G318_list)):
    G318_df_list.append(pd.read_csv(G318_list[i], sep=','))
G318_name_list = np.array(G318_list, dtype=str)

fig = go.Figure()
for i in range(len(G318_df_list)):
    name_split = re.split('[.]', G318_name_list[i])[0]
    name = name_split[-10:]

    fig.add_trace(go.Scatter(x = G318_df_list[i].loc[:,'Epoch'].to_numpy(), y = G318_df_list[i].loc[:,'Intensity'].to_numpy(),
        mode='lines+markers',
        name=name
        ))
fig.update_layout(
    title=str('G318_6GHz'),
    yaxis_range=[0, 1]
    )
fig.show()

#G323
#~~~~~~~~~~~~~~~~~~~~~~
G323_list = glob.glob('*323*')

G323_df_list = []
for i in range(len(G323_list)):
    G323_df_list.append(pd.read_csv(G323_list[i], sep=','))
G323_name_list = np.array(G323_list, dtype=str)

fig = go.Figure()
for i in range(len(G323_df_list)):
    name_split = re.split('[.]', G323_name_list[i])[0]
    name = name_split[-10:]
    fig.add_trace(go.Scatter(x = G323_df_list[i].loc[:,'Epoch'].to_numpy(), y = G323_df_list[i].loc[:,'Intensity'].to_numpy(),
        mode='lines+markers',
        name=name
        ))
fig.update_layout(
    title=str('G323_6GHz'),
    yaxis_range=[0, 1]
    )
fig.show()

#G328.2
#~~~~~~~~~~~~~~~~~~~~~~
G328_2_list = glob.glob('*328.2*')

G328_2_df_list = []
for i in range(len(G328_2_list)):
    G328_2_df_list.append(pd.read_csv(G328_2_list[i], sep=','))
G328_2_name_list = np.array(G328_2_list, dtype=str)


fig = go.Figure()
for i in range(len(G328_2_df_list)):
    name_split = re.split('[.]', G328_2_name_list[i])[1]
    name = name_split[-10:]
    fig.add_trace(go.Scatter(x = G328_2_df_list[i].loc[:,'Epoch'].to_numpy(), y = G328_2_df_list[i].loc[:,'Intensity'].to_numpy(),
        mode='lines+markers',
        name=name
        ))
fig.update_layout(
    title=str('G328.2_6GHz'),
    yaxis_range=[0, 1]
    )
fig.show()


#G328.8
#~~~~~~~~~~~~~~~~~~~~~~
G328_8_list = glob.glob('*328.8*')

G328_8_df_list = []
for i in range(len(G328_8_list)):
    G328_8_df_list.append(pd.read_csv(G328_8_list[i], sep=','))
G328_8_name_list = np.array(G328_8_list, dtype=str)

fig = go.Figure()
for i in range(len(G328_8_df_list)):
    name_split = re.split('[.]', G328_8_name_list[i])[1]
    name = name_split[-10:]
    fig.add_trace(go.Scatter(x = G328_8_df_list[i].loc[:,'Epoch'].to_numpy(), y = G328_8_df_list[i].loc[:,'Intensity'].to_numpy(),
        mode='lines+markers',
        name=name
        ))
fig.update_layout(
    title=str('G328.8_6GHz'),
    yaxis_range=[0, 1]
    )
fig.show()

#G351
#~~~~~~~~~~~~~~~~~~~~~~
G351_list = glob.glob('*351*')

G351_df_list = []
for i in range(len(G351_list)):
    G351_df_list.append(pd.read_csv(G351_list[i], sep=','))
G351_name_list = np.array(G351_list, dtype=str)

fig = go.Figure()
for i in range(len(G351_df_list)):
    name_split = re.split('[.]', G351_name_list[i])[0]
    name = name_split[-10:]
    fig.add_trace(go.Scatter(x = G351_df_list[i].loc[:,'Epoch'].to_numpy(), y = G351_df_list[i].loc[:,'Intensity'].to_numpy(),
        mode='lines+markers',
        name=name
        ))
fig.update_layout(
    title=str('G351_6GHz'),
    yaxis_range=[0, 1],
    xaxis_range=['2021-04-25','2021-11-08']
    )
fig.show()



#Propper processing
#G9 demonstrate averaging
#~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~
G9_neg_027 = glob.glob('*G9*4[11-13]*')

G9_epoch_6GHz = pd.read_csv(G9_neg_027[0], sep=',').loc[:,'Epoch'].to_numpy()
G9_neg_027_list = []
for i in range(len(G9_neg_027)):
    G9_neg_027_list.append(pd.read_csv(G9_neg_027[i], sep=',').loc[:,'Intensity'].to_numpy())

G9_neg_027_mean = np.mean(np.array(G9_neg_027_list), axis = 0)
fig = go.Figure()
fig.add_trace(go.Scatter(x = G9_epoch_6GHz, y = G9_neg_027_list[2],
        mode='lines+markers',
        name='Peak Value',
        marker = {'color' : '#3366CC'},
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.001, 0],
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G9_epoch_6GHz, y = G9_neg_027_list[0],
        mode='lines+markers',
        name='+ Velocity Channel',
        line=dict(dash='dash'),
        marker = {'color' : '#109618'},
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.001, 0],
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G9_epoch_6GHz, y = G9_neg_027_list[1],
        mode='lines+markers',
        name='- Velocity Channel',
        line=dict(dash='dash'),
        marker = {'color' : '#FF9900'},
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.001, 0],
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G9_epoch_6GHz, y = G9_neg_027_mean,
        mode='lines+markers',
        name='Mean Value',
        marker = {'color' : '#DC3912'},
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.001, 0],
            visible=True)
        ))
fig.update_layout(
    title=str('Peak Averaging Example: G9.621+0.196 6.7GHz v=-0.27 km/s'),
    title_x=0.5,
    yaxis_range=[0, 1],
    xaxis_range=['2021-04-25','2021-11-08'],
    xaxis_title="Epoch (Day of 2021)",
    xaxis=dict(tickformat="%j"),
    yaxis_title="Relative Intensity",
)
fig.show()

#G9 full
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#6 Ghz
G9_neg_027 = glob.glob('*G9*4[11-13]*')

G9_epoch_6GHz = pd.read_csv(G9_neg_027[0], sep=',').loc[:,'Epoch'].to_numpy()
G9_neg_027_list = []
for i in range(len(G9_neg_027)):
    G9_neg_027_list.append(pd.read_csv(G9_neg_027[i], sep=',').loc[:,'Intensity'].to_numpy())

G9_neg_027_mean = np.mean(np.array(G9_neg_027_list), axis = 0)


G9_neg_085_1 = glob.glob('*G9*399*')
G9_neg_085_2 = glob.glob('*G9*400*')

G9_neg_085_list = [pd.read_csv(G9_neg_085_1[0], sep=',').loc[:,'Intensity'].to_numpy(), pd.read_csv(G9_neg_085_2[0], sep=',').loc[:,'Intensity'].to_numpy()]

G9_neg_085_mean = np.mean(np.array(G9_neg_085_list), axis = 0)



G9_54 = glob.glob('*G9*53[7-9]*')

G9_54_list = []
for i in range(len(G9_54)):
    G9_54_list.append(pd.read_csv(G9_54[i], sep=',').loc[:,'Intensity'].to_numpy())

G9_54_mean = np.mean(np.array(G9_54_list), axis = 0)

#12 GHz
G9_neg_09 = glob.glob('*G9*16[0-2]*')
G9_epoch_12GHz = pd.read_csv(G9_neg_09[0], sep=',').loc[:,'Epoch'].to_numpy()

G9_neg_09_list = []
for i in range(len(G9_neg_09)):
    G9_neg_09_list.append(pd.read_csv(G9_neg_09[i], sep=',').loc[:,'Intensity'].to_numpy())

G9_neg_09_mean = np.mean(np.array(G9_neg_09_list), axis = 0)

"""
#PLot
#6 GHz
fig = go.Figure()
fig.add_trace(go.Scatter(x = G9_epoch_6GHz, y = G9_54_mean,
        mode='lines+markers',
        name='5.4 km/s',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.001, 0],
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G9_epoch_6GHz, y = G9_neg_027_mean,
        mode='lines+markers',
        name='-0.27 km/s',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.001, 0],
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G9_epoch_6GHz, y = G9_neg_085_mean,
        mode='lines+markers',
        name='-0.85 km/s',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.001, 0],
            visible=True)
        ))
fig.update_layout(
    title=str('G9.621+0.196 6.7 GHz Relative Intensity'),
    title_x=0.5,
    xaxis_title="Epoch (Day of 2021)",
    xaxis=dict(tickformat="%j"),
    yaxis_title="Relative Intensity",
    yaxis_range=[0, 1],
    xaxis_range=['2021-04-25','2021-11-08']
)
fig.show()

#12 GHz
fig = go.Figure()
fig.add_trace(go.Scatter(x = G9_epoch_12GHz, y = G9_neg_09_mean,
        mode='lines+markers',
        name='-0.85 km/s',
        marker = {'color' : '#109618'},
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.03, 0],
            visible=True)
        ))
fig.update_layout(
    title=str('G9.621+0.196 12.1 GHz Relative Intensity'),
    title_x=0.5,
    xaxis_title="Epoch (Day of 2021)",
    xaxis=dict(tickformat="%j"),
    yaxis_title="Relative Intensity",
    yaxis_range=[0, 1],
    xaxis_range=['2021-04-25','2021-11-08']
)
fig.show()
"""
#Comparison
fig = go.Figure()
fig.add_trace(go.Scatter(x = G9_epoch_6GHz, y = G9_54_mean,
        mode='lines+markers',
        legendgroup="6",
        legendgrouptitle_text="6.7 GHz<br>Velocity (kms<sup>-1</sup>)",
        name='5.4',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.001, 0]*G9_54_mean[0:2]*2,
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G9_epoch_6GHz, y = G9_neg_027_mean,
        mode='lines+markers',
        legendgroup="6",
        name='-0.3',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.001, 0]*G9_neg_027_mean[0:2]*2,
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G9_epoch_6GHz, y = G9_neg_085_mean,
        mode='lines+markers',
        legendgroup="6",
        name='-0.9',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.001, 0]*G9_neg_085_mean[0:2]*2,
            visible=True)
        ))
#12 GHz
fig.add_trace(go.Scatter(x = G9_epoch_12GHz, y = G9_neg_09_mean,
        mode='lines+markers',
        legendgroup="12",
        legendgrouptitle_text="12.1 GHz<br>Velocity (kms<sup>-1</sup>)",
        name='-0.9',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.03, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.03, 0.03, 0]*G9_neg_09_mean*2,
            visible=True)
        ))
fig.update_layout(
    title=str('G9.621+0.196 6.7 and 12.1 GHz Relative Intensity'),
    title_x=0.5,
    xaxis_title="Epoch (Day of 2021)",
    xaxis=dict(tickformat="%j"),
    yaxis_title="Relative Intensity",
    yaxis_range=[0, 1],
    xaxis_range=['2021-04-25','2021-11-08']
)
fig.show()

#318 full
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#6GHz
G318_neg_38_2 = glob.glob('*G318*36[67]*')

G318_epoch_6GHz = pd.read_csv(G318_neg_38_2[0], sep=',').loc[:,'Epoch'].to_numpy()

G318_neg_38_2_list = []
for i in range(len(G318_neg_38_2)):
    G318_neg_38_2_list.append(pd.read_csv(G318_neg_38_2[i], sep=',').loc[:,'Intensity'].to_numpy())

G318_neg_38_2_average = np.mean(np.array(G318_neg_38_2_list), axis = 0)


G318_neg_36_2 = glob.glob('*G318*41[12]*')
G318_neg_36_2_list = []
for i in range(len(G318_neg_36_2)):
    G318_neg_36_2_list.append(pd.read_csv(G318_neg_36_2[i], sep=',').loc[:,'Intensity'].to_numpy())

G318_neg_36_2_average = np.mean(np.array(G318_neg_36_2_list), axis = 0)


G318_neg_37 = glob.glob('*G318*39*')
G318_neg_37_list = []
for i in range(len(G318_neg_37)):
    G318_neg_37_list.append(pd.read_csv(G318_neg_37[i], sep=',').loc[:,'Intensity'].to_numpy())

G318_neg_37_average = np.mean(np.array(G318_neg_37_list), axis = 0)

#PLot
#6 GHz
fig = go.Figure()
fig.add_trace(go.Scatter(x = G318_epoch_6GHz, y = G318_neg_36_2_average,
        mode='lines+markers',
        name='-36.2',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.014, 0]*G318_neg_36_2_average[0:2]*2,
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G318_epoch_6GHz, y = G318_neg_37_average,
        mode='lines+markers',
        name='-37.0',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.014, 0]*G318_neg_37_average[0:2]*2,
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G318_epoch_6GHz, y = G318_neg_38_2_average,
        mode='lines+markers',
        name='-38.2',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.014, 0]*G318_neg_38_2_average[0:2]*2,
            visible=True)
        ))
fig.update_layout(
    title=str('G318.948-0.196 6.7 GHz Relative Intensity'),
    title_x=0.5,
    xaxis_title="Epoch (Day of 2021)",
    xaxis=dict(tickformat="%j"),
    yaxis_title="Relative Intensity",
    yaxis_range=[0, 1],
    xaxis_range=['2021-04-25','2021-11-08'],
    legend_title_text="6.7 GHz<br>Velocity (kms<sup>-1</sup>)",
)
fig.show()

#323 full
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
G323_neg_49_27 = glob.glob('*G323*475*')

G323_epoch_6GHz = pd.read_csv(G323_neg_49_27[0], sep=',').loc[:,'Epoch'].to_numpy()
G323_neg_49_27_array = pd.read_csv(G323_neg_49_27[0], sep=',').loc[:,'Intensity'].to_numpy()

G323_neg_48_65 = glob.glob('*G323*48[89]*')
G323_neg_48_65_list = []
for i in range(len(G323_neg_48_65)):
    G323_neg_48_65_list.append(pd.read_csv(G323_neg_48_65[i], sep=',').loc[:,'Intensity'].to_numpy())

G323_neg_48_65_average = np.mean(np.array(G323_neg_48_65_list), axis = 0)

G323_neg_51 = glob.glob('*G323*43[2-4]*')
G323_neg_51_list = []
for i in range(len(G323_neg_51)):
    G323_neg_51_list.append(pd.read_csv(G323_neg_51[i], sep=',').loc[:,'Intensity'].to_numpy())

G323_neg_51_average = np.mean(np.array(G323_neg_51_list), axis = 0)

#12 GHz
G323_neg_48_6 = glob.glob('*G323*24[2-4]*')
G323_epoch_12GHz = pd.read_csv(G323_neg_48_6[0], sep=',').loc[:,'Epoch'].to_numpy()

G323_neg_48_6_list = []
for i in range(len(G323_neg_48_6)):
    G323_neg_48_6_list.append(pd.read_csv(G323_neg_48_6[i], sep=',').loc[:,'Intensity'].to_numpy())

G323_neg_48_6_average = np.mean(np.array(G323_neg_48_6_list), axis = 0)


G323_neg_50_85 = glob.glob('*G323*19[7-9]*')
G323_neg_50_85_list = []
for i in range(len(G323_neg_50_85)):
    G323_neg_50_85_list.append(pd.read_csv(G323_neg_50_85[i], sep=',').loc[:,'Intensity'].to_numpy())

G323_neg_50_85_average = np.mean(np.array(G323_neg_50_85_list), axis = 0)


G323_neg_51_6 = glob.glob('*G323*18[23]*')
G323_neg_51_6_list = []
for i in range(len(G323_neg_51_6)):
    G323_neg_51_6_list.append(pd.read_csv(G323_neg_51_6[i], sep=',').loc[:,'Intensity'].to_numpy())

G323_neg_51_6_average = np.mean(np.array(G323_neg_51_6_list), axis = 0)
#6
fig = go.Figure()
fig.add_trace(go.Scatter(x = G323_epoch_6GHz, y = G323_neg_48_65_average,
        mode='lines+markers',
        legendgroup="6",
        legendgrouptitle_text="6.7 GHz<br>Velocity (kms<sup>-1</sup>)",
        name='-48.5',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.003, 0]*G323_neg_48_65_average[0:2]*2,
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G323_epoch_6GHz, y = G323_neg_49_27_array,
        mode='lines+markers',
        legendgroup="6",
        name='-49.3',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.003, 0]*G323_neg_49_27_array[0:2]*2,
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G323_epoch_6GHz, y = G323_neg_51_average,
        mode='lines+markers',
        legendgroup="6",
        name='-51.2',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.003, 0]*G323_neg_51_average[0:2]*2,
            visible=True)
        ))

#12
fig.add_trace(go.Scatter(x = G323_epoch_12GHz, y = G323_neg_48_6_average,
        mode='lines+markers',
        legendgroup="12",
        legendgrouptitle_text="12.1 GHz<br>Velocity (kms<sup>-1</sup>)",
        name='-48.2',
        marker = {'color' : '#109618'},
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.025, 0]*G323_neg_48_6_average[0:2]*2,
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G323_epoch_12GHz, y = G323_neg_50_85_average,
        mode='lines+markers',
        legendgroup="12",
        name='-50.5',
        marker = {'color' : '#990099'},
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0, 0.025, 0]*G323_neg_50_85_average[0:3]*2,
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G323_epoch_12GHz, y = G323_neg_51_6_average,
        mode='lines+markers',
        legendgroup="12",
        name='-51.2',
        marker = {'color' : '#0099C6'},
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.025, 0]*G323_neg_51_6_average[0:2]*2,
            visible=True)
        ))
fig.update_layout(
    title=str('G323.740-0.263 6.7 and 12.1 GHz Relative Intensity'),
    title_x=0.5,
    xaxis_title="Epoch (Day of 2021)",
    xaxis=dict(tickformat="%j"),
    yaxis_title="Relative Intensity",
    yaxis_range=[0, 1],
    xaxis_range=['2021-04-25','2021-11-08']
)
fig.show()

"""
#Comparison shifted
#6
fig = go.Figure()
fig.add_trace(go.Scatter(x = G323_epoch_6GHz, y = G323_neg_48_65_average-0.23,
        mode='lines+markers',
        legendgroup="6",
        legendgrouptitle_text="6.7 GHz<br>Velocity (kms<sup>-1</sup>)",
        name='-48.5',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.003, 0],
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G323_epoch_6GHz, y = G323_neg_49_27_array-0.23, #12 gig norm on this one
        mode='lines+markers',
        legendgroup="6",
        name='-49.3',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.003, 0],
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G323_epoch_6GHz, y = G323_neg_51_average-0.23,
        mode='lines+markers',
        legendgroup="6",
        name='-51.2',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.003, 0],
            visible=True)
        ))
#12
fig.add_trace(go.Scatter(x = G323_epoch_12GHz, y = G323_neg_48_6_average-0.38-0.23,
        mode='lines+markers',
        legendgroup="12",
        legendgrouptitle_text="12.1 GHz<br>Velocity (kms<sup>-1</sup>)",
        name='-48.2',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.025, 0],
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G323_epoch_12GHz, y = G323_neg_50_85_average-0.31-0.23, #6 Norm this one
        mode='lines+markers',
        legendgroup="12",
        name='-50.5',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.025, 0],
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G323_epoch_12GHz, y = G323_neg_51_6_average+0.13-0.23,
        mode='lines+markers',
        legendgroup="12",
        name='-51.2',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0, 0.025, 0],
            visible=True)
        ))

fig.update_layout(
    title=str('G323.740-0.263 12.1 GHz and 6.7 GHz Relative Intensity (Values Shifted For Comparison)'),
    title_x=0.5,
    xaxis_title="Epoch (Day of 2021)",
    xaxis=dict(tickformat="%j"),
    yaxis_title="Relative Intensity",
    yaxis_range=[0, 1],
    xaxis_range=['2021-04-25','2021-11-08'],
)
fig.show()
"""

#328.2 full
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#6GHz
G328_2_neg_37_5 = glob.glob('*G328.2*608*')

G328_2_epoch_6GHz = pd.read_csv(G328_2_neg_37_5[0], sep=',').loc[:,'Epoch'].to_numpy()

G328_2_neg_37_5_list = []
for i in range(len(G328_2_neg_37_5)):
    G328_2_neg_37_5_list.append(pd.read_csv(G328_2_neg_37_5[i], sep=',').loc[:,'Intensity'].to_numpy())

G328_2_neg_37_5_average = np.mean(np.array(G328_2_neg_37_5_list), axis = 0)


G328_2_neg_43_3 = glob.glob('*G328.2*4[78]*')
G328_2_neg_43_3_list = []
for i in range(len(G328_2_neg_43_3)):
    G328_2_neg_43_3_list.append(pd.read_csv(G328_2_neg_43_3[i], sep=',').loc[:,'Intensity'].to_numpy())

G328_2_neg_43_3_average = np.mean(np.array(G328_2_neg_43_3_list), axis = 0)


G328_2_neg_45_4 = glob.glob('*G328.2*43[123]*')
G328_2_neg_45_4_list = []
for i in range(len(G328_2_neg_45_4)):
    G328_2_neg_45_4_list.append(pd.read_csv(G328_2_neg_45_4[i], sep=',').loc[:,'Intensity'].to_numpy())

G328_2_neg_45_4_average = np.mean(np.array(G328_2_neg_45_4_list), axis = 0)


G328_2_neg_32_5 = glob.glob('*G328.2*7[12]*')
G328_2_neg_32_5_list = []
for i in range(len(G328_2_neg_32_5)):
    G328_2_neg_32_5_list.append(pd.read_csv(G328_2_neg_32_5[i], sep=',').loc[:,'Intensity'].to_numpy())

G328_2_neg_32_5_average = np.mean(np.array(G328_2_neg_32_5_list), axis = 0)


G328_2_neg_48_6 = glob.glob('*G328.2*36[012]*')
G328_2_neg_48_6_list = []
for i in range(len(G328_2_neg_48_6)):
    G328_2_neg_48_6_list.append(pd.read_csv(G328_2_neg_48_6[i], sep=',').loc[:,'Intensity'].to_numpy())

G328_2_neg_48_6_average = np.mean(np.array(G328_2_neg_48_6_list), axis = 0)


G328_2_neg_49_5 = glob.glob('*G328.2*34[12]*')
G328_2_neg_49_5_list = []
for i in range(len(G328_2_neg_49_5)):
    G328_2_neg_49_5_list.append(pd.read_csv(G328_2_neg_49_5[i], sep=',').loc[:,'Intensity'].to_numpy())

G328_2_neg_49_5_average = np.mean(np.array(G328_2_neg_49_5_list), axis = 0)

#PLot
#6 GHz
fig = go.Figure()
fig.add_trace(go.Scatter(x = G328_2_epoch_6GHz, y = G328_2_neg_37_5_average,
        mode='lines+markers',
        name='-37.5',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.01, 0]*G328_2_neg_37_5_average[0:2]*2,
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G328_2_epoch_6GHz, y = G328_2_neg_43_3_average,
        mode='lines+markers',
        name='-43.3',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.01, 0]*G328_2_neg_43_3_average[0:2]*2,
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G328_2_epoch_6GHz, y = G328_2_neg_45_4_average,
        mode='lines+markers',
        name='-45.4',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.01, 0]*G328_2_neg_45_4_average[0:2]*2,
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G328_2_epoch_6GHz, y = G328_2_neg_48_6_average,
        mode='lines+markers',
        name='-48.6',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.01, 0]*G328_2_neg_48_6_average[0:2]*2,
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G328_2_epoch_6GHz, y = G328_2_neg_49_5_average,
        mode='lines+markers',
        name='-49.5',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.01, 0]*G328_2_neg_49_5_average[0:2]*2,
            visible=True)
        ))
fig.update_layout(
    title=str('G328.237-0.547 6.7 GHz Relative Intensity'),
    title_x=0.5,
    xaxis_title="Epoch (Day of 2021)",
    xaxis=dict(tickformat="%j"),
    yaxis_title="Relative Intensity",
    yaxis_range=[0, 1],
    xaxis_range=['2021-04-25','2021-11-08'],
    legend_title_text="6.7 GHz<br>Velocity (kms<sup>-1</sup>)",
)
fig.show()



#328.8 full
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#6GHz
G328_8_neg_43_8 = glob.glob('*G328.8*46[123]*')

G328_8_epoch_6GHz = pd.read_csv(G328_8_neg_43_8[0], sep=',').loc[:,'Epoch'].to_numpy()

G328_8_neg_43_8_list = []
for i in range(len(G328_8_neg_43_8)):
    G328_8_neg_43_8_list.append(pd.read_csv(G328_8_neg_43_8[i], sep=',').loc[:,'Intensity'].to_numpy())

G328_8_neg_43_8_average = np.mean(np.array(G328_8_neg_43_8_list), axis = 0)


G328_8_neg_45_2 = glob.glob('*G328.8*43[012]*')
G328_8_neg_45_2_list = []
for i in range(len(G328_8_neg_45_2)):
    G328_8_neg_45_2_list.append(pd.read_csv(G328_8_neg_45_2[i], sep=',').loc[:,'Intensity'].to_numpy())

G328_8_neg_45_2_average = np.mean(np.array(G328_8_neg_45_2_list), axis = 0)


G328_8_neg_46_3 = glob.glob('*G328.8*40[789]*')
G328_8_neg_46_3_list = []
for i in range(len(G328_8_neg_46_3)):
    G328_8_neg_46_3_list.append(pd.read_csv(G328_8_neg_46_3[i], sep=',').loc[:,'Intensity'].to_numpy())

G328_8_neg_46_3_average = np.mean(np.array(G328_8_neg_46_3_list), axis = 0)


G328_8_neg_46_6_1 = glob.glob('*G328.8*399*')
G328_8_neg_46_6_2 = glob.glob('*G328.8*40[01]*')
G328_8_neg_46_6_2.append(G328_8_neg_46_6_1[0])
G328_8_neg_46_6_list = []
for i in range(len(G328_8_neg_46_6_2)):
    G328_8_neg_46_6_list.append(pd.read_csv(G328_8_neg_46_6_2[i], sep=',').loc[:,'Intensity'].to_numpy())

G328_8_neg_46_6_average = np.mean(np.array(G328_8_neg_46_6_list), axis = 0)


#PLot
#6 GHz
fig = go.Figure()
fig.add_trace(go.Scatter(x = G328_8_epoch_6GHz, y = G328_8_neg_43_8_average,
        mode='lines+markers',
        name='-43.8',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.05, 0.05, 0.025, 0]*G328_8_neg_43_8_average[0:4]*2,
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G328_8_epoch_6GHz, y = G328_8_neg_45_2_average,
        mode='lines+markers',
        name='-45.2',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.05, 0.05, 0.025, 0]*G328_8_neg_45_2_average[0:4]*2,
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G328_8_epoch_6GHz, y = G328_8_neg_46_3_average,
        mode='lines+markers',
        name='-46.2',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.05, 0.05, 0.025, 0]*G328_8_neg_46_3_average[0:4]*2,
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G328_8_epoch_6GHz, y = G328_8_neg_46_6_average,
        mode='lines+markers',
        name='-46.6',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.05, 0.05, 0.025, 0]*G328_8_neg_46_6_average[0:4]*2,
            visible=True)
        ))
fig.update_layout(
    title=str('G328.809+0.633 6.7 GHz Relative Intensity'),
    title_x=0.5,
    xaxis_title="Epoch (Day of 2021)",
    xaxis=dict(tickformat="%j"),
    yaxis_title="Relative Intensity",
    yaxis_range=[0, 1],
    xaxis_range=['2021-04-25','2021-11-08'],
    legend_title_text="6.7 GHz<br>Velocity (kms<sup>-1</sup>)",
)
fig.show()



#351 full
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#6GHz
G351_neg_8_6_1 = glob.glob('*G351*48[789]*')
G351_neg_8_6_2 = glob.glob('*G328.8*490*')
G351_neg_8_6_1.append(G351_neg_8_6_2[0])
G351_neg_8_6_list = []
for i in range(len(G351_neg_8_6_1)):
    G351_neg_8_6_list.append(pd.read_csv(G351_neg_8_6_1[i], sep=',').loc[:,'Intensity'].to_numpy())

G351_neg_8_6_average = np.mean(np.array(G351_neg_8_6_list), axis = 0)

G351_epoch_6GHz = pd.read_csv(G351_neg_8_6_1[0], sep=',').loc[:,'Epoch'].to_numpy()


G351_neg_8_9_1 = glob.glob('*G351*48[012]*')
G351_neg_8_9_2 = glob.glob('*G328.8*479*')
G351_neg_8_9_1.append(G351_neg_8_9_2[0])
G351_neg_8_9_list = []
for i in range(len(G351_neg_8_9_1)):
    G351_neg_8_9_list.append(pd.read_csv(G351_neg_8_9_1[i], sep=',').loc[:,'Intensity'].to_numpy())

G351_neg_8_9_average = np.mean(np.array(G351_neg_8_9_list), axis = 0)


G351_neg_9_4 = glob.glob('*G351*47[0123]*')
G351_neg_9_4_list = []
for i in range(len(G351_neg_9_4)):
    G351_neg_9_4_list.append(pd.read_csv(G351_neg_9_4[i], sep=',').loc[:,'Intensity'].to_numpy())

G351_neg_9_4_average = np.mean(np.array(G351_neg_9_4_list), axis = 0)


G351_neg_6_6 = glob.glob('*G351*53[234]*')
G351_neg_6_6_list = []
for i in range(len(G351_neg_6_6)):
    G351_neg_6_6_list.append(pd.read_csv(G351_neg_6_6[i], sep=',').loc[:,'Intensity'].to_numpy())

G351_neg_6_6_average = np.mean(np.array(G351_neg_6_6_list), axis = 0)


G351_neg_7_4 = glob.glob('*G351*51[456]*')
G351_neg_7_4_list = []
for i in range(len(G351_neg_7_4)):
    G351_neg_7_4_list.append(pd.read_csv(G351_neg_7_4[i], sep=',').loc[:,'Intensity'].to_numpy())

G351_neg_7_4_average = np.mean(np.array(G351_neg_7_4_list), axis = 0)


G351_neg_9_7 = glob.glob('*G351*46[23]*')
G351_neg_9_7_list = []
for i in range(len(G351_neg_9_7)):
    G351_neg_9_7_list.append(pd.read_csv(G351_neg_9_7[i], sep=',').loc[:,'Intensity'].to_numpy())

G351_neg_9_7_average = np.mean(np.array(G351_neg_9_7_list), axis = 0)


G351_neg_11_2 = glob.glob('*G351*43[012]*')
G351_neg_11_2_list = []
for i in range(len(G351_neg_11_2)):
    G351_neg_11_2_list.append(pd.read_csv(G351_neg_11_2[i], sep=',').loc[:,'Intensity'].to_numpy())

G351_neg_11_2_average = np.mean(np.array(G351_neg_11_2_list), axis = 0)


#12 GHz
G351_neg_9_6 = glob.glob('*G351*22[012]*')
G351_neg_9_6_list = []
for i in range(len(G351_neg_9_6)):
    G351_neg_9_6_list.append(pd.read_csv(G351_neg_9_6[i], sep=',').loc[:,'Intensity'].to_numpy())

G351_neg_9_6_average = np.mean(np.array(G351_neg_9_6_list), axis = 0)

G351_epoch_12GHz = pd.read_csv(G351_neg_9_6[0], sep=',').loc[:,'Epoch'].to_numpy()


G351_neg_10_9 = glob.glob('*G351*19[56]*')
G351_neg_10_9_list = []
for i in range(len(G351_neg_10_9)):
    G351_neg_10_9_list.append(pd.read_csv(G351_neg_10_9[i], sep=',').loc[:,'Intensity'].to_numpy())

G351_neg_10_9_average = np.mean(np.array(G351_neg_10_9_list), axis = 0)


G351_neg_11_1 = glob.glob('*G351*1[89][901]*')
G351_neg_11_1_list = []
for i in range(len(G351_neg_11_1)):
    G351_neg_11_1_list.append(pd.read_csv(G351_neg_11_1[i], sep=',').loc[:,'Intensity'].to_numpy())

G351_neg_11_1_average = np.mean(np.array(G351_neg_11_1_list), axis = 0)


#PLot
#6 GHz
fig = go.Figure()
fig.add_trace(go.Scatter(x = G351_epoch_6GHz, y = G351_neg_6_6_average,
        mode='lines+markers',
        legendgroup="6",
        legendgrouptitle_text="6.7 GHz<br>Velocity (kms<sup>-1</sup>)",
        name='-6.6',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.002, 0]*G351_neg_6_6_average[0:2]*2,
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G351_epoch_6GHz, y = G351_neg_7_4_average,
        mode='lines+markers',
        legendgroup="6",
        name='-7.4',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.002, 0]*G351_neg_7_4_average[0:2]*2,
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G351_epoch_6GHz, y = G351_neg_8_6_average,
        mode='lines+markers',
        legendgroup="6",
        name='-8.6',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.002, 0]*G351_neg_8_6_average[0:2]*2,
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G351_epoch_6GHz, y = G351_neg_8_9_average,
        mode='lines+markers',
        legendgroup="6",
        name='-8.9',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.002, 0]*G351_neg_8_9_average[0:2]*2,
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G351_epoch_6GHz, y = G351_neg_9_4_average,
        mode='lines+markers',
        legendgroup="6",
        name='-9.4',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.002, 0]*G351_neg_9_4_average[0:2]*2,
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G351_epoch_6GHz, y = G351_neg_9_7_average,
        mode='lines+markers',
        legendgroup="6",
        name='-9.7',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.002, 0]*G351_neg_9_7_average[0:2]*2,
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G351_epoch_6GHz, y = G351_neg_11_2_average,
        mode='lines+markers',
        legendgroup="6",
        name='-11.2',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.002, 0]*G351_neg_11_2_average[0:2]*2,
            visible=True)
        ))

#12 GHz
fig.add_trace(go.Scatter(x = G351_epoch_12GHz, y = G351_neg_9_6_average,
        mode='lines+markers',
        legendgroup="12",
        legendgrouptitle_text="12.1 GHz<br>Velocity (kms<sup>-1</sup>)",
        name='-9.7',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.01, 0]*G351_neg_9_6_average[0:2]*2,
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G351_epoch_12GHz, y = G351_neg_10_9_average,
        mode='lines+markers',
        legendgroup="12",
        name='-10.8',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.01, 0]*G351_neg_10_9_average[0:2]*2,
            visible=True)
        ))
fig.add_trace(go.Scatter(x = G351_epoch_12GHz, y = G351_neg_11_1_average,
        mode='lines+markers',
        legendgroup="12",
        name='-11.2',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=[0.01, 0]*G351_neg_11_1_average[0:2]*2,
            visible=True)
        ))
fig.update_layout(
    title=str('G351.417+0.645 6.7 and 12.1 GHz Relative Intensity'),
    title_x=0.5,
    xaxis_title="Epoch (Day of 2021)",
    xaxis=dict(tickformat="%j"),
    yaxis_title="Relative Intensity",
    yaxis_range=[0, 1],
    xaxis_range=['2021-04-25','2021-11-08']
)
fig.show()