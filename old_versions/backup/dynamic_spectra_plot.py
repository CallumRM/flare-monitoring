import numpy as np
import pandas as pd
import plotly.graph_objs as go
import argparse

#Author: Callum Macdonald
#Direct any questions to either Callum.Macdonald@utas.edu.au or CallumRMacdonald@hotmail.com
#If another author contributes new methods to the script, please add your name and contact here

"""
A script that will take several csv input spectra and produce a dynamic time spectra using plotly's heatmap
The csvs will need to have the date in the metadata at the start of the file (maser_plot does this automatically fo csv files)
The velocity axis will be taken from the first spectra in the input list

"""
#TODO: Check if all of the velocities passed in are the same, and print out a warning if they are not
#TODO: Add a feature to reduce the plotting velocity width, or just specify a window
class dynamic_spectra:
    def __init__(self):
        
        #Use argparse to to accept the spectra and any customisation commands
        parser = argparse.ArgumentParser()

        #The input spectra
        #Should be a list of csv files, each with 
        parser.add_argument(
            "-spec",
            dest="spectra",
            nargs="+"
        )

        #Max z value to plot
        parser.add_argument(
            "-maxz",
            dest="max_z",
            type=float
        )

        #Min z value to plot
        parser.add_argument(
            "-minz",
            dest="min_z",
            type=float
        )

        #Log plot
        #Trueth string used to determine if a log plot should be generated instead of a linear plot
        #Will use log base 10
        parser.add_argument(
            "-log",
            dest="log_plot",
            type=str,
            default='False'
        )

        #Title
        #The title of the plot
        parser.add_argument(
            "-title",
            dest="title",
            type=str,
            default='Title'
        )

        #Extract the input arguments
        args = parser.parse_args()

        self.spectra = args.spectra
        self.max_z = args.max_z
        self.min_z = args.min_z
        self.log_plot = str(args.log_plot)
        self.title = str(args.title)

        #Extract the information from the spectra
        self.epochs = []
        self.velocities = []
        self.intensities = []

        for i in range(len(self.spectra)):
            #Get the intensity and velicity data out for each spectra. Skiprows used to avoid the epoch data
            csv_pandas_df = pd.read_csv(self.spectra[i], sep=',', skiprows=1)
            self.velocities.append(csv_pandas_df.loc[:,'vel(km/s)'].to_numpy())

            #If a logt plot is to be produced, take the log of the abslute intensity
            if self.log_plot == 'True':
                #Need the absolute intensity, so that log doesn't return nan
                abs_intensity = np.abs(csv_pandas_df.loc[:,'intensity(Jy)'].to_numpy())
                log_10_intensity = np.log10(abs_intensity)
                self.intensities.append(log_10_intensity)
            else:
                self.intensities.append(csv_pandas_df.loc[:,'intensity(Jy)'].to_numpy())

            #Get the epoch data out
            temp_file = open(self.spectra[i], mode = 'r')
            epoch_temp_string = temp_file.readline() #readline() reads the first line in a file
            epoch_temp = np.datetime64(epoch_temp_string[1:]) #Use string[1:] to drop the first char, which is the # used to indicate metadata
            self.epochs.append(epoch_temp)

    def plot_spectra(self):

        #To plot the spectra, need to just use one velocity. Just use the one from the first file
        velocity = self.velocities[0]
        figure = go.Figure(data = go.Heatmap(
                    z = self.intensities,
                    x = velocity,
                    y = self.epochs,
                    ygap=0.5,
                    #zsmooth = 'best',
                    colorscale='Jet'
                    ))
        
        figure.update_layout(
            title=str(self.title)
        )

        #If a minimum and maximum value for z have been passed in, then use them to scale the plot
        if (self.min_z != None) and (self.max_z != None):
            figure.data[0].update(zmin=self.min_z, zmax=self.max_z)

        figure.show()

dynspec = dynamic_spectra()
dynspec.plot_spectra()

