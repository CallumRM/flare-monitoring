import argparse
import numpy as np
import pandas as pd

#Author: Callum Macdonald
#Direct any questions to either Callum.Macdonald@utas.edu.au or CallumRMacdonald@hotmail.com
#If another author contributes new methods to the script, please add your name and contact here

"""
A script to extract a single channel from a list spectra csvs.
If a list of channel numbers is entered, a csv will be created for each channel
Using the maser_plot.py csv format with the date metadata in the first row, will create a new csv file with
    intensity and date as the columns.
"""
#TODO: I wrote this intiial to just extarct a single channel, then realised it would be way more useful if it could
    #extract multiple. As such it is terribly optimised, and could be made to run much faster with a bit of rearranging
    #Only takes of order 2 seconds to run, but still.
#TODO: Ambitious, make a script that will automatically extract and plot. Currently doing plotting manually,
    #but for a large number of sources have a list of velocities and an ability to just run a script on a
    #directory and have a bunch of plots produced would be good.

class channel_extract:
    def __init__(self):
        
        #Use argparse to to accept the spectra and any customisation commands
        parser = argparse.ArgumentParser()

        #The input spectra
        #Should be a list of csv files, each with date metadata in the first row 
        parser.add_argument(
            "-spec",
            dest="spectra",
            nargs="+"
        )

        #The name to be used for the output csv file.
        #Will have channel number added to it automatically
        #Sensible name would be source_frequency
        parser.add_argument(
            "-name",
            dest="name",
            type=str,
            default='no_name'
        )

        #The channel numbers to extract
        #A list of integers, that can be determine beforehand by inspecting the csv
        #TODO: Quality of life, add an option to find the channel nearest a given velocity, then extract that
        parser.add_argument(
            "-channels",
            dest="channels",
            nargs="+",
        )

        #Channel number correction
        #The channel number specified above corresponds to the position within the python array, rather 
            #than the line number I would see looking at a csv in a notepad
        #The difference is 2 (for metadata and headers) + 1 (python starts at 0, notepad at 1)
        #This correction will add subtract 3 from the channel number specified, so the number read
            #in a csv can be entered
        parser.add_argument(
            "-chan_cor",
            dest="channel_correction",
            type=str,
            default='True'
        )

        #Extract the input information
        args = parser.parse_args()
        self.spectra = args.spectra
        self.name = str(args.name)
        self.channels = np.array(args.channels, dtype=int)

        #If the channel correction is to be applied, subtract 3 from the channel number
        if str(args.channel_correction) == 'True':
            self.channels = self.channels - 3

    #Method to extract the date and intensity from a single csv file
    #Will return a tuple with the intensity and epoch
    def extract_one_scan(self, csv_file, channel):

        #Open the csv file using read, to extract the date metadata
        epoch_csv_file = open(csv_file, mode = 'r')
        epoch_string = epoch_csv_file.readline() #readline() reads the first line in a file
        epoch = np.datetime64(epoch_string[1:]) #Use string[1:] to drop the first char, which is the # used to indicate metadata

        #Open the data of the csv file using pandas. No need for velocity, unless that is implemented later
        csv_pandas_df = pd.read_csv(csv_file, sep=',', skiprows=1)
        intensity = csv_pandas_df.loc[:,'intensity(Jy)'].to_numpy()
        velocity = csv_pandas_df.loc[:,'vel(km/s)'].to_numpy()

        #Extract the single intensity value
        intensity_value = intensity[channel]

        #Return the value and epoch as a tuple
        return (intensity_value, epoch)
    
    #Method to write a csv file
    #Intensities will be the y_values, epochs the x_values
    def write_csv(self, name, intensities, epochs, x_name, y_name):

        #Create and open the csv file
        data = open(str(name) + ".csv", "x")

        #Put in the headers
        #Eah value needs to have quotation marks around it
        data.write('"' + str(x_name) + '"' + ',' + '"' + str(y_name) + '"' + '\n')

        #Write out the dates with the intensity on that date
        for i in range(len(epochs)):
            data.write('"' + str(epochs[i]) + '"' + ',' + '"' + str(intensities[i]) + '"' + '\n')

        data.close()
    
    #Method to extract the a single channel
    #Will loop through each scan, and use extract_one_scan to get the information out
    #Will then use write_csv to write the output
    def extract_one_channel(self, channel):
        
        print("Now extracting channel: " + str(channel + 3))
        #Lists to hold the epochs and intensities
        epoch_list = []
        intensity_list = []

        #Loop through every value in the spectra list, and extract the epoch and intensity
        for i in range(len(self.spectra)):

            output = self.extract_one_scan(self.spectra[i], channel)
            
            #Output has the form (intensity, epoch)
            intensity_list.append(output[0])
            epoch_list.append(output[1])
        
        #Now create the full name and write the csv file
        name = 'single_channel_' + str(self.name) + '_ch' + str(channel+3) 

        self.write_csv(name = name, intensities = intensity_list, epochs = epoch_list, x_name = 'Epoch', y_name = 'Intensity')
    
    #Method to extract each channel entered
    #Will run extract_one_channel for each channel
    def extract_all(self):

        #Loop through each channel
        for i in range(len(self.channels)):
            self.extract_one_channel(channel = self.channels[i])
    

extractor = channel_extract()
extractor.extract_all()
            


            
