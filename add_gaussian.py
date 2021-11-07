import argparse
import numpy as np
import pandas as pd

#Author: Callum Macdonald
#Direct any questions to either Callum.Macdonald@utas.edu.au or CallumRMacdonald@hotmail.com
#If another author contributes new methods to the script, please add your name and contact here

"""
A script deigned to add a gaussian profile to a csv file.
A csv file of the form produced by maser_plot.py (with date metadata in the first row) can
    be passed in along with the parameters for a gaussian
The script will save a new csv file that has that gaussian added to the original intensity data
"""

class add_gaussian:
    def __init__(self):
        
        #Use argparse to to accept the spectra and any customisation commands
        parser = argparse.ArgumentParser()

        #TODO: Change from nargs="+" to single or triple, not sure how to do it currently
            #Not a huge issue, just neater
        
        #The input spectrum
        #Single csv file, with date metadata in the first row, an intensity column and a velocity column
        parser.add_argument(
            "-spec",
            dest="spectrum",
            nargs="+"
        )

        #The Gaussian parameters
        #Three parameters, in order, amp, centre, stdv
        parser.add_argument(
            "-params",
            dest="parameters",
            nargs="+"
        )

        #Extract the input arguments
        args = parser.parse_args()
        self.spectrum = args.spectrum[0]
        self.parameters = args.parameters
    
    #Short function to create a gaussian profile
    #Velocity is the x values, params is a list of the three required parameters: amp, centre, stdv
    def gaussian(self, velocity, params):
        gaussian_profile = float(params[0]) * np.exp(-1 * (velocity - float(params[1]))**2 / (2 * float(params[2])**2))
        return gaussian_profile
    
    #Function to open the csv file, use gaussian to create the profile, add it to the intensity, and save the new csv
    def add_gaussian(self):

        #First opne the csv file, extracting the date information as well
        temp_file = open(self.spectrum, mode = 'r')
        epoch_string = temp_file.readline() #readline() reads the first line in a file

        csv_pandas_df = pd.read_csv(self.spectrum, sep=',', skiprows=1)
        velocity = csv_pandas_df.loc[:,'vel(km/s)'].to_numpy()
        intensity = csv_pandas_df.loc[:,'intensity(Jy)'].to_numpy()

        #Make a gaussian profile using self.gaussian
        gaussian_profile = self.gaussian(velocity= velocity, params = self.parameters)
        #Add the profile to the intensity array
        intensity_with_g = intensity + gaussian_profile

        #Now save everything to a new csv
        #Need to generate the name, by first dropping the .csv, then adding gaussian_{params}
        file_name_drop_csv = str(self.spectrum).replace('.csv','')

        file_name_new = file_name_drop_csv + '_gaussian_' + str(self.parameters[0]) + '_' + str(self.parameters[1]) + '_' + str(self.parameters[2]) + '.csv'

        #Create and open the csv file, then add the epoch data in the first row
        csv_file = open(str(file_name_new), "x")
        csv_file.write(epoch_string)

        #Put in the headers
        #Eah value needs to have quotation marks around it
        csv_file.write('"vel(km/s)"' + ',' + '"intensity(Jy)"' + '\n')

        #Write out the velocity and the new intensity data
        for i in range(len(velocity)):
            csv_file.write('"' + str(velocity[i]) + '"' + ',' + '"' + str(intensity_with_g[i]) + '"' + '\n')

        #Close the file, and everything is done
        csv_file.close()
    
adder = add_gaussian()
adder.add_gaussian()



