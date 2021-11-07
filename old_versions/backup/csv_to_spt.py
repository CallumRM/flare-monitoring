from astropy.coordinates import spectral_coordinate
from astropy.units.equivalencies import spectral
import numpy as np
import os, argparse, re
#import astropy.time
import pandas as pd


#Will take in spectra from MMB in .dat form and convert them to .spt for easy compirosn of data.
#Will give the new .spt files the same name as the input file.
#Code largly copyed from my hart_to_spt script.

class csv_to_spt:
    def __init__(self):

        #Use argparse to accept the input arguments
        #To pass an input arugment, enter the flag followed by the value.
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "-spec",
            dest="spec_files",
            help="The input .csv files",
            nargs="+" 
            #List of spectra
        )

        #Velocity range
        #The range of velocities to generate the plot for. First value is lower, second is upper limit
        #These will be centred on the source velocity.
        #eg entering [-20, 30] for a source centred at 55 will give a range of 35 to 85
        #TODO: Not implemented
        """
        parser.add_argument(
            "-vr",
            dest="velocity_range",
            nargs="+",
            default=[-20,20] #defualt set to standard maser zoom
        )
        """

        args = parser.parse_args()
        self.spec_files = args.spec_files
        self.termination_message = 'Script terminating'

        if self.spec_files == None:
            print('No spectra file passed.\n' + self.termination_message)
            quit()
        
        #Extract the spectra from the files using pandas read_csv
        #Create empty lists for storing the names, velocity and intensity data for each file passed in.
        self.name_list = []
        self.velocities = []
        self.intensities = []
        for i in range(len(self.spec_files)):

            #Get the name of the file. As the file might be in another directory, split at '/'
            full_name_split = re.split('[/]', str(self.spec_files[i]))

            #Extract the list string (the name) and remove the .csv at the end
            temp_name = str(full_name_split[-1]).replace('.csv','')
            self.name_list.append(temp_name)

            #Open the csv file and get the velocity and intensity information out
            csv_pandas_df = pd.read_csv(self.spec_files[i], sep=',', skiprows=1)
            self.velocities.append(csv_pandas_df.loc[:,'vel(km/s)'].to_numpy())
            self.intensities.append(csv_pandas_df.loc[:,'intensity(Jy)'].to_numpy())

        #Now have all of the lists of data
    
    def write_plots(self):
        #Write the splots out in the .spt format using the spectra data provided.

        #also extract the maximum and velocity
        for i in range(len(self.name_list)):
            print('Converting ' + str(self.name_list[i]))

            #Create and open the spt file
            plot = open(str(self.name_list[i]) + ".spt", "x")

            #First, insert the name and axis label components
            plot.write("heading \n " + str(self.name_list[i]) + "\n xlabel \n Velocity w.r.t. LSR (km/s) \n ylabel \n Flux (Jy) \n ")

            #Determine the largest flux value in the scan
            largest_flux_val = 0
            for j in range(len(self. intensities[i])):
                if self. intensities[i][j] > largest_flux_val:
                    largest_flux_val = self. intensities[i][j]

            
            #Extract the number of data points
            data_point_count = len(self. intensities[i])

            #Use the largest flux value to scale the plot, and if it is oo small increase to a sensible minimum
            if largest_flux_val < 30:
                largest_flux_val = 30
            flux_min = -0.2 * largest_flux_val
            flux_max = 1.2 * largest_flux_val

            #Get the upper and lower velocities from the end and start of the velocity array
            vel_lower = self.velocities[i][0]
            vel_upper = self.velocities[i][-1]

            #Now write these sections into the file
            plot.write("baxis\n " + str(vel_lower) + "\n " + str(vel_upper) + "\n taxis\n " + str(data_point_count) + "\n 1\n laxis\n " + str(flux_min) + "\n " + str(flux_max) + "\n data\n " + str(data_point_count) + "\n ")
            
            #Write the data into the file in the .spt format
            for j in range(len(self. intensities[i])):
                plot.write(str(self.velocities[i][j]) + "    " + str(self. intensities[i][j]) + " \n ")

            #Now add an end comment and close the file
            plot.write("end")
            plot.close()
            #Should be everything


converter = csv_to_spt()
converter.write_plots()
