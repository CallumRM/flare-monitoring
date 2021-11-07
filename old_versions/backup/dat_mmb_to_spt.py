from astropy.coordinates import spectral_coordinate
from astropy.units.equivalencies import spectral
import numpy as np
import os, argparse, re
#import astropy.time
import pandas as pd


#Will take in spectra from MMB in .dat form and convert them to .spt for easy compirosn of data.
#Will give the new .spt files the same name as the input file.
#Code largly copyed from my hart_to_spt script.

class mmb_to_spt:
    def __init__(self):

        #Use argparse to accept the input arguments
        #To pass an input arugment, enter the flag followed by the value.
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "-spec",
            dest="spec_files",
            help="The input .txt spectra file from hart", 
            nargs="+" 
            #List of spectra
        )

        #Velocity range
        #The range of velocities to generate the plot for. First value is lower, second is upper limit
        #These will be centred on the source velocity.
        #eg entering [-20, 30] for a source centred at 55 will give a range of 35 to 85
        #TODO: Add an error catcher for case where the inputs are not used correctly
        parser.add_argument(
            "-vr",
            dest="velocity_range",
            nargs="+",
            default=[-20,20] #defualt set to standard maser zoom
        )

        args = parser.parse_args()
        self.spec_files = args.spec_files
        self.termination_message = 'Script terminating'

        if self.spec_files == None:
            print('No spectra file passed.\n' + self.termination_message)
            quit()
        
        #Extract the spectra from the files using pandas read_csv
        #Create empty lists for storing the names, velocity and intensity data for each file passed in.
        self.name_list = []
        self.velocity_data_list = []
        self.intensity_data_list = []
        for i in range(len(self.spec_files)):

            #Open the file using pandas 
            df_temp = pd.read_csv(self.spec_files[i], sep=' ', header=None)

            #Get the name of the file, need to drop the .dat
            temp_name = str(self.spec_files[i]).replace('.dat','')

            #Get the intesity and data values out. The first column is a row number.
            temp_vel = df_temp.loc[:,1].to_numpy()
            temp_intensity = df_temp.loc[:,2].to_numpy()

            self.name_list.append(temp_name)
            self.velocity_data_list.append(temp_vel)
            self.intensity_data_list.append(temp_intensity)

        #Now have all of the lists of data
    
    def write_plots(self):
        #Write the splots out in the .spt format using the spectra data provided.

        #also extract the maximum and velocity
        for i in range(len(self.name_list)):
            plot = open(str(self.name_list[i]) + ".spt", "x")
            #First, insert the name and axis label components
            plot.write("heading \n " + str(self.name_list[i]) + "\n xlabel \n Velocity w.r.t. LSR (km/s) \n ylabel \n Flux (Jy) \n ")

            #Determine the largest flux value in the scan
            largest_flux_val = 0
            for j in range(len(self.intensity_data_list[i])):
                if self.intensity_data_list[i][j] > largest_flux_val:
                    largest_flux_val = self.intensity_data_list[i][j]

            
            #Extract the number of data points
            data_point_count = len(self.intensity_data_list[i])

            #Use the largest flux value to scale the plot, and if it is oo small increase to a sensible minimum
            if largest_flux_val < 30:
                largest_flux_val = 30
            flux_min = -0.2 * largest_flux_val
            flux_max = 1.2 * largest_flux_val

            vel_lower = self.velocity_data_list[i][-1]
            vel_upper = self.velocity_data_list[i][0]

            #Now write these sections into the file
            plot.write("baxis\n " + str(vel_lower) + "\n " + str(vel_upper) + "\n taxis\n " + str(data_point_count) + "\n 1\n laxis\n " + str(flux_min) + "\n " + str(flux_max) + "\n data\n " + str(data_point_count) + "\n ")
            
            #Write the data into the file in the .spt format
            for j in range(len(self.intensity_data_list[i])):
                plot.write(str(self.velocity_data_list[i][j]) + "    " + str(self.intensity_data_list[i][j]) + " \n ")

            #Now add an end comment and close the file
            plot.write("end")
            plot.close()
            #Should be everything


converter= mmb_to_spt()
converter.write_plots()
