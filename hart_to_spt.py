from astropy.coordinates import spectral_coordinate
from astropy.units.equivalencies import spectral
import numpy as np
import os, argparse, re
import astropy.time


#This script will take in the .txt spectra files from the Hart team and convert it to a .spt for easy compirosn of data.
#Will be able to produce as many spectra as there are in the input .txt
#Name format expected is Spec_SOURCENAME_FREQMHZ.txt

class hart_to_spt:
    def __init__(self):

        #Use argparse to accept the input arguments
        #To pass an input arugment, enter the flag followed by the value.
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "-spec",
            dest="spec_file",
            help="The input .txt spectra file from hart", 
            default=None #If no spectra is passed in, the program will quit
        )

        args = parser.parse_args()
        self.spec_file = args.spec_file
        self.termination_message = 'Script terminating'

        if self.spec_file == None:
            print('No spectra file passed.\n' + self.termination_message)
            quit()
        
        #Extract the spectra from the file using numpy loadtxt
        raw_text = np.loadtxt(self.spec_file)

        #raw_text[0] is two place holder values and then the velocity array
        #All subsequent raw_text[i] are the spectra contained in the file, with two points of data data at the start.

        #Extract the velocity, starting at array position 3 (python [2]) as the first two data points are placeholder text
        self.velocity_array = raw_text[0][2:]
        
        #Now go through the rest of the observations and extract the relevant data
        self.date_list = []
        self.spec_list = []
        for i in range(1, len(raw_text) - 1):
            #Convert the julian date to yyymmdd using astropy
            date_jd = astropy.time.Time(raw_text[i][0], format='mjd')
            date = date_jd.fits
            #Add the converted date to the list of dates
            self.date_list.append(date)
            #Add the spectra to the list of spectra
            self.spec_list.append(raw_text[i][2:])
        
        #Now have all spectra and date information extracted
    
    def write_plots(self):
        #Write the splots out in the .spt format using the spectra data provided.

        #Extract the required sections of the file name to generate the file name
        #Spec_Source_freq.txt, so filler _ source _ freq.txt (can't split .txt easily here because of 9.6 GHz)
        space_filler, source, freq_txt = re.split('[__]', str(self.spec_file))
        #Now split freq_txt into the frequency and the .txt
        freq, space_filler = re.split('[.]', str(freq_txt))

        #also extract the maximum and velocity
        vel_lower = self.velocity_array[0]
        vel_upper = self.velocity_array[-1] #-1 loops back around the last value, the max velocity

        for i in range(len(self.spec_list)):
            name = 'Hart_' + str(source) + '_' + str(freq) + '_' + str(self.date_list[i])
            plot = open(name + ".spt", "x")
            #First, insert the name and axis label components
            plot.write("heading \n " + name + "\n xlabel \n Velocity w.r.t. LSR (km/s) \n ylabel \n Flux (Jy) \n ")

            #Determine the largest flux value in the scan
            largest_flux_val = 0
            for j in range(len(self.spec_list[i])):
                if self.spec_list[i][j] > largest_flux_val:
                    largest_flux_val = self.spec_list[i][j]
            
            #Extract the number of data points
            data_point_count = len(self.spec_list[i])

            #Use the largest flux value to scale the plot, and if it is oo small increase to a sensible minimum
            if largest_flux_val < 30:
                largest_flux_val = 30
            flux_min = -0.2 * largest_flux_val
            flux_max = 1.2 * largest_flux_val

            #Now write these sections into the file
            plot.write("baxis\n " + str(vel_lower) + "\n " + str(vel_upper) + "\n taxis\n " + str(data_point_count) + "\n 1\n laxis\n " + str(flux_min) + "\n " + str(flux_max) + "\n data\n " + str(data_point_count) + "\n ")
            
            #Write the data into the file in the .spt format
            for j in range(len(self.spec_list[i])):
                plot.write(str(self.velocity_array[j]) + "    " + str(self.spec_list[i][j]) + " \n ")

            #Now add an end comment and close the file
            plot.write("end")
            plot.close()
            #Should be everything


converter= hart_to_spt()
converter.write_plots()