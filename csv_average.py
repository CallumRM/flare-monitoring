import argparse
import numpy as np
import pandas as pd

#Author: Callum Macdonald
#Direct any questions to either Callum.Macdonald@utas.edu.au or CallumRMacdonald@hotmail.com
#If another author contributes new methods to the script, please add your name and contact here


"""
A script designed to be able to average the four (or eight with ke) csv scans per source per obs.
Will cycle through each source and frequency and find all the files passed in with that combination, and then average them.
Output will be a new csv file, using the fm name and epoch from the first scan of that type encountered.
The date of the epoch will be fine, but the hour time will obviously not be relevant to an averaged scan.
If no scans for a particular source/frequencie combo are encountered, no 
"""

#TODO: Same as dynamic spectra script, need to add a checker to ensure that the velocities for the spectra
    #being averaged over are the same.

class csv_average:
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

        #The name to be used for the outputs and log file
        #For averaging within an obs, something like fm*** would make sense
        #For averaging accross obs, something like fm004_to_fm030 would make sense
        parser.add_argument(
            "-name",
            dest="name",
            type=str,
            default='no_name'
        )

        #The frequencies that will be search for and averaged over
        #As any combinations that are not found will be skipped, any future users of this script can 
            #add as many frequencies as they like
        parser.add_argument(
            "-freqs",
            dest="frequencies",
            nargs="+",
            default=[6181.128, 6668.5192, 7283.449, 7682.232, 7830.864, 11964.007, 12178.595, 12229.348, 12329.666] #7682 was replaced by 7283, but still have 7682 data
        )

        #The sources that will be searched for to average over
        #Same as frequencies, as many sources can be added as the user wants
        parser.add_argument(
            "-sources",
            dest="sources",
            nargs="+",
            default=['G328.237', 'G328.809', 'G351.417', 'G323.740', 'G318.948', 'G9.621'] 
        )

        #Extract the input arguments
        args = parser.parse_args()
        self.spectra = args.spectra
        self.name = str(args.name)
        self.frequencies = args.frequencies
        self.sources = args.sources

    def average_and_save(self, spectra, name, epoch, x_name = 'vel(km/s)', y_name = 'intensity(Jy)'):
        
        #Open the spectra and extract the vlocity and intensity data
        velocities = []
        intensities = []
        for i in range(len(spectra)):
            #Get the intensity and velicity data out for each spectra. Skiprows used to avoid the epoch data
            csv_pandas_df = pd.read_csv(spectra[i], sep=',', skiprows=1)
            velocities.append(csv_pandas_df.loc[:,'vel(km/s)'].to_numpy())
            intensities.append(csv_pandas_df.loc[:,'intensity(Jy)'].to_numpy())

        #Get the average intensity of the spectra
        average_intensity = np.zeros(len(intensities[0]))
        for i in range(len(intensities)):
            average_intensity += intensities[i]
        average_intensity = average_intensity / len(intensities)

        #All of the velocities should be the same, so just use the first one (which will always be present)
        velocity = velocities[0]

        #Create and open the csv file
        data = open(str(name) + ".csv", "x")

        #Put in the date metadata
        data.write('#' + str(epoch) + '\n')

        #Put in the headers
        #Eah value needs to have quotation marks around it
        data.write('"' + str(x_name) + '"' + ',' + '"' + str(y_name) + '"' + '\n')

        #Write out the velocity from the first spectra, and the average intensity
        for i in range(len(velocity)):
            data.write('"' + str(velocity[i]) + '"' + ',' + '"' + str(average_intensity[i]) + '"' + '\n')

        data.close()
    
    def proc_all(self):
        #Open a log file to keep track of which spectra were averaged
        log_file = open("averaging_log_" + str(self.name) + ".txt", "x")
        log_file.write("The source/frequency combinations that were not found are listed at the end of the file. \n\n")

        #Make a list to store the source/frequency combinations that were not used
        not_found_list = []
        #Cycle through each source and frequency
        for source in self.sources:
            for frequency in self.frequencies:

                #Make an empty list to store the spectra in, if any are found
                matched_spectra = []

                #Cycle through each spectra and see if they contain the source and frequency
                for i in range(len(self.spectra)):

                    #Use in to see if a string is contained within another, will return true if it is
                    if (str(source) in self.spectra[i]) and (str(frequency) in self.spectra[i]):
                        matched_spectra.append(self.spectra[i])
                
                #If matched_spectra is not empty, then run average_and_save on the spectra add information to the log file
                if len(matched_spectra) != 0:

                    #Extract the epoch from the first scan to be used in the averaged spectra
                    temp_file = open(matched_spectra[0], mode = 'r')
                    epoch_temp_string = temp_file.readline() #readline() reads the first line in a file
                    epoch = np.datetime64(epoch_temp_string[1:]) #Use string[1:] to drop the first char, which is the # used to indicate metadata
                    
                    #Make the name string for the averaged spectrum
                    spectrum_name = str(self.name) + '_' + str(source) + '_' +str(frequency)

                    #Print out the name and number of spectra
                    print('Now processing: ' + str(spectrum_name) + ' with ' + str(len(matched_spectra)) + ' spectra')

                    #Run average and save
                    self.average_and_save(spectra = matched_spectra, name = spectrum_name, epoch = epoch)

                    #Add the information to the log file
                    log_file.write(str(frequency) + ' ' + str(source) + ':' + '\n')
                    for i in range(len(matched_spectra)):
                        log_file.write('    ' + str(matched_spectra[i]) + '\n')
                    log_file.write('\n')
                
                else:
                    #If no spectra matching the source/frequency combination were found, add that combination to 
                        #the list of unfound combinations
                    not_found_list.append(str(frequency) + ' ' + str(source))
        
        #Once all of the combinations have been searched through, and those found processed, add unprocessed to the log file
        log_file.write('\nThe folowing frequency/source combinations were not found \n\n')
        for i in range(len(not_found_list)):
            log_file.write(not_found_list[i] + '\n')
        
        log_file.close()

averager = csv_average()
averager.proc_all()




