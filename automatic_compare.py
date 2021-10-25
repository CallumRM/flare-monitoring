import argparse
import numpy as np
import pandas as pd

#Author: Callum Macdonald
#Direct any questions to either Callum.Macdonald@utas.edu.au or CallumRMacdonald@hotmail.com
#If another author contributes new methods to the script, please add your name and contact here

"""
A script designed to automatically compare two sets of flare monitoring spectra.
The spectra will need to have standardised velocity samples (-sdv True in maser plot)
Only one spectra of each frequency / source combination will be used, so if multiple spectra exist, 
    they should be averaged
"""

#TODO: Add error catchers
#TODO: It might be good to be able to compare un-averaged spectra, but this will require some thought and 
    #time to implement cleanly. How will the output be structured?
#TODO: RMS for normalised changes depending on if there is a strong maser or not. Multiply RMS by normalisationfactor,
    #so that it can be compared within the epoch easier.

class automatic_compare:
    def __init__(self):

        #Use argparse to to accept the spectra and any customisation commands
        parser = argparse.ArgumentParser()


        #The new input spectra, to be compared to the older spectra passed in
        #Should be a list of csv files, named with at least the source and frequency 
        parser.add_argument(
            "-spec_new",
            dest="new_spectra",
            nargs="+"
        )

        #The old input spectra, used as a reference which the new spectra are compared to
        #Should be a list of csv files, named with at least the source and frequency 
        parser.add_argument(
            "-spec_old",
            dest="old_spectra",
            nargs="+"
        )

        #The name to be used for the log file
        #Sensible naming would include the obs code of the two
        parser.add_argument(
            "-name",
            dest="name",
            type=str,
            default='no_name'
        )

        #Sigma limit
        parser.add_argument(
            "-sig_lim",
            dest="sigma_limit",
            type=float,
            default=2
        )

        #Consecutive
        #The number of consecutive bins that need to be above the sigma limit for a detection to be flagged
        #e.i. if this is set to 4, and 3 bins in a row are above the limit, but the 4th isn't, then no detection will be flagged
        #TODO: Test what is a reasonable number. Can be balanced with sigma to get false positive down.
            #Main issue will be phasecal tones if they are left in.
        parser.add_argument(
            "-consec",
            dest="consecutive",
            type=int,
            default=10
        )

        #Velocity range
        #The velocity range around the centre velocity which will be searched for detections
        #Tuple of floats, to allow for asymetric searching
        #First value should be negative, unless you want to have a window offset from the centre velocity
        parser.add_argument(
            "-vr",
            dest="velocity_range",
            nargs="+",
            default=[-10,10] #Default set for +-10 km/s around centre
        )

        #Normalise
        #A bool used to determine if the spectra should be normalised before being compared
        #This will likely reduce the flase positive rate
        parser.add_argument(
            "-norm",
            dest="normalise",
            type=str,
            default='False' #Not normalised by default
        )

        #The frequencies that will be search for and compared for each source
        #As any combinations that are not found will be skipped, any future users of this script can 
            #add as many frequencies as they like
        parser.add_argument(
            "-freqs",
            dest="frequencies",
            nargs="+",
            default=[6181.128, 6668.5192, 7283.449, 7682.232, 7830.864, 11964.007, 12178.595, 12229.348, 12329.666] #7682 was replaced by 7283, but still have 7682 data
        )

        #The sources that will be searched for and compared at each frequency
        #Same as frequencies, as many sources can be added as the user wants
        parser.add_argument(
            "-sources",
            dest="sources",
            nargs="+",
            default=['G328.237', 'G328.809', 'G351.417', 'G323.740', 'G318.948', 'G9.621'] 
        )

        #Extract the input arguments
        args = parser.parse_args()
        self.new_spectra = args.new_spectra
        self.old_spectra = args.old_spectra
        self.name = str(args.name)
        self.sigma_limit = float(args.sigma_limit)
        self.consecutive = int(args.consecutive)
        self.velocity_range = args.velocity_range
        self.normalise = str(args.normalise) #Not a real bool, truth inputs are strange, so strings are easier
        self.frequencies = args.frequencies
        self.sources = args.sources
        self.termination_message = 'Script terminating'


    #Method to generate the RMS noise level for data
    #This function is a slightly modified version of a function I wrote for maser_plot.py
    #Will take an intensity array, velocity array, centre velocity, exclusion size and a inclusion size
    #The exlusion min and max correspond to the velocities between which the data should be excluded (where the feature is)
    #The inclusion size will be centred on the centre velocity and everything within it other than the excluded region will
        #be used to determine the RMS noise level
    #Will also return the mean intensity from the calculation region, which can be used to centre the spectra on 0
    #Values should be given in km/s
    #The intensity and velocity arrays should be the same length
    #The exclusion size should be smaller than the inclusion size
    def calculate_rms_and_mean(self, intensity_array, velocity_array, centre_velocity, exclusion_min, exclusion_max, inclusion_size):
        
        #Check that the intensity and velocity arrays are the same size
        if len(intensity_array) != len(velocity_array):
            print('Length of intensity and velocity arrays did not match while calculating RMS noise \n' + self.termination_message)
            quit()
        
        #Check that the exclusion area is smaller than the inclusion area
        if (centre_velocity - exclusion_min >= inclusion_size) or (exclusion_max - centre_velocity) >= inclusion_size:
            print('Exclusion size must be smaller than inclusion size for calculating RMS noise \n' + self.termination_message)
            quit()
        
        #Calculate the velocity bounds for inclusion
        vel_inc_low = centre_velocity - inclusion_size
        vel_inc_up = centre_velocity + inclusion_size

        #Create a new array with only the relevant intensity points in it
        #Order is not important for this calculation
        flagged_intensity = []
        for i in range(len(intensity_array)):
            #Check that the velocity is within inclusion and not exclusion
            if (velocity_array[i] > vel_inc_low and velocity_array[i] < vel_inc_up) and not(velocity_array[i] > exclusion_min and velocity_array[i] < exclusion_max):
                #Add the corresponding intensity value to the new intensity array
                flagged_intensity.append(intensity_array[i])

        #Centre the intensity values on zero
        mean_intensity = np.mean(flagged_intensity)
        flagged_intensity_zeroed = flagged_intensity - mean_intensity

        #Use the zeroed intensity array to calculate the RMS
        squared_intensity = flagged_intensity_zeroed**2
        mean_square = np.mean(squared_intensity)
        RMS = np.sqrt(mean_square)

        return (RMS, mean_intensity)
    
    #Look if anything is above x sigma
    #Then use the list of detections to print them out, followed by non-detections
    #Can have a search range, where window is centrered around vel[len/2]
    #Can have consecutive requirement, need to have x bins in a row before it calls it a detection
    #With consecutive, have it keep running, and when it stops, can list the first and last velocity to give an idea of what was found
        #That way -5 to 5 lets someone know that it was probably a fuckup
    #TODO: This method is meaty, might be good to break it up into some smaller methods
    def compare_all(self):

        #Open the log file to store the output information. List the parameters
        log_file = open("comparison_log_" + str(self.name) + ".txt", "x")
        log_file.write("Comparison parameters:\n")
        log_file.write("    Sigma multiple threshold: " + str(self.sigma_limit) + "\n")
        log_file.write("    Consecutive requirement number: " + str(self.consecutive) + "\n")
        log_file.write("    Velocity range around centre: " + str(self.velocity_range) + "\n")
        log_file.write("    Normalised: " + str(self.normalise) + "\n")
        log_file.write("    Sources searched for: " + str(self.sources) + "\n")
        log_file.write("    Frequencies searched for: " + str(self.frequencies) + "\n\n")

        #Create lists to hold the detected and undectected combinations, as well as the unfound combinations
        #The detected list will be tuples, with the freq/source as one entry, and a list of touples of found veloicty ranges as the second
        detected = []
        #undetected and unfound lists will just have the freq/source
        undetected = []
        unfound =[]
        

        #Cycle through each source and frequency
        for source in self.sources:
            for frequency in self.frequencies:
                
                #Bools to test if a spectra has been in both lists for a frequency / source combination
                #Also allows checking of single spectra condition
                found_new = False
                found_old = False
                
                #See if the frequency source combination can be found in the new spectra. If it is, save it and set the bool to true
                    #If mutliple are found, quit the program, as it will break later steps.
                for i in range(len(self.new_spectra)):
                    
                    #Use in to see if a string is contained within another, will return true if it is
                    if (str(source) in self.new_spectra[i]) and (str(frequency) in self.new_spectra[i]):
                        
                        #Check if another spectra at this frequency and source has been found already
                        if found_new == False:
                            new_spectrum = (self.new_spectra[i])
                            found_new = True
                        else:
                            print('Two spectra with source: ' + str(source) + ' and frequency: ' + str(frequency) + ' found in new spectra list.')
                            print(self.termination_message)
                            quit()
                
                #Repeat the above process for the old spectra
                #This could probably be made a function, and then old and new passed in, but just repeating works too.
                for i in range(len(self.old_spectra)):

                    #Use in to see if a string is contained within another, will return true if it is
                    if (str(source) in self.old_spectra[i]) and (str(frequency) in self.old_spectra[i]):
            
                        #Check if another spectra at this frequency and source has been found already
                        if found_old == False:
                            old_spectrum = (self.old_spectra[i])
                            found_old = True
                        else:
                            print('Two spectra with source: ' + str(source) + ' and frequency: ' + str(frequency) + ' found in new spectra list.')
                            print(self.termination_message)
                            quit()

                #Now, if both have been found, comparison can begin. Otherwise, note that the spectra wasn't sound for at least one
                if (found_new == True) and (found_old == True):

                    print('Spcetra found for source: ' + str(source) + ' and frequency: ' + str(frequency))

                    #Open the two spectra. Skiprows for the epoch metadata
                    new_csv_pandas_df = pd.read_csv(new_spectrum, sep=',', skiprows=1)
                    old_csv_pandas_df = pd.read_csv(old_spectrum, sep=',', skiprows=1)

                    #Extract the velocity and intensity data
                    #Need to use np.around() because the script was randomaly adding 1e-8 to some values, so caomparison would break later
                    new_velocity = np.around(new_csv_pandas_df.loc[:,'vel(km/s)'].to_numpy(), decimals = 5)
                    new_intensity = new_csv_pandas_df.loc[:,'intensity(Jy)'].to_numpy()

                    old_velocity = np.around(old_csv_pandas_df.loc[:,'vel(km/s)'].to_numpy(), decimals = 5)
                    old_intensity = old_csv_pandas_df.loc[:,'intensity(Jy)'].to_numpy()

                    #Check that both the velocity arrays are the same, otherwise comparison cannot take place
                    #Use some cool numpy shit, where == makes a truth array, then .all does a complete and of the array
                    velocity_comparison = new_velocity == old_velocity
                    if velocity_comparison.all() == False:
                        print('Velocity arrays did not match.')
                        print(self.termination_message)
                        print(new_velocity - old_velocity)
                        quit()
                    
                    #If the spectra are to be normalised, do it here
                    if self.normalise == 'True':

                        #Find the maximum intensity values for both spectra
                        largest_flux_new = 0
                        largest_flux_old = 0
                        for i in range(len(new_intensity)):
                            if abs(new_intensity[i]) > largest_flux_new:
                                largest_flux_new = abs(new_intensity[i])
                            if abs(old_intensity[i]) > largest_flux_old:
                                largest_flux_old = abs(old_intensity[i])
                        
                        #Normalise the spectra
                        new_intensity = new_intensity / largest_flux_new
                        old_intensity = old_intensity / largest_flux_old
                        
            
                    #Now take the difference of the two spectra.
                    difference_spectrum = new_intensity - old_intensity
                    
                    #Calculate the max and min velocities to be searched through based on provided range
                    #The search velocities will be excluded from the rms and mean calculations
                    centre_velocity = new_velocity[int(len(new_velocity)/2)] #Centre velocity value from middle of vel array
                    lower_vel = centre_velocity + float(self.velocity_range[0])
                    upper_vel = centre_velocity + float(self.velocity_range[1])
                    
                    #Now calculate the rms of the difference spectrum, which will be used to see 
                        #if any differences are significant
                    #Do this using the calculate rms function, as copied from maser_plot
                    #Inclusion size for this script is just the whole velocity array,
                        #function uses +-, so inclusion size is half of velocity
                    inclusion_size = abs((new_velocity[-1] - new_velocity[0])/2)

                    difference_rms, intensity_mean = self.calculate_rms_and_mean(intensity_array = difference_spectrum, velocity_array = new_velocity, centre_velocity = centre_velocity, exclusion_min = lower_vel, exclusion_max = upper_vel, inclusion_size = inclusion_size)
                    
                    #Subtract the mean from the difference spectrum, so that the rms is representative of the noise in it
                    difference_spectrum_zeroed = difference_spectrum - intensity_mean

                    #Get the cuttoff value by multiplying the rms by the sigma limit
                    detection_threshold = self.sigma_limit * difference_rms

                    #Bool to track if a change has been detected
                    found_change = False

                    #Counter for the number of consecutive detections made, and a list to store tuples of velocities in
                    conscutive_counter = 0
                    detected_velocity_ranges = []

                    #Now go through the array and see if each point is above the rms multiple value
                    #If >= consecutive are found, then it will be noted
                    #For notes on the process for consecutive, see the bottom of this script 
                        #Complicated to explain in parts
                    #TODO: IF RUN OUT OF VEL AND CONSEC IS > MIN, NEED TO SAVE
                    for i in range(len(new_velocity)):
                        
                        #Only interested in values within the checking range
                        if (new_velocity[i] > lower_vel) and (new_velocity[i] < upper_vel):
                            #If it is above the threshold value, need to note on the consecutive counter
                            #print(str(new_velocity[i]) + '  ' + str(abs(difference_spectrum_zeroed[i])) + '  ' + str(conscutive_counter))
                            if abs(difference_spectrum_zeroed[i]) > detection_threshold:
                                #If this is the first detection since a non-detection, need to mark the velocity as
                                    #a potential sarting velocity
                                if conscutive_counter == 0:
                                    temp_starting_velocity = new_velocity[i]

                                #Regardless of if it is the first or not, need to iterate the consecutive counter
                                conscutive_counter += 1
                            
                            #If it is not above the threshold, need to check if this is the end of a 
                                #chain of detections above the minimum consecutive number
                            else:
                                #If this is the first low value after a row of detections, save the previous velocity
                                    #and change the detection bool to true
                                if conscutive_counter >= self.consecutive:
                                    temp_ending_velocity = new_velocity[i-1] #-1 is fine because this cannot occur for i = 0

                                    #Add the velocity range to the detected velocity range list
                                    temp_velocity_range_string = '(' + str(temp_starting_velocity) + ', ' + str(temp_ending_velocity) + ')'
                                    detected_velocity_ranges.append(temp_velocity_range_string)

                                    found_change = True
                            
                                #Regardless of whether or not this is the end of a chain longer than the minimum, need to
                                    #reset the consecutive counter to 0
                                conscutive_counter = 0
                    
                    #Now that the loop is done, check if a detection was made, and add the result to the appropriate list
                    if found_change == True:
                        detected.append(str(frequency) + ' ' + str(source) + '(RMS: ' + str(difference_rms) + ', centre velocity: ' + str(centre_velocity) + ') detected at the following velocities: ' + str(detected_velocity_ranges))
                        print('     Difference found for ' + str(frequency) + ' ' + str(source))
                    #No detection, just note combination
                    else:
                        undetected.append(str(frequency) + ' ' + str(source) + '(RMS: ' + str(difference_rms) + ')')
                
                #If no spectra matching the source/frequency combination were found, add that combination to 
                        #the list of unfound combinations
                else:
                    unfound.append(str(frequency) + ' ' + str(source))

        #Once all combinations have been cycled through, can write out the log file.
        log_file.write('Number of spectra where a change matching the above criteria was found: ' + str(len(detected)) + '\n')
        log_file.write('Number of spectra where no change matching the above criteria was found: ' + str(len(undetected)) + '\n')
        log_file.write('Number of combinations where a spectrum was not found for both old and new: ' + str(len(unfound)) + '\n\n')

        log_file.write('The following combinations had significant detections: \n')
        for i in range(len(detected)):
            log_file.write(detected[i] + '\n\n')
        
        log_file.write('The following combinations did not have significant detections: \n')
        for i in range(len(undetected)):
            log_file.write(undetected[i] + '\n')
        
        log_file.write('\nThe following combinations were not found: \n')
        for i in range(len(unfound)):
            log_file.write(unfound[i] + '\n')
        
        log_file.close()

comparer = automatic_compare()
comparer.compare_all()


"""
Consecutive implementation note
Want to be able to give the start and stop velocity for a range of values consecutively above the rms theshold value
To do this, when a large value is found, the consecutive counter is iterated. If the counter was on zero, it will also
    note the current velocity value as a potential starting velocity
When a value below the theshold is encountered, it will check if the current consecutive counter is above the minimum
If it is not, then it is reset to zero and the process continues
If it is above the minimum, then the previous velocity is paired with the starting velocity for that consecutive run and
    the two values are added as a tuple to the list of detections for this combo
"""