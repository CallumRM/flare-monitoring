from astropy.time.formats import TimeStardate
import numpy as np
import os, argparse, re
import errno
import math
#import sys
import PyAstronomy.pyasl
from astropy import units
import astropy.time
import astropy.coordinates
from numpy.lib.polynomial import poly
#from numpy.core.fromnumeric import mean
#from numpy.lib.function_base import _calculate_shapes

#Author: Callum Macdonald
#Direct any questions to either Callum.Macdonald@utas.edu.au or CallumRMacdonald@hotmail.com
#If another author makes contributions of new methods, please add your name and contact here

"""
A script that will take some parameters and a list of spectra .bin files, and then produce a plot
Will either produce a simple on-off plot, where the on and off spectra are selected from the list
    or, will produce an spectra normalised by the average of the other spectra
    or, will just plot a single spectra.

"""
#TODO: Add proper help print out
#TODO: Perl output, Do on-off or baseline for multiple inputs. Combine scans, should be easy, just need to convert vel before merging
#Will need to take in the start time of the observing epcoh in numpy YYYY-MM-DD'T'HH:MM:SS form, and then have scan length, use numpy datetime, and then use timedelta to add extra
#This way all data can be velocity corrected before before being plotted.
#TODO: Need to set intensity, will need to fit guassian and then use that to find conversion factor
    #Or use other intensity calibration that Simmon talked about
#TODO: Set up flexiblity for different observation patterns. Use naming to identify intenisty calibration scans.
#TODO: Make a script that can read a schedule file and produce the inputs required from it (list of source names with on off stuff)
#TODO: Make option to print out the baseline fit. Have it name it properly. Can also have naming change if baseline removed
#TODO: Add overide variable that can be used to have easy 'profiles' for cd and ef under standard fm obs.
#TODO: Need to add a warning message if the number of scans entered is not a multiple of the obs list
class maser_plot:
    def __init__(self):

        #Use argparse to accept the input arguments
        #To pass an input arugment, enter the flag followed by the value.
        #For example, to pass in a bandwidth, maser_plot.py -bw 32 ...
        #For arguments with multiple components, like -spec, just add a space between each input.
        #EG -spec file_one file_two file_three ...
        #Arguments like this can accept wildcard in linux, EG -spec file_*

        #Parser will hold the argument data, will be extracted into other variables below.
        parser = argparse.ArgumentParser()

        #Bandwidth
        #Bandwidth in MHz
        #If this is updated in future, change the processAll.py script to add bandwidth information into the 
        #   name, and then extract it for each scan as is done for many other parameters
        parser.add_argument(
            "-bw",
            dest="bandwidth",
            help="observing bandwidth in MHz", 
            type=float, 
            default=32 #default set to 32 for the AuScope 12m antennas
        )

        #Output Format
        #The format for the spectra to be saved to.
        #Currently have
            #spt: .spt format for use with Simon's Splot program.
            #csv: .csv format for easy further processing. csv file will have a header '#epoch' so that the date can be used for further processing 
        parser.add_argument(
            "-format",
            dest="format",
            help="Format of the output data", 
            type=str, 
            default='spt' #default set to base inifile name for e and f
        )

        #List of spectra
        #Will accept a single spectrum or a list of spectra
        #If a single spectrum is to be plotted, only one spectrum needs to be passed in
        #For a simple on-off, both the on and off spectra must be passed in
        #For the on averaged-off, at least two spectra need to be passed in (for sensible result, need at least ~5)
        parser.add_argument(
            "-spec",
            dest="spec_list",
            nargs="+",
            default=[] #This argument obviously needs to be passed in, and if it isn't then an error will be returned
        )

        #Intensity scale
        #The scale factor to convert on-off quotient values to intensity
        #Default set roughly for G9.621 in May
        parser.add_argument(
            "-int_s",
            dest="intensity_scale",
            type=float,
            default=3500 / 0.9
        )

        #Epoch of the observation
        #The date of the observation as a string in fits format'YYYY-MM-DDTHH:MM:SS'
        #Example '2021-03-31T01:30:00'
        #This will be used to remove the relative velocity of the earth with respect to LSR
        #TODO: add erorr catcher
        parser.add_argument(
            "-epoch",
            dest="epoch",
            default=None
        )

        #Observation Length
        #The length of each observation in seconds
        #Default set for 10 minutes plus 30 sec for switch over time
        #Used to keep track of source velocities as earth rotates, as the time provided with 
        #   epoch is only for the first scan
        parser.add_argument(
            "-length",
            dest="observation_length",
            default=630 #10 mins and 30 sec between observations
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
            default=[-15,15] #defualt set to standard maser zoom
        )

        #Velocity inclusion size
        #This is the width around the central velocity used for calculating the rms and baseline, in km/s
        #In these calculations, the exclusion window specified in the source dictionary will be excluded, 
        #   and everything else within the inclusion velocty will be included in the calculations
        parser.add_argument(
            "-inc",
            dest="velocity_inclusion_size",
            type=float,
            default=20 #defualt set to a reasonable window. Larger windows will require a higher order polynomial
        )
        
        #Baseline polynomial order
        #The order of the polynomial used in the baseline fitting process
        #The size of the inclusion winodw should be considered when choosing an order. If a larger window is 
        #   used, a higher order polynomial will be required for a good fit
        parser.add_argument(
            "-bpo",
            dest="baseline_polynomial_order",
            type=int,
            default=7 #defualt set to 7 for a reasonable fit over 40 km/s
        )

        #Remove baseline
        #A bool to used to determine if the baseline should be removed
        #Baseline will be removed if true, ortherwise the spectra will be plotted as is
        parser.add_argument(
            "-rbl",
            dest="remove_baseline",
            type=bool,
            default=True #Baseline removed by default
        )

        #Rounding option
        #A correction used for the early fm observations where the proc file was set up without full frequency specification
        #To get the correct velocity results, need to round the maser frequency to match what was entered in the proc file
        #TODO: Can this be fixed for future obs? 
        parser.add_argument(
            "-round",
            dest="frequency_rounding",
            type=bool,
            default=False #defualt set to False, as this will hopefully only be needed for legacy processing
        )

        #Observation order
        #The order of the sources observed, using names which are found in the dictionary below
        #TODO: Make a nice external reference to this dictionary.
        #The names of the observations in a list (list of strings)
        #Default set for fm observations of 6 sources in on off format
        #TODO Add something like group name, so that fm*** can be added to the end of each name.
        #CHECK THAT THIS IS THE WAY IT WORKS: Scan's named 'off' will not be plotted in single plot version, 
        parser.add_argument(
            "-order",
            dest="observation_order",
            nargs="+",
            default=['G328.237', 'off', 'G328.809', 'off', 'G351.417', 'off', 'G323.740', 'off', 'G318.948', 'off', 'G9.621', 'off'] 
        )

        #Plotting option argument
        #-po "number". 1 will just plot a single spectra, 2 will plot on-off, 3 will plot on-off
            #with the average of other spectra used as off
        parser.add_argument(
            "-po",
            dest="plot_option",
            help="The plotting option",
            type=int,
            default=2 #default to 2, as fm obs are by default on-off
        )

        #If Channel Frequencies
        #The frequencies in MHz of the lines being observed in the different IF and channel combinations
        #A tuple of tuples. Each of the secondary tuples will contain the frequency values for the channels of the IF associated with that tuple
        #For example, if IF A and B are used, both observing 2 channels, with A having 1 and 2 MHz, and 3 and 4 MHz for B
        #    then the structure would be ((1,2), (3,4)).
        #If some IFs are not used (Like A and B for FM), then just use None values in their place.
        #Need to include spacing for all IFs up to the ones used, as this structure will be searched to find the frequencies of interest
        parser.add_argument(
            "-freqs",
            dest="IF_chan_freqs",
            help="The frequencie of transitions for each IF and channel. See the script for further details",
            nargs="+",
            default=(None, None,(6181.128, 6668.5192, 7283.449, 7830.864),(6181.128, 6668.5192, 7283.449, 7830.864),(11964.007, 12178.595, 12229.348, 12329.666),(11964.007, 12178.595, 12229.348, 12329.666), None, None) 
            #Set up for FM observations (A,B not used, four channels in each of cdef, cd and ef used with same freqs, as they are complamentary polaristions), no higher used (left in as example though)
            #For obs before fm023, need to exchnage 7.2 with 7.68 GHz transition)
        )

        #End of input arguments
        #Now extract variables out of the parser for simplicity
        #TODO: Have a print out of the entered parameters
        args = parser.parse_args() #args now is a list of all of the input variables
        self.bandwidth = int(args.bandwidth) * 1e6 #Convert bandwidth to Hz
        self.spec_list = args.spec_list
        self.format = str(args.format)
        self.intensity_scale = float(args.intensity_scale)
        self.epoch = str(args.epoch)
        self.observation_length = int(args.observation_length)
        self.observation_order = args.observation_order
        self.plot_option = int(args.plot_option)
        self.IF_chan_frequencies_MHz = args.IF_chan_freqs
        self.min_plot_vel = int(args.velocity_range[0]) #extract the minimum and maximum velocity to plot. Need int, as it reads as string
        self.max_plot_vel = int(args.velocity_range[1])
        self.velocity_inclusion_size = float(args.velocity_inclusion_size)
        self.baseline_polynomial_order = int(args.baseline_polynomial_order)
        self.remove_baseline_bool = bool(args.remove_baseline)
        self.frequency_rounding = bool(args.frequency_rounding)
        self.termination_message = 'Script terminating' #Message to be used when the script encounters a potential error and quits. Can be changed later to something more interesting.

        #TODO: Add error catchers for all of the other values (ie plot option, spec_pos out of bounds)
        #If the spec_list is empty, there is no data to be plotted, and the program quits
        if len(self.spec_list) < 1:
            print('Must have at least one spectrum passed using -spec')
            print(self.termination_message)
            quit()

        #Define some dictionaries for use in some of the functions used to plot
        
        
        #Dictionary of the observatories around Australia that are being used or may be used for this project
        #Each entry is a tuple of the form (lat, long, alt), with the lat and long in degrees, and the alt in metres
        #Location data from (J. Lovell et.al 2012, The AuScope geodetic VLBI array, http://www.phys.utas.edu.au/physics/mt_pleasant_observatory.html and https://www.exploroz.com/places/103554/sa+ceduna-radio-astronomy-observatory)
        #Converted into decimal degrees
        self.observatories = {
            'hb' : ( -42.80557, 147.4381, 40.967), #Hobart 12m
            'ke' : (-14.37546, 132.15237, 189.262), #Katherine 12m
            'yg' : (-29.04715, 115.34563, 248.236), #Yarragadee 12m
            'ho' : (-42.805, 147.439, 43), #Hobart 26m (info not specific, might be worth finding good source)
            'cd' : (-31.86883, 133.809799, 154), #Ceduna 30m
            'parkes' : (-32.99839, 148.26352, 415), #Parkes 64m (Murriyang) (For vel calibration tests using ATNF online tool, sourced from Parkes user guide)
        }


        #TODO: For multiple intensity calibrations will likely need to include them here.
        #Maybe not, might be able to do that in a seperate way and just use the specified positions to determine
            #intensity for a particular scan

        #For ease of re-ordering observations, each source name can be entered as a list, and all of the details 
            #will be extracted using the below dictionary
        #If another naming convention makes more sense, this can be changed easily
        #Using the information here, a list of velocities and positions will be made, which is easier to reference later
        #Each entry in the dictionary has the name as a label, and a tuple of the form (vel_centre, vel_low, vel_high, RA, Dec, full_name)
        #The upper and lower velocities will be used to for baseline fitting and RMS calculations to exclude the spectral features
        #The RA is in the form **h**m**.*s, and the Dec is in the form ***d**m**.*s (hour, min, sec and deg, min, sec respectively)
        #All information from Methanol Multibeam survey
        #Some exlcusion windows narrowed
        self.source_params = {
            'off' : (0, -5, 5, '0h0m0s', '0d0m0s', 'none'),
            'G328.237' : (-44.7, -47.0, -31.5, '15h57m58.28s', '-53d59m22.7s', 'G328.237-0.547'),
            'G328.809' : (-44.5, -47.5, -42.0, '15h55m48.70s', '-52d43m05.5s', 'G328.809+0.633'),
            'G351.417' : (-10.4, -12.0, -4.5, '17h20m53.37s', '-35d47m01.2s', 'G351.417+0.645'),
            'G323.740' : (-50.5, -59.0, -42.0, '15h31m45.45s', '-56d30m50.1s', 'G323.740-0.263'),
            'G318.948' : (-34.6, -39.0, -31.0, '15h00m55.40s', '-58d58m52.1s', 'G318.948-0.196'),
            'G9.621' : (1.3, -3, 7, '18h06m14.67s', '-20d31m32.4s', 'G9.621+0.196')
        }


        #Check if the input .bin spectra files exist, and create a list of file names and a list of the paths to these spectra
        #The name will be used to determine which scan, IF and channel the spectra is of.
        #This process is largely based on Guifre's code
        self.file_name_list = []
        self.path_list = []
        for i in range(len(self.spec_list)):
            if os.path.exists(self.spec_list[i]):
                path, file_name = os.path.split(self.spec_list[i]) #split the path and file name
                self.file_name_list.append(file_name)
                self.path_list.append("{}/".format(os.path.abspath(path)))
            else:
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), self.spec_list[i])
    

    #The following function was NOT written by me (Callum Macdonald) and was sourced from github.
    #This function returns the barycentric velocity correct, which is used in the velocity_correction function below
    #Function from https://gist.github.com/StuartLittlefair/5aaf476c5d7b52d20aa9544cfaa936a1 on github
    def bary_vel_cor(self, time, skycoord, location=None):
        """Barycentric velocity correction.
  
        Uses the ephemeris set with  ``astropy.coordinates.solar_system_ephemeris.set`` for corrections. 
        For more information see `~astropy.coordinates.solar_system_ephemeris`.
  
        Parameters
        ----------
        time : `~astropy.time.Time`
        The time of observation.
        skycoord: `~astropy.coordinates.SkyCoord`
        The sky location to calculate the correction for.
        location: `~astropy.coordinates.EarthLocation`, optional
        The location of the observatory to calculate the correction for.
        If no location is given, the ``location`` attribute of the Time
        object is used
    
        Returns
        -------
        vel_corr : `~astropy.units.Quantity`
        The velocity correction to convert to Barycentric velocities. Should be added to the original
        velocity.
        """
  
        if location is None:
            if time.location is None:
                raise ValueError('An EarthLocation needs to be set or passed '
                         'in to calculate bary- or heliocentric '
                         'corrections')
            location = time.location
    
        # ensure sky location is ICRS compatible
        if not skycoord.is_transformable_to(astropy.coordinates.ICRS()):
            raise ValueError("Given skycoord is not transformable to the ICRS")
  
        ep, ev = astropy.coordinates.solar_system.get_body_barycentric_posvel('earth', time) # ICRS position and velocity of Earth's geocenter
        op, ov = location.get_gcrs_posvel(time) # GCRS position and velocity of observatory
        # ICRS and GCRS are axes-aligned. Can add the velocities
        velocity = ev + ov # relies on PR5434 being merged
  
        # get unit ICRS vector in direction of SkyCoord
        sc_cartesian = skycoord.icrs.represent_as(astropy.coordinates.UnitSphericalRepresentation).represent_as(astropy.coordinates.CartesianRepresentation)
        return sc_cartesian.dot(velocity).to(units.km/units.s) # similarly requires PR5434



    #Method to calculate the velocity correction given an observatory location, epoch and Ra and Dec
    #lat long in degrees, alt in metres
    #epoch numpy datetime64 object
    #RA string in **h**m**.*s format
    #Dec string in ***d**m**.*s
    #Will return a float of the velocity correction
    def velocity_correction(self, lat, long, alt, epoch, RA, dec):
        #Convert the epoch, sky coordinates and observation location into the formats required for bary_vel-cor
        epoch_astro = astropy.time.Time(epoch)
        epoch_julian = epoch_astro.jd
        #print(epoch_astro)
        #print(epoch_julian)
        #print(RA)
        #print(dec)
        sky_coords = astropy.coordinates.SkyCoord(RA, dec, frame = 'icrs')
        obs_location = astropy.coordinates.EarthLocation.from_geodetic(long, lat, height=alt)

        #Use these values to calculate the barycentric velocity correction using the github sourced function bary_vel_cor
        bary_vel_output = self.bary_vel_cor(time=epoch_astro, skycoord = sky_coords, location = obs_location)

        #The formatting of this output is not just a float value, and it needs to be converted to one for later
        #First, get the velocity string out of the output
        bary_vel_string, space_filler , space_filler , space_filler = re.split('[ ]', str(bary_vel_output))

        #Now convert the string to a float, and change it to recession velocity (negative)
        bary_velocity = -float(bary_vel_string)
        #print('Old Bary ' + str(bary_velocity))




        #TEMPORARY CODE TO COMPARE BARY VEL
        #Convert RA and dec into deg
        d_dec, m_dec, s_dec, space_filler = re.split('[dms]', str(dec))
        dec_deg = (float(d_dec) + float(m_dec)/60 + float(s_dec)/(60*60))

        h_ra, m_ra, s_ra, space_filler = re.split('[hms]', str(RA))
        ra_deg = (float(h_ra) + float(m_ra)/60 + float(s_ra)/(60*60)) * 360 / 24

        helio_velocity = PyAstronomy.pyasl.helcorr(obs_long = long, obs_lat = lat, obs_alt = alt, ra2000 = ra_deg, dec2000 = dec_deg, jd = epoch_julian)
        bary_hel_velocity = PyAstronomy.pyasl.baryCorr(ra = ra_deg, dec = dec_deg, jd = epoch_julian)

        #print('Helio ' + str(helio_velocity))
        #print('Bary corr ' + str(bary_hel_velocity))
        #print(epoch_julian)

        bary_correction = bary_hel_velocity[1] - bary_hel_velocity[0]
        corrected_bary_velocity = bary_velocity + bary_correction
        #print('Bary ' + str(corrected_bary_velocity))
        ###########


        #To convert to LSR, define the direction of the LSR, then take the dot product of the unit vectors
        #This result can be multiplied by 20 km/s to get LSR velocity conversion factor
        #TODO: Find a source for these values, currently just using the script that Jim wrote, given to me by Simon
        #Also given on ATNF sky freq calculator, 20 km/s in direction RA 18 hour (270 deg), Dec + 30 deg, needs b1900 equinox, so use fk4
        lsrk_epoch = astropy.time.Time('b1900')
        lsrk_sky_coords = astropy.coordinates.SkyCoord(270*units.deg, 30*units.deg, frame = 'fk4', equinox=lsrk_epoch, obstime = epoch_astro)

        #Skycoords.separation returns the separation between two skycoord objects in degrees        
        dms_separation = sky_coords.separation(lsrk_sky_coords)

        #Convert this to degrees, then radians
        #Use re.split to get the degree, min, sec components, then add them, with the mins and seconds converted to deg, then convert to rad
        d_sep, m_sep, s_sep, space_filler = re.split('[dms]', str(dms_separation))
        separation = (float(d_sep) + float(m_sep)/60 + float(s_sep)/(60*60)) * np.pi / 180
        #print(separation * 180 / np.pi)

        #Now the separation can be used to calculate the dot product between the two, and this can be multiplied by 20 km/s to give LSR vel correction
        lsr_correction = -20 * np.cos(separation)
        #print('LSR ' + str(lsr_correction))

        #Combine the two velocities to get the LSR corrected velocity
        LSR_velocity = bary_velocity + lsr_correction
        #print(LSR_velocity)

        return LSR_velocity
    

    #Method to generate the RMS noise level for data
    #Will take an intensity array, velocity array, centre velocity, exclusion size and a inclusion size
    #The exlusion min and max correspond to the velocities between which the data should be excluded (where the feature is)
    #The inclusion size will be centred on the centre velocity and everything within it other than the excluded region will
        #be used to determine the RMS noise level
    #Values should be given in km/s
    #The intensity and velocity arrays should be the same length
    #The exclusion size should be smaller than the inclusion size
    #TODO: Implement a way for different exclusions to be used based on the source
    def calculate_rms(self, intensity_array, velocity_array, centre_velocity, exclusion_min, exclusion_max, inclusion_size):
        
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

        return RMS


    #Some of the code in the following section is based on the code written by Dr Yates for KYA320 
    #Method to calculate the coefficients of a polynomial baseline of specified order
    #These coefficients can then be used in apply_baseline() to return the leveled spectra
    #Will exclude the section of the spectra between exclusion min and exclusion max (to remove the spectral features)
    #The inclusion_size is the velocity range on either side of the centre_velocity to be fit (with the exclusion region excluded from the fit)
    def fit_baseline(self, intensity_array, velocity_array, exclusion_min, exclusion_max, centre_velocity, inclusion_size, polynomial_order):
        
        #Check that the intensity and velocity arrays are the same size
        if len(intensity_array) != len(velocity_array):
            print('Length of intensity and velocity arrays did not match while calculating RMS noise \n' + self.termination_message)
            quit()
        
        #Check that the exclusion area is smaller than the inclusion area
        if (centre_velocity - exclusion_min >= inclusion_size) or (exclusion_max - centre_velocity) >= inclusion_size:
            print('Exclusion size must be smaller than inclusion size for fitting the baseline \n' + self.termination_message)
            quit()
        
        #Calculate the velocity bounds for inclusion
        vel_inc_low = centre_velocity - inclusion_size
        vel_inc_up = centre_velocity + inclusion_size

        #Use numpy logical arrays to exclude the sections of the data that are not to be fit (spectral feature, and data outside of inclusion range)
        logic_inclusion_data = np.logical_and(velocity_array > vel_inc_low, velocity_array < vel_inc_up)
        logic_not_excluded = np.logical_or(velocity_array < exclusion_min, velocity_array > exclusion_max)

        #Combine these two logical arrays with an and operation (must be within inclusion bounds, and outside exclusionm bounds)
        logic_selection = np.logical_and(logic_inclusion_data, logic_not_excluded)
                
        #Can use np.polyfit to fit the polynomial, only using the values specified by the logic array
        polynomial_coefficients = np.polyfit(velocity_array[logic_selection], intensity_array[logic_selection], deg=polynomial_order)
        
        return(polynomial_coefficients)
        
        
    #Simple method to remove the baseline from an intensity array
    #Will generate the values of a polynomial from the coefficients provided using np.polinomial.polynomial.polyval
    #Will then subtract them from the intensity array provided
    #polynomial coefficients should be in order of decreasing power
    #TODO: Error catchers
    def remove_baseline(self, intensity_array, velocity_array, polynomial_coefficients):
        
        #First, use polyval to generate an array of polynomial data
        #An annoying feature of .polyval, it requires coefficients to be in reverse order to other np stuff
        #Fip the order of the polynomial coefficients
        flipped_coefficients = np.flip(polynomial_coefficients)

        #Can now calculate the polynomial values
        #Will use the velocity array as 'x' and the coefficients entered (dervied in fit_baseline)
        polynomial_array = np.polynomial.polynomial.polyval(velocity_array, flipped_coefficients)

        #Now subtract these values from the intensity array
        baseline_removed_intensity_array = intensity_array - polynomial_array
        return(baseline_removed_intensity_array)
        #return(polynomial_array) #used to test what the polynomial looks like

        
    #Method to wrtie a .spt file given the parameters of the file.
    #Some values that will be the same for all plots in this project are entered as default values
    #This method is just formatting, details of it would be as long as it is, read in line for further help
    #Will write out a file with the name file_name.spt
    #TODO: Add a label option to add things like rms, baseline properties, transition ect
    def write_file(self, file_name, min_vel, max_vel, data_point_count, min_flux, max_flux, velocity_array, intensity_array, x_label = 'Velocity w.r.t. LSR (km/s)',y_label = 'Flux (Jy)'):
        
        #Create the file
        plot = open(str(file_name) + ".spt", "x")

        #Insert the name and axis label components
        plot.write(" heading \n " + str(file_name) + " \n xlabel \n " + str(x_label) + " \n ylabel \n " + str(y_label) + " \n ")
        #Instert the velocity range on the bottom axis
        plot.write("baxis\n " + str(min_vel) + "\n " + str(max_vel) + "\n")
        #Insert the top axis, which gives information about the data point count
        plot.write("taxis\n " + str(data_point_count) + "\n 1\n")
        #Insert the flux range on the left axis
        plot.write("laxis\n " + str(min_flux) + "\n " + str(max_flux) + "\n")

        #Now enter the number of data points as the first section of the data
        plot.write("data\n " + str(data_point_count) + "\n ")

        #Now insert the data of the spectrum
        for i in range(len(velocity_array)):
            if (velocity_array[i] > min_vel and velocity_array[i] < max_vel): #If the velocity is within the range, write the flux/vel data into the spectra
                plot.write(str(velocity_array[i]) + "    " + str(intensity_array[i]) + " \n ")

        #Now add an end comment and close the file
        plot.write("end")
        plot.close()
    

    #TODO: ADD OPTION TO USE THIS
    #Method to produce a csv data output of the velocity and intensity.
    #This format will be useful for doing composite plot generation, and for doing Fourier transform stuff on
    def write_data(self, file_name, min_vel, max_vel, velocity_array, intensity_array, epoch, x_name = 'vel(km/s)', y_name = 'intensity(Jy)'):
        
        #Create and open the csv file
        #For csv, just have value,value. so only thing other than input will be + ',' + between each value, with '\n' for a new line
        data = open("data_" + str(file_name) + ".csv", "x")

        #Put in the date metadata
        data.write('#' + str(epoch) + '\n')

        #Put in the headers
        #Eah value needs to have quotation marks around it
        data.write('"' + str(x_name) + '"' + ',' + '"' + str(y_name) + '"' + '\n')

        #Go through the data and add it to the csv
        for i in range(len(velocity_array)):
            if (velocity_array[i] > min_vel and velocity_array[i] < max_vel): #If the velocity is within the range, write the flux/vel data into the csv file
                data.write('"' + str(velocity_array[i]) + '"' + ',' + '"' + str(intensity_array[i]) + '"' + '\n')

        data.close()


    #Method to produce a plot of a spctra
    #Will wrtie out a file using write_file
    #velocity array and the spectra arrays hold the data to be plotted (for plots not using off, can just be left as None)
    #plot name is the string that will be used for the name of the file and plot
    #Vel source is the plot the spectra will be centered on, also used to center the inclusion window for baseline and RMS calculations
    #Vel source low and high used to exclude spectral features from RMS and baseline. Should be tightly bound to spectral features
    #Inclusion size is the width either side of vel_source to be considered for RMS and baseline calculations. Polynomial order specifies what order polynomial to use for baseline
    def write_plot(self, velocity_array, spectra_array_on, spectra_array_off, intensity_scale, plot_name, vel_source, vel_source_low, vel_source_high, epoch, format, inclusion_size, polynomial_order):
        #Write out the file to a .spt file, to be used with Simon's script 'splot.pl'

        #Different plotting options are done seperately           

        if self.plot_option == 1:
            print('Plot option 1 selected for ' + plot_name)

            #Need to dervie values that will be used to define the plot parameters
            #Also need to produce the intensity array for on/off, and then remove the baseline

            #Specifiy the plotting window base on the velocity of the spectral feature
            vel_plot_lower = vel_source + self.min_plot_vel
            vel_plot_upper = vel_source + self.max_plot_vel

            intensity_array = spectra_array_on * intensity_scale

            #Calculate the polynomial coefficients for the baseline 
            
            baseline_coefficients = self.fit_baseline(intensity_array = intensity_array, velocity_array=velocity_array, exclusion_min = vel_source_low, exclusion_max = vel_source_high, centre_velocity = vel_source, inclusion_size = inclusion_size, polynomial_order = polynomial_order)

            #Use the coefficients to remove the baseline using remove_baseline(). Skip this step if remove_baseline bool is False
            #rb for removed baseline
            if self.remove_baseline_bool == True:
                intensity_array_rb = self.remove_baseline(intensity_array=intensity_array, velocity_array=velocity_array, polynomial_coefficients=baseline_coefficients)
            else:
                intensity_array_rb = intensity_array

            #,spt form needs data \n 'number of points'
            #Run through the list of velocities, if they are within the range, add to the count
            #Also use this loop to find the maximum value of flux to use in scaling the y axis
            #TODO: Come up with a more efficient way to do this (maybe calculate the number of points base on vel)
            data_point_count = 0
            largest_flux_val = 0
            for i in range(len(velocity_array)):
                if (velocity_array[i]>vel_plot_lower and velocity_array[i]<vel_plot_upper):
                    data_point_count += 1
                    if intensity_array_rb[i] > largest_flux_val:
                        largest_flux_val = intensity_array_rb[i]
            

            #Use the largest flux value found the scale the y axis
            #In the case of very low flux, need to have a minimum scaling size so the plot isn't full of noise
            flux_min = -0.2 * largest_flux_val
            flux_max = 1.2 * largest_flux_val

            if largest_flux_val < 30:
                flux_min = -30
                flux_max = 30

            #Print out the RMS value
            #TODO: Create a way for this to be stored
            #TODO: UPDATE THIS WITH PROPER EXCLUSION
            #Exclusion set to 7 km/s, inclusion set to 25
            RMS = self.calculate_rms(intensity_array = intensity_array_rb, velocity_array = velocity_array, centre_velocity = vel_source, exclusion_max = vel_source_high, exclusion_min = vel_source_low, inclusion_size = inclusion_size)
        
            #Write the file out to the chosen data format
            if format == 'spt':
                self.write_file(file_name = plot_name, min_vel = vel_plot_lower, max_vel = vel_plot_upper, data_point_count = data_point_count, min_flux = flux_min, max_flux = flux_max, velocity_array = velocity_array, intensity_array = intensity_array_rb)
            elif format == 'csv':
                self.write_data(file_name = plot_name, min_vel = vel_plot_lower, max_vel = vel_plot_upper, velocity_array = velocity_array, intensity_array = intensity_array, epoch = epoch)
            else:
                print('Output data format not recognised')
                print(self.termination_message)
                quit()

            #Return the RMS and baseline coefficients so that they can be recorded in the output log
            return(RMS, baseline_coefficients)
            

        if self.plot_option == 2:
            print('Plot option 2 selected for ' + plot_name) 
            
            #Need to dervie values that will be used to define the plot parameters
            #Also need to produce the intensity array for on/off, and then remove the baseline

            #Specifiy the plotting window base on the velocity of the spectral feature
            vel_plot_lower = vel_source + self.min_plot_vel
            vel_plot_upper = vel_source + self.max_plot_vel
            

            #Calculate the quotient between the on and off plots,
            #print(self.spec_data_list[on_plot_position])
            #print(self.spec_data_list[off_plot_position])
            quotient_array = np.zeros(len(spectra_array_on)) #Use on position for size as it must be within list length anyway
            
            for i in range(len(quotient_array)):
                #divide the ith value from the on data by the ith value from the on data, then subtract 1 (quotient formula)
                quotient_array[i] = spectra_array_on[i] / spectra_array_off[i] - 1 

            #Need to convert this ratio of on and off into a flux density.
            #This is done by calibrating the telescope
            intensity_array = quotient_array * intensity_scale

            #Calculate the polynomial coefficients for the baseline 
            
            baseline_coefficients = self.fit_baseline(intensity_array = intensity_array, velocity_array=velocity_array, exclusion_min = vel_source_low, exclusion_max = vel_source_high, centre_velocity = vel_source, inclusion_size = inclusion_size, polynomial_order = polynomial_order)


            #Use the coefficients to remove the baseline using remove_baseline()
            #rb for removed baseline
            if self.remove_baseline_bool == True:
                intensity_array_rb = self.remove_baseline(intensity_array=intensity_array, velocity_array=velocity_array, polynomial_coefficients=baseline_coefficients)
            else:
                intensity_array_rb = intensity_array
            
            #,spt form needs data \n 'number of points'
            #Run through the list of velocities, if they are within the range, add to the count
            #Also use this loop to find the maximum value of flux to use in scaling the y axis
            #TODO: Come up with a more efficient way to do this (maybe calculate the number of points base on vel)
            data_point_count = 0
            largest_flux_val = 0
            for i in range(len(velocity_array)):
                if (velocity_array[i]>vel_plot_lower and velocity_array[i]<vel_plot_upper):
                    data_point_count += 1
                    if intensity_array_rb[i] > largest_flux_val:
                        largest_flux_val = intensity_array_rb[i]
            

            #Use the largest flux value found the scale the y axis
            #In the case of very low flux, need to have a minimum scaling size so the plot isn't full of noise
            flux_min = -0.2 * largest_flux_val
            flux_max = 1.2 * largest_flux_val

            if largest_flux_val < 30:
                flux_min = -30
                flux_max = 30

            #Print out the RMS value
            #TODO: Create a way for this to be stored
            #TODO: UPDATE THIS WITH PROPER EXCLUSION
            #Exclusion set to 7 km/s, inclusion set to 25
            RMS = self.calculate_rms(intensity_array = intensity_array_rb, velocity_array = velocity_array, centre_velocity = vel_source, exclusion_max = vel_source_high, exclusion_min = vel_source_low, inclusion_size = inclusion_size)
            
            print(plot_name)
            #Write the file out to the chosen data format
            if format == 'spt':
                self.write_file(file_name = plot_name, min_vel = vel_plot_lower, max_vel = vel_plot_upper, data_point_count = data_point_count, min_flux = flux_min, max_flux = flux_max, velocity_array = velocity_array, intensity_array = intensity_array_rb)
            elif format == 'csv':
                self.write_data(file_name = plot_name, min_vel = vel_plot_lower, max_vel = vel_plot_upper, velocity_array = velocity_array, intensity_array = intensity_array, epoch = epoch)
            
            #Return the RMS and baseline coefficients so that they can be recorded in the output log
            return(RMS, baseline_coefficients)


    #TODO: Fill this out properly, currently just set for on off
    #TODO: VELOCITY CURRENTLY CALCULATED FROM START OF OBS, WOULD MIDDLE BE BETTER?
    #TODO: Use name to work out IF, chan, and position 
    #   This way, can just pass in a big chunk of spectra and go, do this
    #   Might be strange for on-off if off isn't passed.
    #   Might work if first step is to organise everything and check for stuff like that. That way, can just premtively change from on-off to on
    #   With arrays, need to define entries before values can be assigned, ei can't have test[4] = 1 if test[4] hasn't been created before
    #Want to have it just cycle through all file names and use ths scan number to find name data on the order list. Assume that each on isfollowed by an off in number
    #   If off doesn't exist, will be skpped with warning raised, if on doesn't exist, won't be detected by the sweep through.
    #TODO: If intenisty is done as part of obs and processed with swspec, give them a different name and then pass them in seperately
    #TODO: Add error catcher to allow all plot files to be deleted from within the program, so that new ones can be produced
    def plot_all(self):
        #Two different procedures. One for on / on-off plotting, and another for plotting using every other spectra as baseline (NOT DONE YET)
        
        #Both procedures will have an associated log file, which will store information on the rms, baseline and general parameters used
        log_file = open("log_" + str(self.format) + ".txt", "x")
        log_file.write("Epoch: " + str(self.epoch) + "\n")
        log_file.write("Scan duration (including source change time): " + str(self.observation_length) + " seconds\n")
        log_file.write("Inclusion size: " + str(self.velocity_inclusion_size) + " km/s\n")
        log_file.write("Baseline polynomial order: " + str(self.baseline_polynomial_order) + "\n\n")

        #Go through the full list of file names, using the information within to determine which plots need to be made
        if self.plot_option == 1 or 2:
            #TODO: Add a warning with the on-off to print out that spextra are expected in on-off alternating order, more specifically, each on should have a following off, but multiple offs is fine (ie missing scan)
            #TODO: Implament a way of naming off scans based on preceding scan, like there was with old format
            for i in range(len(self.file_name_list)):
                
                #Get the file name and path out of the lists
                file_name = self.file_name_list[i]
                file_path = self.path_list[i]

                #Need to break apart the string to find the relevant pieces of information, such as IF, channel and scan number
                #In this step, the expected file name format is 'obs name'_IF:'IF'_ch:'channel_number'_st:'station'_'scan_number'_'number_of_data_points'pt_'integration_time's_'processing_format'_swspec.bin
                split_name = re.split('[_:_:_:_____]', file_name)

                #To determine if the file needs to be plotted, extract the scan number
                #This can be used to find the type of observation, and determine if the scan is off (which means it won't be plotted) 
                scan_number = split_name[7] 

                #Need to remove the leading zeroes from the scan number
                scan_number.lstrip('0')

                #Now convert it into and int, and subtract 1 so that it can be used to find the corresponding obs order value from the list
                scan_number_int = int(scan_number) - 1 #Subtract 1 because python starts at 0, while the first scan is scan 1

                #To find the corresponding observation information, take the modulus with the lnegth of the observation order list
                #   This way, if there is a repeat in the schedule, the list doesn't need to be twice as long
                #   For obs with no repeat, the list will be longer than the scan number by 1, so mod will always return the input unchanged
                scan_information_position = np.remainder(scan_number_int, len(self.observation_order))

                #Extract the information for this scan from the observation order list
                scan_tag = self.observation_order[scan_information_position]

                #If the scan tag is 'off', don't need to produce the spectra
                #   For trying to look at off spectra, can just set obs_order to ['test'] and pass in single spectra in question
                if scan_tag != 'off':
                    #Now that it is know that the scan needs to be plotted, the other iformation can be extracted from the name
                    obs_name = split_name[0]
                    IF = split_name[2] #Skipped indice here is for filler_text that is in positions 1,3,5
                    channel_number = int(split_name[4])
                    observatory = split_name[6]
                    num_channels = int(re.split('[p]', split_name[8])[0]) #split at the p in pt, and then just take the first part which is the number 
                    integration_time = int(re.split('[s]', split_name[9])[0]) #Same as above for s after time 
                    processing_format = split_name[10]

                    #If plot option 1 is selected, then there is no need for an off spectra to be extracted
                    #To facilitate this, set off_spectra to None, so that this body of code doesn't need to be repeated or awkwardly defined
                    off_spectra = None

                    #TODO: Add error catcher if off scan picked is going to be insuitable (dif length)
                    #TODO: Add an efficiency optimisation to check if the if, channel and channel numbers are the same as the last scan
                    #   In this case, velocity array will be the same and it doesn't need to be recalculated


                    #Derive parameters that will be used in the creation of the velocity array and later in the spectra data extraction
                    #TODO: Work out exactly what these variables are, so as to better document the code. Currently this is ripped
                        #from Guifre's code
                    #Can incorperate factor of two into both numchannel lines, as it is just due to twice as many channels as there are data points
                    #TODO: Currently bandwidth cannot be changed for 12m, but if in future it can, this can be updated so that the bandwidth can vary between IFs, rather than being fix like it is here
                    
                    #The following five lines of code were based on Guifre's code
                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    SR = int(2 * self.bandwidth) 
                    df = float(SR / num_channels)
                    Nfft = int((num_channels / 2) + 1)
                    jf = np.arange(Nfft, dtype=int)
                    ff = df * jf
                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                    #Use the If and channel number to determine what frequency has been recorded
                    #To get the information out the the tuple of tuples, convert the letter IF to a number
                    #To convert a lower casse letter to it's place number in the alphabet, use the following code
                    #-1 as python starts from 0, so a will be 0, b 1 ect
                    IF_number = ord(IF) - 96 - 1

                    #Check that this IF and and the channel number have a frequency value entered for them. If they don't print an error and quit
                    #IF number starts at 0, channel starts at 1. So IF is >=, chan is > (chan subtraction occurs later)
                    if (IF_number >= len(self.IF_chan_frequencies_MHz)) or channel_number > len(self.IF_chan_frequencies_MHz[IF_number]):
                        print('Error trying to determine frequency associated with file: ' + str(file_name))
                        print('Either the IF number: ' + str(IF_number) + ' or channel number: ' + str(channel_number) + ' could not be found in the frequency 2D array.')
                        print('IF len: ' + str(len(self.IF_chan_frequencies_MHz)) + ' Chan_len: ' + str(len(self.IF_chan_frequencies_MHz[IF_number])))
                        print(self.termination_message)
                        quit()
                    
                    #Now, extract the frequency of the maser from the tuple of tuples
                    #The structure of this 2D array might be better served as a dictionary (so IF letters can be entered)
                    #   but this for makes it easier to adjust for slightly different if/channel structures
                    #Multiply by 1e6 to convert from MHz to Hz
                    maser_freq = self.IF_chan_frequencies_MHz[IF_number][channel_number-1] * 1e6 #-1 because python starts at 0

                    centre_freq = maser_freq
                    #Because of an issue with the early proc files (at least fm025 and earlier), need to adjust the frequency taken from ff
                    if self.frequency_rounding == True:
                        #This would be a simple case of rounding, but Gabor decied to round down for the c and d IFs
                        if IF == 'c' or IF == 'd':
                            centre_freq = math.floor(maser_freq / 1e6) * 1e6 #Divide by 1e6 so that the floor function can be used to round down, then bring back to Hz
                        else:
                            centre_freq = round(maser_freq / 1e6) * 1e6

                    #Convert full frequency to velocity
                    #As ff is centred around 16 MHz, to move that to be centred on centre freq, need to subtract 16 and add centre freq
                    
                    #Little fudge factor used to get velocity to line up if need be
                    fudge_factor = 0
                    
                    ff_cen = ff + (centre_freq - 16.0*1e6 + fudge_factor)

                    #Now convert the frequency to velocity using the frequency of the transition
                    #This velocty does not factor in the Earth's motion realtive to LSR
                    vel_freq_array = 2.998*1e8 * (1 - ff_cen / maser_freq) / 1000 #Use doppler eq, divide by 1000 for Km/s

                    #To calculate the velocity correction, need to extract the coordnates of the observatory used and the source information, as well as calculate the time of the scan
                    #Check to make sure the observatory name is contained in the dictionary
                    if self.observatories.get(observatory) == False: #If the observatory name is not conatined in the dictionary
                        print('Observatory id: ' + observatory + ' not found in observatories dictionary for file: ' + file_name)
                        print(self.termination_message)
                        quit()
                    
                    obs_coords = self.observatories[observatory] #Tuple in the form (lat, long, alt)
                    lat = obs_coords[0]
                    lon = obs_coords[1]
                    alt = obs_coords[2]


                    #Extract the source information using the scan tag
                    if self.source_params.get(scan_tag) == None:
                        print('Name entered in -names is not contained in the dictionary of sources for tag: ' + scan_tag + 'from file: ' + file_name)
                        print(self.termination_message)
                        quit()
                    
                    source_parameters = self.source_params[scan_tag]

                    #Extract the information from the source parameters                   
                    #The format is (vel_centre, vel_low, vel_high, RA, Dec, full_name)
                    vel_centre = source_parameters[0]
                    vel_low = source_parameters[1]
                    vel_high = source_parameters[2]
                    RA = source_parameters[3]
                    dec = source_parameters[4]
                    full_name = source_parameters[5]


                    #Calculate the time of the scan
                    #To do this, multiply the scan number -1 by the length of a single scan, and add that to the starting epoch
                    #(Scan number -1 because the the first scan starts at the starting time, and the second one scan later ect)
                    #Also add on a centre time, which can be used to shift the time at which the velocity is calculated during each observation. Currently 0
                    centre_time = int(self.observation_length / 2) #Use half the observation length to place the time at the middle of the observation.
                    epoch_numpy = np.datetime64(self.epoch) + np.timedelta64(scan_number_int * self.observation_length + centre_time, 's')
                    

                    #Calculate the appropriate velocity correction and apply it to the velocity array
                    velocity_correction = self.velocity_correction(lat = lat, long = lon, alt = alt, epoch = epoch_numpy, RA = RA, dec=dec)
                    print(velocity_correction)
                    #Apply the velocity correction
                    corrected_velocity_array = vel_freq_array - velocity_correction

                    #The following lines of code are based upon code written by Guifre, without much change
                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    #Extract the spectra information from the .bin file
                    #The following three lines of code were based on Guifre's code
                    file_temp_on = open(f"{file_path}{file_name}", "rb") #Open the .bin file
                    spectra_array = np.fromfile(file=file_temp_on, dtype="float32", count=Nfft) #Extract the data from the file
                    file_temp_on.close #Close the file
                    
                    #If an on-off plot is to be produced (plot option 2), then extract the scan that comes next in the list of files
                    if self.plot_option == 2:
                        file_name_off = self.file_name_list[i+1]
                        file_path_off = self.path_list[i+1]
                        file_temp_off = open(f"{file_path_off}{file_name_off}", "rb") #Open the .bin file
                        off_spectra = np.fromfile(file=file_temp_off, dtype="float32", count=Nfft) #Extract the data from the file
                        file_temp_off.close #Close the file
                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


                    #Use the information above the construct the name for the plot
                    plot_name = str(full_name) + "_" + str(obs_name) + "_IF:" + str(IF) + "_" + str(maser_freq / 1e6) + "MHz_scan:" + str(scan_number) + "_" + str(observatory)

                    #If the plot is on-off, add that to the end of the name
                    if self.plot_option == 2:
                        plot_name = plot_name + str("_on_off")

                    #Can now use write_plot to turn this information into a spectra
                    #Write plot will return the rms and polynomial baseline coefficients
                    RMS, baseline_coefficients = self.write_plot(velocity_array = corrected_velocity_array, spectra_array_on = spectra_array, spectra_array_off = off_spectra, intensity_scale = self.intensity_scale, plot_name = plot_name, vel_source = vel_centre, vel_source_low = vel_low, vel_source_high = vel_high, epoch = epoch_numpy, format = self.format, inclusion_size = self.velocity_inclusion_size, polynomial_order = self.baseline_polynomial_order)

                    #Now the plot will have been made. Only thing left is to add an entry into the log file with the full information on the plot
                    log_file.write(plot_name + "\n")
                    log_file.write("    Time of observation: " + str(epoch_numpy) + "\n")
                    log_file.write("    Integration time: " + str(integration_time) + "s\n")
                    log_file.write("    Spectrometer format: " + str(processing_format) + "\n")
                    log_file.write("    Number of channels: " + str(num_channels) + "\n")
                    log_file.write("    RMS: " + str(RMS) + "\n")
                    log_file.write("    Baseline polynomial coefficients: " + str(baseline_coefficients) + "\n")
                    log_file.write("\n") #leave a gap between entries in the log file to improve readability


plot_maser = maser_plot()
plot_maser.plot_all()