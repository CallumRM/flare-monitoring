from astropy.time.formats import TimeStardate
import numpy as np
import os, argparse, re
import errno
import sys
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
#TODO: Perl output, Do on-off or baseline for multiple inputs. Combine scans, should be easy, just need to convert vel before merging
#Will need to take in the start time of the observing epcoh in numpy YYYY-MM-DD'T'HH:MM:SS form, and then have scan length, use numpy datetime, and then use timedelta to add extra
#This way all data can be velocity corrected before before being plotted.
#TODO: Need to set intensity, will need to fit guassian and then use that to find conversion factor
    #Or use other intensity calibration that Simmon talked about
#TODO: Set up flexiblity for different observation patterns. Use naming to identify intenisty calibration scans.
#TODO: Make a script that can read a schedule file and produce the inputs required from it (list of source names with on off stuff)
#TODO: Make option to print out the baseline fit. Have it name it properly. Can also have naming change if baseline removed
#TODO: Make a log file that will hold input params and values like rms
#TODO Add overide variable that can be used to have easy 'profiles' for cd and ef under standard fm obs.
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
        #Bandwidth as an int in MHz
        parser.add_argument(
            "-bw",
            dest="bandwidth",
            help="observing bandwidth in MHz", 
            type=float, 
            default=32 #default set to 32 for the AuScope 12m antennas
        )

        #Number of channels
        #For formatting, see bandwidth above. Use -nc 'number of channels'
        #Must be an integer (there are an integer number of channels)
        #TODO: SET THIS UP FOR 12 and 6 INDIVIDUALLY (Could set it up for each channel-if combo individually)
        parser.add_argument(
            "-nc",
            dest="num_channels",
            help="The number of channels recorded",
            type=int,
            default=64000 #default set to 64000 which is the number used for 6.7 GHz obs, for 12.2 use 32000
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

        #Location of observatory
        #Either to initials of the included observaatories, or the lat, lon and altitude
        #Give coordinates in degrees, altitude in metres
        #Default set to Hobart 12 m
        #TODO: Allow DMS form to be given
        #Included observatories are:
        #   hb: Hobart 12m
        #   ke: Katherine 12m
        parser.add_argument(
            "-obs",
            dest="observatory",
            type=str,
            default='hb'
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

        #Observation order
        #The order of the sources observed, using names which are found in the dictionary below
        #TODO: Make a nice external reference to this dictionary.
        #The names of the observations in a list (list of strings)
        #Default set for fm observations of 6 sources in on off format
        #TODO Add something like group name, so that fm*** can be added to the end of each name.
        parser.add_argument(
            "-order",
            dest="observation_order",
            nargs="+",
            default=['G328.237', 'off', 'G328.809', 'off', 'G351.417', 'off', 'G323.740', 'off', 'G318.948', 'off', 'G9.621', 'off'] 
        )
        
        #Show spectrum
        #Boolean value. If true, the spectrum will be displayed.
        #If false, the spectrum will be saved as a .svg
        #Default false for ease of use with automation
        #TODO: Separate this and save plot, so that both can be done simultaneously
        #Currently does nothing
        parser.add_argument(
            "-show_spec",
            dest="show_spec",
            default=False
        )


        #The following variables are used for managing the number of spectra to be plotted
        #IFs
        #The IFs to be processed. List of chars
        parser.add_argument(
            "-ifs",
            dest="IFs",
            help="The IFs to be processed. List of chars",
            nargs="+",
            default=['c', 'd', 'e', 'f'] #Default set to all IFs used in fm observations
        )

        #Channels
        #The channels to be processed. List of ints
        parser.add_argument(
            "-chan",
            dest="channels",
            help="The channels to be processed. List of ints",
            nargs="+",
            default=[1, 2, 3, 4] #Default set to all four channels used in fm observations
        )
        
        #The following inputs are for the different plotting options
        #If just one spectra is to be plotted, the number position in the list must be sepcified
        #If an on-off is to be plotted, the position of the on and off must be specified
        #If a plot of one spectra using the others as a reference is to be produce, the psition of
            #the one to be plotted must be specified

        #Plotting option argument
        #-po "number". 1 will just plot a single spectra, 2 will plot on-off, 3 will plot on-off
            #with the average of other spectra used as off
        parser.add_argument(
            "-po",
            dest="plot_option",
            help="The plotting option",
            type=int,
            default=1 #default to 1, to just plot a spectra without any quotient
        )

        #Maser Frequency
        #The frequency in MHz of the maser transition being observed
        parser.add_argument(
            "-freq",
            dest="maser_freq",
            help="The frequency of the maser, in MHz",
            type=float,
            default=6668.5192 #The 6.7 GHz maser as the default, as this is the most commonly observed transition
        )

        #Single spectra position
        #Will be used with either option 1 or 3
        #Position in the list, starting from 1 for simplicity, as scans start with 1
        parser.add_argument(
            "-s_pos",
            dest="single_spec_position",
            help="The position in the list of spectra for a signle spectrum plot (op 1 or 3), starting from 1",
            type=int,
            default=1 #default to 1, as at least 1 spectrum is alwasy passed in
        )

        #On-off spectra positions
        #These are both defaulted to 1 to avoid  crashes, but if both are left as 1, a trivial plot will be produced
        #TODO: Add warning in the on-off plot section if both are set to same value
        #Use as -on_off_pos "on_list_position" "off_list_position" (accepts two numbers separated by a space)
        parser.add_argument(
            "-on_off_pos",
            dest="on_off_positions",
            nargs="+",
            type=int,
            default=(1,1) #Default to 1 and 1, as there will always be at least one spectrum
        )

        #End of input arguments
        #Now extract variables out of the parser for simplicity
        #TODO: Have a print out of the entered parameters
        args = parser.parse_args() #args now is a list of all of the input variables
        self.bandwidth = int(args.bandwidth) * 1e6 #Convert bandwidth to Hz
        self.num_channels = int(args.num_channels)
        self.spec_list = args.spec_list
        self.intensity_scale = float(args.intensity_scale)
        self.epoch = str(args.epoch)
        self.observatory = str(args.observatory)
        self.observation_length = int(args.observation_length)
        self.observation_order = args.observation_order
        self.plot_option = int(args.plot_option)
        self.maser_freq = float(args.maser_freq * 1e6) #Convert to Hz
        self.centre_freq = int(args.maser_freq * 1e6) #integer frequency for window centering.
        self.show_spec = bool(args.show_spec)
        self.min_plot_vel = int(args.velocity_range[0]) #extract the minimum and maximum velocity to plot. Need int, as it reads as string
        self.max_plot_vel = int(args.velocity_range[1])
        self.single_spec_position = int(args.single_spec_position) - 1 #Subract 1 as python lists start at 0
        self.on_position = int(args.on_off_positions[0]) - 1 #Get the on (first) value from the tuple, subract 1 as ptyhon lists start at 0
        self.off_position = int(args.on_off_positions[1]) - 1 #Get the off (second) value from the tuple, subract 1 as ptyhon lists start at 0
        self.termination_message = 'Script terminating' #Message to be used when the script encounters a potential error and quits. Can be changed later to something more interesting.

        #Extract the observatory location data
        #If a name has been provided, use the associated coords
        #Location data from (J. Lovell et.al 2012, The AuScope geodetic VLBI array, http://www.phys.utas.edu.au/physics/mt_pleasant_observatory.html and https://www.exploroz.com/places/103554/sa+ceduna-radio-astronomy-observatory)
        #Converted into decimal degrees
        
        #Dictionary of the observatories around Australia that are being used or may be used for this project
        #Each entry is a tuple of the form (lat, long, alt), with the lat and long in degrees, and the alt in metres
        self.observatories = {
            'hb' : ( -42.80557, 147.4381, 40.967), #Hobart 12m
            'ke' : (-14.37546, 132.15237, 189.262), #Katherine 12m
            'yg' : (-29.04715, 115.34563, 248.236), #Yarragadee 12m
            'ho' : (-42.805, 147.439, 43), #Hobart 26m (info not specific, might be worth finding good source)
            'cd' : (-31.86883, 133.809799, 154), #Ceduna 30m
            'parkes' : (-32.99839, 148.26352, 415), #Parkes 64m (For vel calibration tests using ATNF online tool, sourced from Parkes user guide)
        }


        #If the observatory name entered is not within the dictionary (.get() returns None if it cannot find it) then quit
        #Otherwise, extract the coord data
        if self.observatories.get(self.observatory) == None:
            print('Observatory data entered using -obs does not match any of the observatory labels')
            print(self.termination_message)
            quit()
        else:
            coords = self.observatories[self.observatory] #Tuple in the form (lat, long, alt)
            self.lat = coords[0]
            self.long = coords[1]
            self.alt = coords[2]

        #TODO: Add error catchers for all of the other values (ie plot option, spec_pos out of bounds)
        #If the spec_list is empty, there is no data to be plotted, and the program quits
        if len(self.spec_list) < 1:
            print('Must have at least one spectrum passed using -spec')
            print(self.termination_message)
            quit()


        #TODO: It might be useful to have some specifications as to the width of the source here. Similar thing to script from in sem work
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
        self.source_params = {
            'off' : (0, -5, 5, '0h0m0s', '0d0m0s', 'none'),
            'G328.237' : (-44.7, -47.0, -31.5, '15h57m58.28s', '-53d59m22.7s', 'G328.237-0.547'),
            'G328.809' : (-44.5, -47.5, -42.0, '15h55m48.70s', '-52d43m05.5s', 'G328.809+0.633'),
            'G351.417' : (-10.4, -12.0, -4.5, '17h20m53.37s', '-35d47m01.2s', 'G351.417+0.645'),
            'G323.740' : (-50.5, -59.0, -42.0, '15h31m45.45s', '-56d30m50.1s', 'G323.740-0.263'),
            'G318.948' : (-34.6, -39.0, -31.0, '15h00m55.40s', '-58d58m52.1s', 'G318.948-0.196'),
            'G9.621' : (1.3, -4.8, 8.9, '18h06m14.67s', '-20d31m32.4s', 'G9.621+0.196') #Note, 8.9 includes small feature, might make baseline fitting poor due to large gap
        }


        #Use the above dictionary to extract the source parameters
        #If there are more input files than sources specified, the list of sources will be looped.
            #i.e. 24 obs, 12 sources 'a', 'off', 'b', 'off' ..., then observation 13 will be 'a' again, and 14 'off'
        #For off scans, the coordinates from the previous scan will be used, but with off appended to the name

        #Define the lists that will hold the observation parameter data in order
        #Observations names is just a list of strings
        #The coordinates will be a list of tuples, where each tuple corresponds to the (RA, Dec) of a single source
        #The velocities will be a tuple of the form (vel_cen, vel_low, vel_high) for a single source
        self.observation_names = []
        self.observation_coordinates = []
        self.source_velocities = []


        for i in range(len(self.spec_list)):
            #Need to ensure that the reference position is within the names list.
            #If it is greater than the length of the list, take the modulus with the list length
            if i < len(self.observation_order):
                scan_num = i
            else:
                scan_num = np.remainder(i, len(self.observation_order))
            
            #Check if the name in observation_order[scan_num] is located in the dictionary.
            #Use .get(), will return number value for valid entries, or None for invalid
            if self.source_params.get(self.observation_order[scan_num]) == None:
                print('Name entered in -names is not contained in the dictionary of sources.')
                print(self.termination_message)
                quit()
            
            #Now the name is either a source or off.
            #If it is off, then the information from the previous scan will be used
                #The name will have a 'off' attached to it
            #If the off scan is in the first position, then 0 values will be used

            #Reset the 'off' tag to an empty string
            name_addition = ''
            if self.observation_order[scan_num] == 'off':
                #If the first scan entered is off, then set values to 0 and 'off'
                if scan_num == 0:
                    print('Off scan entered in first scan position. Assigning 0 vel, 0,0 coords, and \'off\' name.')
                #If it is not the first scan, set the scan number to 1 less
                else:
                    scan_num = scan_num - 1
                    name_addition = '_off'
            
            #Extract the parameters from the dictionary and append them to the corresponding lists.
            params_temp = self.source_params[self.observation_order[scan_num]]

            #Params have the format (vel, RA, Dec, name)
            #Append each parameter to it's corresponding list.
            #Velocities are entered in (vel_centre, vel_low, vel_high) tuples
            self.source_velocities.append((params_temp[0], params_temp[1], params_temp[2]))
            #coords are entered in (RA, Dec) tuples
            self.observation_coordinates.append((params_temp[3], params_temp[4]))
            #Names are entered as strings sourced from the dictionary. The name addition will add the '_off' to the end if needed 
            self.observation_names.append(str(params_temp[5]) + str(name_addition))
        #Now the lists of coordinates, velocitie and names have been created for a flexible schedule order.


        #Extract the data from the spectra .bin files
        #Check if the input .bin spectra files exist, and create a list of file names and a list of the paths to these spectra
        #The splitting of name and path may not be necessary, but I don't know enough about opening files to make a call.
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

        #TODO: Get spectra data out for each spectra, as frequencies will be different
        #Now calculate some variables and lists to be used for the plot
        #TODO: Work out exactly what these variables are, so as to better document the code. Currently this is ripped
            #from Guifre's code
        self.SR = int(2 * self.bandwidth) #Two bits
        self.df = float(self.SR / self.num_channels) #Differential freq unit
        self.Nfft = int((self.num_channels / 2) + 1)
        self.jf = np.arange(self.Nfft, dtype=int)
        self.ff = self.df * self.jf

        #Convert full frequency to velocity
        #As ff is centred around 16 MHz, to move that to be centred on centre freq, need to subtract 16 and add centre freq
        #Little fudge factor used to get velocity to line up
        fudge_factor = -0.0355*1e6
        #TODO: work out why this is necessary
        self.ff_cen = self.ff + (self.centre_freq - 16.0*1e6 + fudge_factor)

        #Now convert the frequency to velocity using the frequency of the transition
        #This velocty does not factor in the Earth's motion realtive to LSR
        #That will be done for each plot individually, as the correction will change over the course of the observation
        self.general_vel = 2.998*1e8 * (1 - self.ff_cen / self.maser_freq) / 1000 #Use doppler eq, divide by 1000 for Km/s
        

        #Make a list of the spectra data extracted from the .bin files
        self.spec_data_list = []
        for i in range(len(self.spec_list)): #For each file name passed in
            file_temp = open(f"{self.path_list[i]}{self.file_name_list[i]}", "rb") #Open the .bin file
            temp_data = np.fromfile(file=file_temp, dtype="float32", count=self.Nfft) #Extract the data from the file
            self.spec_data_list.append(temp_data) #Put array of data in the list
            file_temp.close #Close the file
    

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
    def velocity_correction(self, lat, long, alt, epoch, RA, Dec):
        #Convert the epoch, sky coordinates and observation location into the fpormats required for bary_vel-cor
        epoch_astro = astropy.time.Time(epoch)
        sky_coords = astropy.coordinates.SkyCoord(RA, Dec, frame = 'icrs')
        obs_location = astropy.coordinates.EarthLocation.from_geodetic(long, lat, height=alt)

        #Use these values to calculate the barycentric velocity correction using the github sourced function bary_vel_cor
        bary_vel_output = self.bary_vel_cor(time=epoch_astro, skycoord = sky_coords, location = obs_location)

        #The formatting of this output is not just a float value, and it needs to be converted to one for later
        #First, get the velocity string out of the output
        bary_vel_string, space_filler , space_filler , space_filler = re.split('[ ]', str(bary_vel_output))

        #Now convert the string to a float, and change it to recession velocity (negative)
        bary_velocity = -float(bary_vel_string)
        #print(bary_velocity)


        #To convert to LSR, define the direction of the LSR, then take the dot product of the unit vectors
        #This result can be multiplied by 20 km/s to get LSR velocity conversion factor
        #TODO: Find a source for these values, currently just using the script that Jim wrote, given to me by Simon
        lsr_sky_coords = astropy.coordinates.SkyCoord(270*units.deg, 30*units.deg, frame = 'icrs')

        #Skycoords.separation returns the separation between two skycoord objects in degrees
        dms_separation = sky_coords.separation(lsr_sky_coords)
        
        #Convert this to degrees, then radians
        #Use re.split to get the degree, min, sec components, then add them, with the mins and seconds converted to deg, then convert to rad
        d_sep, m_sep, s_sep, space_filler = re.split('[dms]', str(dms_separation))
        separation = (float(d_sep) + float(m_sep)/60 + float(s_sep)/(60*60)) * np.pi / 180
        
        #Now the separation can be used to calculate the dot product between the two, and this can be multiplied by 20 km/s to give LSR vel correction
        lsr_correction = -20 * np.cos(separation)

        #Combine the two velocities to get the LSR corrected velocity
        LSR_velocity = bary_velocity +lsr_correction
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
            print('Length of intensity and velocity arrays did not match while calculating baseline noise \n' + self.termination_message)
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


    #Method to produce a plot of a spctra
    #If the plot position are left unspecified, they will be set to the input spec position.
    #For the plot_all method, they will be set as required.
    #Will wrtie out a file
    #TODO: Create system for specifying the inclusion size and polynomial order properly, currently hard set
    def write_plot(self, single_plot_position = None, on_plot_position = None, off_plot_position = None, velocity_correction = 0, inclusion_size = 20, polynomial_order = 7):
        #Write out the file to a .spt file, to be used with Simon's script 'splot.pl'
        
        #If input position were not specified, set them to the initialised values
        if single_plot_position == None:
            single_plot_position = self.single_spec_position
        
        if on_plot_position == None:
            on_plot_position = self.on_position
        
        if off_plot_position == None:
          off_plot_position = self.off_position

        #Now write the plot
        #Different plotting options are done seperately           

        if self.plot_option == 2:
            print("Plot option 2 selected, on off plot: " + str(on_plot_position) + " " +str(off_plot_position))
            print("Source name: " + self.observation_names[on_plot_position])
            
            #Need to dervie values that will be used to define the plot parameters
            #Also need to produce the intensity array for on/off, and then remove the baseline
            plot_name = str(self.observation_names[on_plot_position]) + str('_on_off')


            #Extract the velocity information
            #Velocities stored as tuple of the form (vel_centre, vel_low, vel_high)
            #The low and high velocities will be used to fit the baseline and calculate the RMS, they should be at the bounds of the spectra features
            vel_source = self.source_velocities[on_plot_position][0]
            vel_source_low = self.source_velocities[on_plot_position][1]
            vel_source_high = self.source_velocities[on_plot_position][2]

            #Specifiy the plotting window base on the velocity of the spectral feature
            vel_plot_lower = vel_source + self.min_plot_vel
            vel_plot_upper = vel_source + self.max_plot_vel
            
            #Now subtract the velocity correction from the velocity values, which will place the maser at the source velocity
            adjusted_vel = self.general_vel - velocity_correction

            #Temporary code for plotting with frequency as the bottom axis
            #adjusted_vel = self.ff_cen
            #vel_plot_lower = 6667*10**6
            #vel_plot_upper = 6670*10**6



            #Calculate the quotient between the on and off plots,
            #print(self.spec_data_list[on_plot_position])
            #print(self.spec_data_list[off_plot_position])
            quotient_array = np.zeros(len(self.spec_data_list[on_plot_position])) #Use on position for size as it must be within list length anyway
            #print(quotient_array)
            for i in range(len(quotient_array)):
                #divide the ith value from the on data by the ith value from the on data, then subtract 1 (quotient formula)
                quotient_array[i] = self.spec_data_list[on_plot_position][i] / self.spec_data_list[off_plot_position][i] - 1 

            #Need to convert this ratio of on and off into a flux density.
            #This is done by calibrating the telescope
            intensity_array = quotient_array * self.intensity_scale

            #Calculate the polynomial coefficients for the baseline 
            
            baseline_coefficients = self.fit_baseline(intensity_array = intensity_array, velocity_array=adjusted_vel, exclusion_min = vel_source_low, exclusion_max = vel_source_high, centre_velocity = vel_source, inclusion_size = inclusion_size, polynomial_order = polynomial_order)


            #Use the coefficients to remove the baseline using remove_baseline()
            #rb for removed baseline
            intensity_array_rb = self.remove_baseline(intensity_array=intensity_array, velocity_array=adjusted_vel, polynomial_coefficients=baseline_coefficients)
            #intensity_array_rb=intensity_array

            #,spt form needs data \n 'number of points'
            #Run through the list of velocities, if they are within the range, add to the count
            #Also use this loop to find the maximum value of flux to use in scaling the y axis
            #TODO: Come up with a more efficient way to do this (maybe calculate the number of points base on vel)
            data_point_count = 0
            largest_flux_val = 0
            for i in range(len(adjusted_vel)):
                if (adjusted_vel[i]>vel_plot_lower and adjusted_vel[i]<vel_plot_upper):
                    data_point_count += 1
                    if intensity_array_rb[i] > largest_flux_val:
                        largest_flux_val = intensity_array_rb[i]
            

            #Use the largest flux value found the scale the y axis
            #In the case of very low flux, need to have a minimum scaling size so the plot isn't full of noise
            if largest_flux_val < 30:
                largest_flux_val = 30
            flux_min = -0.2 * largest_flux_val
            flux_max = 1.2 * largest_flux_val

            #Print out the RMS value
            #TODO: Create a way for this to be stored
            #TODO: UPDATE THIS WITH PROPER EXCLUSION
            #Exclusion set to 7 km/s, inclusion set to 25
            #RMS = self.calculate_rms(intensity_array = intensity_array_rb, velocity_array = adjusted_vel, centre_velocity = vel_source, exclusion_max = vel_source_high, exclusion_min = vel_source_low, inclusion_size = inclusion_size)
            #print(RMS)
        
            #Use write_file to write the .spt file using the values calculated above
            self.write_file(file_name = plot_name, min_vel = vel_plot_lower, max_vel = vel_plot_upper, data_point_count = data_point_count, min_flux = flux_min, max_flux = flux_max, velocity_array = adjusted_vel, intensity_array = intensity_array_rb)
    
    #TODO: Fill this out properly, currently just set for on off
    #TODO: Need to be able to combine plots (double feature, position + (len/2 + pos))
    #   Probably better achieved later using another program
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
        #Two different procedures for on / on-off plotting, and for plotting using every other spectra as baseline

        #Go through the full list of file names, using the information within to determine which plots need to be made
        if self.plot_option == 1 or 2:
            for i in range(len(self.file_name_list)):
                file_name = self.file_name_list[i]

                #Need to break apart the string to find the relevant pieces of information, such as IF, channel and scan number
                #In this step, the expected file name format is 'obs name'_IF:'IF'_ch:'channel_number'_st:'station'_'scan_number'_'number_of_data_points'_'integration_time'_'processing_format'_swspec.bin
                split_name = re.split('[_:_:_:_____]', file_name)

                #Extract the specific values from the 
                obs_name = split_name[0]
                IF = split_name[2] #Skipped indice here is for filler_text that is in positions 1,3,5
                channel_number = split_name[4]
                station = split_name[6] 
                scan_number = split_name[7] 
                num_data_points = split_name[8], 
                integration_time = split_name[9], 
                processing_format = split_name[10]

                #Now use this information to determine what values need to be calculated to allow plotting
                #The first thing to do is to use the scan number do determine what this scan is, using the observation order

                #Need to remove the leading zeroes from the scan number
                scan_number.lstrip('0')

                #Now convert it into and int, and subtract 1 so that it can be used to find the corresponding obs order value from the list
                scan_number_int = int(scan_number) - 1 #Subtract 1 because python starts at 0, while the first scan is scan 1

                scan_tag = self.observation_order[scan_number_int]

                #If the scan tag is 'off', don't need to produce the spectra
                #   For trying to look at off spectra, can just set obs_order to ['test'] and pass in single spectra in question
                if scan_tag != 'off':
                    #Need 
                    #Informaiton needed by plot (currently from pos) name, velocity stuff, on and off spectra
                    test=True




                #TODO: GOING TO NEED TO DELETESTUFF IN __init__ and instead just define those things up there.
                #Then just get info out of dictionaries here. Maybe slightly slower, but with computing power available will only add ~order minutes to processing time


    def plot_all_old(self, start_position = 1):
        #Will take in a start position, and then produce .spt for each on off pair
        #Will add time of previous scans to the epoch to properly correct veloity
        self.plot_option = 2
        #This is the starting time.
        epoch_numpy = np.datetime64(self.epoch)

        #Convert from scan number to python list indicee
        i = start_position - 1
        #Use -1, as the off scans will be i + 1, and the last scan will be an off scan.
        while i+1 < len(self.spec_data_list):
            #print(i)
            #Calculate the velocity correction
            RA_temp = self.observation_coordinates[i][0]
            Dec_temp = self.observation_coordinates[i][1]
            #print(RA_temp)
            #print(Dec_temp)

            vel_cor_temp = self.velocity_correction(self.lat, self.long, self.alt, epoch_numpy, RA_temp, Dec_temp)

            #Write the .spt plot
            self.write_plot(on_plot_position = i, off_plot_position = i+1, velocity_correction = vel_cor_temp)

            #Set the epoch to two scan times forwards (as an on and off have occured between previous and current)
            on_off_time_delta = np.timedelta64(2 * self.observation_length, 's')
            epoch_numpy = epoch_numpy + on_off_time_delta

            #Interate i to th next on scan, which is after the follow off scan, so +2
            i = i + 2


plot_maser = maser_plot()
plot_maser.plot_all_old()


#TODO: Proccess entire observation
#Also want option to use other scans to form baseline
#Need it to be flexible
#Lots of the scans are cotemporal, so there is probably a smart way to avoid excesive redefinition of datetime
