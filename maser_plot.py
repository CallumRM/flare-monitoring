import numpy as np
import os
import argparse
import re
import errno
import math
from astropy import units
import astropy.time
import astropy.coordinates
import scipy.fft
import scipy.optimize
#from astropy.coordinates.builtin_frames.ecliptic import BaseEclipticFrame
#from astropy.time.formats import TimeStardate
#import PyAstronomy.pyasl
#from numpy.lib.polynomial import poly
#from erfa.core import cr


#Author: Callum Macdonald
#Direct any questions to either Callum.Macdonald@utas.edu.au or CallumRMacdonald@hotmail.com
#If another author contributes new methods to the script, please add your name and contact here

"""
A script that will take some parameters and a list of spectra .bin files, and then produce processed spectra
Will either produce a simple on-off plot, where the on and off spectra are selected from the list
    or, will produce an spectra normalised by the average of the other spectra (not implemented yet)
    or, will just plot a single spectra. #TODO: FIX THIS BLURB
"""
#TODO: Implement median window filter, plot option 3. Will require a bit of work, but nothing too difficult
    #Needs to check if any of the 'off' scans have spectral features in the same frequency space as the on scan
    #Will need to calculate every velocity correction first, and then convert these to frequency shifts to check
    #Will look a bit different in implementation to on-off scans.
    #On-off will still be useful for intensive observations, like monitoring a single source during a flare.
#TODO: Add proper help print out
#TODO: Make a script that can read a schedule file and produce the inputs required from it (list of source names with on off stuff)
#TODO: Make option to print out the baseline fit. Have it name it properly.
#TODO: Need to add a warning message if the number of scans entered is not a multiple of the obs list
#TODO: Potentially implement internal intensity calibration. Given stability of 12m, might only need intensity ~twice a month
    #In that case, probably not worth it.
#TODO: Add error catcher to allow all plot files to be deleted from within the program, so that new ones can be produced
    #if the script is called in a directory with plot files in it.
#TODO: Finish error catchers for inputs. More informative than just quiting
    #Need: Epoch, velocity, spectra list, file format, intensity scale,
#TODO: Might be good to have custom inclusion sizes for baseline and rms fitting for each source
    #Would also be really good to implement a way for custome boundries, so that spectra with two distinct
    #spectral features can have fitting in between, otherwise have huge gap and baseline is poor.
#TODO: Implement a more modular form of data input. Currently works well for SWSpec, but it shouldn't be too hard to define
    #A general function that produces the outputs for different input formats. Only need intensity, frequency, and a few other bits
    #Same way as there is spt and csv output, have .bin, ect... input method
#TODO: Polyfit gives warnings for poor optimisation in the order of 10+. Not sure if any other package might be able to do better,
    #but if one exists, it would be good to replace it.
#TODO: Peakutils (package) might have some good things to add here. Might be able to do a better job with the baseline
    #Than my manual job
#TODO: Would be good to have averaging over IFs and duplicate scans in this script, so that baselines can be fit to the 
    #resulting spectrum, producing better fits.
#TODO: Implement a gain curve parameter. Corrects for gain change cause by telescope distortion, I didn't use it during
    #my project (oops)
#TODO: Add observation order to the epoch dictionary. This way, if there are several different orders used for future monitoring,
    #it will be easy to reprocess data without having to try and determine which order was used.
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
            default='csv'
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

        #Intensity scales
        #The scale factors to convert on-off quotient values to intensity
        #One factor for each IF, a,b,c,d,e,f
        #Default all set to 1 (can be set to a rough approximation as scale doesn't change by much over time)
        parser.add_argument(
            "-int_s",
            dest="intensity_scales",
            nargs="+",
            default=(1, 1, 1, 1, 1, 1,)
        )

        #Epoch of the observation
        #The date of the observation as a string in fits format'YYYY-MM-DDTHH:MM:SS'
        #Example '2021-03-31T01:30:00'
        #This will be used to remove the relative velocity of the earth with respect to LSR
        parser.add_argument(
            "-epoch",
            dest="epoch",
            default=None
        )

        #Name of observation
        #Name of the observation used to name the log file
        parser.add_argument(
            "-name",
            dest="name",
            type=str,
            default='' #Default set to no name
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
        parser.add_argument(
            "-vr",
            dest="velocity_range",
            nargs="+",
            default=[-20,20] #defualt set to standard maser zoom
        )

        #Velocity inclusion size
        #This is the width around the central velocity used for calculating the rms and baseline, in km/s
        #In these calculations, the exclusion window specified in the source dictionary will be excluded, 
        #   and everything else within the inclusion velocty will be included in the calculations
        parser.add_argument(
            "-inc",
            dest="velocity_inclusion_size",
            type=float,
            default=21 #defualt set to a reasonable window. Larger windows will require a higher order polynomial
        )
        
        #Baseline polynomial order
        #The order of the polynomial used in the baseline fitting process
        #The size of the inclusion winodw should be considered when choosing an order. If a larger window is 
        #   used, a higher order polynomial will be required for a good fit
        parser.add_argument(
            "-bpo",
            dest="baseline_polynomial_order",
            type=int,
            default=9 #defualt set to 7 for a reasonable fit over 40 km/s
        )

        #The following bool type arguments are read as strings. This is becasue bool measures the trueth value, so anythong other than '' is true
            #Will just use a string value for if statements in the code
        
        #Remove baseline
        #A bool to used to determine if the baseline should be removed
        #Baseline will be removed if true, ortherwise the spectra will be plotted as is
        parser.add_argument(
            "-rbl",
            dest="remove_baseline",
            type=str,
            default='True' #Baseline removed by default
        )

        #Normalise
        #A bool used to determine if the spectra is to be normalised to a peak intensity of 1
        #Will overwrite any intensity scaling and set the intensity of the strongest peak in the spectra to 1
        #Usefull for investigation of relative intensity changes between peaks over time.
        parser.add_argument(
            "-norm",
            dest="normalise",
            type=str,
            default='False' #Not normalised by default
        )

        #Standardise velocity
        #A bool to used to determine if the velocity should be standardised using fourier methods and specific sampling
        #If true, one of the velocity bins will be set to 0, and all others will be an integer multiple of the 
            #velocity resolution
        parser.add_argument(
            "-sdv",
            dest="standardise_velocity",
            type=str,
            default='True' #Velocity standardised by default
        )

        #Rounding option
        #A correction used for the early fm observations where the proc file was set up without full frequency specification
        #To get the correct velocity results, need to round the maser frequency to match what was entered in the proc file
        #TODO: This cannot be easily fixed, in future might be good to stardardise all to be rounded or rounded down.
        parser.add_argument(
            "-round",
            dest="frequency_rounding",
            type=str,
            default='True' #defualt set to True, might be fixed in future, but always on for now
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
            #Two defaults used for 7.2 and 7.6 GHz
            default=(None, None,(6181.128, 6668.5192, 7283.449, 7830.864),(6181.128, 6668.5192, 7283.449, 7830.864),(11964.007, 12178.595, 12229.348, 12329.666),(11964.007, 12178.595, 12229.348, 12329.666), None, None) 
            #default=(None, None,(6181.128, 6668.5192, 7682.232, 7830.864),(6181.128, 6668.5192, 7682.232, 7830.864),(11964.007, 12178.595, 12229.348, 12329.666),(11964.007, 12178.595, 12229.348, 12329.666), None, None) 
            #Set up for FM observations (A,B not used, four channels in each of cdef, cd and ef used with same freqs, as they are complamentary polaristions), no higher used (left in as example though)
            #For obs before fm023, need to exchnage 7.2 GHz with 7682.232 MHz transition)
        )

        #End of input arguments
        #Now extract variables out of the parser for simplicity
        #TODO: Have a print out of the entered parameters
        args = parser.parse_args() #args now is a list of all of the input variables
        self.bandwidth = int(args.bandwidth) * 1e6 #Convert bandwidth to Hz
        self.spec_list = args.spec_list
        self.format = str(args.format)
        self.name = str(args.name)
        self.observation_length = int(args.observation_length)
        self.observation_order = args.observation_order
        self.plot_option = int(args.plot_option)
        self.IF_chan_frequencies_MHz = args.IF_chan_freqs
        self.min_plot_vel = int(args.velocity_range[0]) #extract the minimum and maximum velocity to plot. Need int, as it reads as string
        self.max_plot_vel = int(args.velocity_range[1])
        self.velocity_inclusion_size = float(args.velocity_inclusion_size)
        self.baseline_polynomial_order = int(args.baseline_polynomial_order)
        self.remove_baseline_bool = str(args.remove_baseline)
        self.normalise = str(args.normalise)
        self.standardise_velocity = str(args.standardise_velocity)
        self.frequency_rounding = str(args.frequency_rounding)
        self.termination_message = 'Script terminating' #Message to be used when the script encounters a potential error and quits. Can be changed later to something more interesting.

        #If the spec_list is empty, there is no data to be plotted, and the program quits
        if len(self.spec_list) < 1:
            print('Must have at least one spectrum passed using -spec')
            print(self.termination_message)
            quit()

        #Define some dictionaries for use in some of the functions used to plot
        
        #A dictionary to hold the epochs of the fm observations
        #I added this to make it easier to process the data for analysis in my thesis, as otherwise I would need to look
            #up and copy the data each time I processed each epoch
        #Contains fm*** for all succesful observations, with the date and intensity as the dictionary entry
        self.fm_epochs = {
            'fm004' : ('2021-05-04T13:19:39', [1, 1, 5790, 5790, 1, 1]),
            'fm005' : ('2021-05-17T12:28:32', [1, 1, 7124, 7124, 1, 1]),
            'fm006' : ('2021-05-22T12:08:52', [1, 1, 5928, 5928, 1, 1]),
            #'fm007' : ('2021-05-26T11:53:09', [1, 1, 1, 1, 1, 1]),
            'fm008' : ('2021-05-30T11:37:25', [1, 1, 7437, 7437, 1, 1]),
            'fm009' : ('2021-06-04T11:17:45', [1, 1, 7029, 7029, 1, 1]),
            #'fm010' : ('2021-06-08T11:02:02', [1, 1, 1, 1, 1, 1]),
            #'fm011' : ('2021-06-09T10:58:06', [1, 1, 1, 1, 1, 1]),
            'fm012' : ('2021-06-13T10:42:22', [1, 1, 6640, 6640, 10681, 10681]),
            'fm014' : ('2021-06-27T09:47:20', [1, 1, 6262, 6262, 8224, 8224]),
            'fm015' : ('2021-06-30T09:35:32', [1, 1, 6698, 6698, 7960, 7960]),
            'fm016' : ('2021-07-04T09:19:48', [1, 1, 6382, 6382, 7766, 7766]),
            'fm017' : ('2021-07-07T09:08:00', [1, 1, 5251, 5251, 6923, 6923]),
            #Test'fm017' : ('2021-07-07T09:08:00', [1, 1, 1, 1, 6923, 6923]),
            'fm018' : ('2021-07-12T08:48:21', [1, 1, 5899, 5899, 7711, 7711]),
            'fm019' : ('2021-07-17T08:28:41', [1, 1, 6056, 6056, 1, 1]),
            'fm021' : ('2021-07-21T08:12:58', [1, 1, 5949, 5949, 1, 1]),
            #'fm022' : ('2021-07-26T07:53:18', [1, 1, 1, 1, 1, 1]),
            #ROUGH INTENSITIES BELOW, FIX THIS
            'fm023' : ('2021-08-02T07:25:47', [1, 1, 6100, 6100, 7500, 7500]),
            'fm024' : ('2021-08-09T07:28:11', [1, 1, 6100, 6100, 7500, 7500]),#Started at sidereal 14:30 rather than usual 14:00
            'fm029' : ('2021-10-02T03:25:56', [1, 1, 4300, 4300, 5847, 5847]),
            'fm030' : ('2021-10-05T03:14:09', [1, 1, 4700, 4700, 5343, 5343]),
            #Two following intensities calulcated through extrapolation
            'fm031' : ('2021-10-13T03:42:32', [1, 1, 4480, 4480, 5039, 5039]), #Started at sidereal 15:00
            'fm032' : ('2021-10-17T02:26:58', [1, 1, 5080, 5080, 5343, 5343]), 
            'fm035' : ('2021-10-24T01:59:26', [1, 1, 1, 1, 1, 1]),

            #G9 Intense monitoring
            #Intensities calulated by fixing intensity of secondary feature in G9
            'fm033' : ('2021-10-22T04:06:58', [1, 1, 1, 1, 6563, 6563]), #6 Broken
            #'fm034' : ('2021-10-23T05:40:47', [1, 1, 1, 1, 1, 1]), #Broken
            'fm036' : ('2021-10-25T03:57:10', [1, 1, 4064, 4064, 5962, 5962]),
            'fm037' : ('2021-10-26T03:51:15', [1, 1, 4260, 4260, 6387, 6387]),

            #fi intensity calibration epochs
            'fi001' : ('2021-10-13T01:12:56',  [1, 1, 1, 1, 1, 1]),
            'fi002' : ('2021-10-14T00:09:10', [1, 1, 1, 1, 1, 1]),
            'fi003' : ('2021-10-16T23:57:22', [1, 1, 1, 1, 1, 1]),
        }

        #If the epoch passed in is just a fm*** code, then use the dictionary entry. Otherwise, use it as a datetime string
        if self.fm_epochs.get(str(args.epoch)) != False:
            self.epoch = self.fm_epochs[str(args.epoch)][0] #Epoch is the first entry (0th psoition)
            self.intensity_scales = self.fm_epochs[str(args.epoch)][1] #Intensities is the second item
        else:
            self.epoch = str(args.epoch)
            self.intensity_scales = args.intensity_scales
        
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
        

        #For ease of re-ordering observations, each source name can be entered as a list, and all of the details 
            #will be extracted using the below dictionary
        #If another naming convention makes more sense, this can be changed easily
        #Each entry in the dictionary has the name as a label, and a tuple of the form (vel_centre, (vel_low_6GHz, vel_high_6GHz), (vel_low_12Ghz, vel_high_12GHz), RA, Dec, full_name)
        #The upper and lower velocities will be used to for baseline fitting and RMS calculations to exclude the spectral features
            #There are seperate entries for the different frequencies, as the resolution of the 12 GHz makes fitting
            #The wider baseline more difficult, and it is often noisier to begin with. Limits from 12 GHz mmb, narrowed for some sources
            #These might need to be changed for a higher sensitivity antenna, as it may exclude features we have no hope of 
            #seeing with the 12m telescopes
        #The RA is in the form **h**m**.*s, and the Dec is in the form ***d**m**.*s (hour, min, sec and deg, min, sec respectively)
        #All information from Methanol Multibeam survey, with some exlcusion windows narrowed as we have lower sensitivity (they included small features that are completely missed by us)
        #Any other information required for future features can be added easily.
        self.source_params = {
            'off' : (0, (-5, 5), (-5, 5), '0h0m0s', '0d0m0s', 'none'),
            'G328.237' : (-44.7, (-47.0, -31.5), (-46.2, -42.5),  '15h57m58.28s', '-53d59m22.7s', 'G328.237-0.547'), #6 GHz includes 328.254, which is offset from 237. Cannot resolve in 12 GHz, so not included
            'G328.809' : (-44.5, (-47.5, -42.0), (-47.5, -42), '15h55m48.70s', '-52d43m05.5s', 'G328.809+0.633'),
            'G351.417' : (-10.4, (-12.0, -4.5), (-12.1, -7.8), '17h20m53.37s', '-35d47m01.2s', 'G351.417+0.645'),
            'G323.740' : (-50.5, (-54.0, -47.0),(-53.0, -47.0), '15h31m45.45s', '-56d30m50.1s', 'G323.740-0.263'), #MMB has ver large exclusion size, might be a typo. I was getting fitting errors with the baseline
            'G318.948' : (-34.6, (-39.0, -31.0),(-37.0, -31.0), '15h00m55.40s', '-58d58m52.1s', 'G318.948-0.196'),
            'G9.621'   : (1.3, (-3.0, 6.0), (-2.2, 3.2), '18h06m14.67s', '-20d31m32.4s', 'G9.621+0.196'),
            'Virgo'    : (0, (-0.1, 0.1), (-0.1, 0.1), '12h30m49.42s', '12d23m28.04s', 'Virgo_A_Calibrator'), #Calibrator. Small exclusion window, not sure if the program would crash with 0 window
            #Expanded list, not used for the original fm observations   
        }


        #Check if the input .bin spectra files exist, and create a list of file names and a list of the paths to these spectra
        #The name will be used to determine which scan, IF and channel the spectra is of.
        #The following lines of code are largely based on Guifre's code
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.file_name_list = []
        self.path_list = []
        for i in range(len(self.spec_list)):
            if os.path.exists(self.spec_list[i]):
                path, file_name = os.path.split(self.spec_list[i]) #split the path and file name
                self.file_name_list.append(file_name)
                self.path_list.append("{}/".format(os.path.abspath(path)))
            else:
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), self.spec_list[i])
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    

    #The following function was NOT written by me (Callum Macdonald) and was sourced from github.
    #This function returns the barycentric velocity correct, which is used in the velocity_correction function below
    #Function from https://gist.github.com/StuartLittlefair/5aaf476c5d7b52d20aa9544cfaa936a1 on github
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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




        #The following code was used to experiment with different types of velocity correction.
        #After some testing, none of these tests were consistenly better than the bary_velocity from above
        """
        print('Old Bary ' + str(bary_velocity))
        #Convert RA and dec into deg
        d_dec, m_dec, s_dec, space_filler = re.split('[dms]', str(dec))
        dec_deg = (float(d_dec) + float(m_dec)/60 + float(s_dec)/(60*60))

        h_ra, m_ra, s_ra, space_filler = re.split('[hms]', str(RA))
        ra_deg = (float(h_ra) + float(m_ra)/60 + float(s_ra)/(60*60)) * 360 / 24

        helio_velocity = PyAstronomy.pyasl.helcorr(obs_long = long, obs_lat = lat, obs_alt = alt, ra2000 = ra_deg, dec2000 = dec_deg, jd = epoch_julian)
        bary_hel_velocity = PyAstronomy.pyasl.baryCorr(ra = ra_deg, dec = dec_deg, jd = epoch_julian)

        print('Helio ' + str(helio_velocity))
        print('Bary corr ' + str(bary_hel_velocity))
        #print(epoch_julian)

        bary_correction = bary_hel_velocity[1] - bary_hel_velocity[0]
        corrected_bary_velocity = bary_velocity + bary_correction
        print('Bary corrected ' + str(corrected_bary_velocity))
        bary_velocity = -float(bary_hel_velocity[1])
        """

        #To convert to LSR, define the direction of the LSR, then take the dot product of the unit vectors
        #This result can be multiplied by 20 km/s to get LSR velocity conversion factor
        #Given on ATNF sky freq calculator, 20 km/s in direction RA 18 hour (270 deg), Dec + 30 deg, needs b1900 equinox, so use fk4
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
        #print('LSRK correction ' + str(lsr_correction))

        #Combine the two velocities to get the LSR corrected velocity
        LSR_velocity = bary_velocity + lsr_correction
        #print('Final ' + str(LSR_velocity))

        return LSR_velocity
    

    #Method to generate the RMS noise level for data
    #Will take an intensity array, velocity array, centre velocity, exclusion size and a inclusion size
    #The exlusion min and max correspond to the velocities between which the data should be excluded (where the feature is)
    #The inclusion size will be centred on the centre velocity and everything within it other than the excluded region will
        #be used to determine the RMS noise level
    #Values should be given in km/s
    #The intensity and velocity arrays should be the same length
    #The exclusion size should be smaller than the inclusion size
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


    #Some of the code in the following section is based on the code written by Dr Patrick Yates for KYA320 
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
    def write_spt(self, file_name, min_vel, max_vel, data_point_count, min_flux, max_flux, velocity_array, intensity_array, x_label = 'Velocity w.r.t. LSRK (kms\\u-1\\d)',y_label = 'Flux Density (Jy)'):
        
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
    def write_csv(self, file_name, min_vel, max_vel, velocity_array, intensity_array, epoch, x_name = 'vel(km/s)', y_name = 'intensity(Jy)'):
        
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
    #TODO: This method is meaty, might be good to break it up into some smaller methods. Could then have PO 1 2 3 not repeated
    def plot_spectrum(self, velocity_array, spectra_array_on, spectra_array_off, fourier_counts, fourier_velocity, intensity_scale, plot_name, vel_source, vel_source_low, vel_source_high, epoch, format, inclusion_size, polynomial_order):
        #Write out the file to a .spt file, to be used with Simon's script 'splot.pl'

        #Different plotting options are done seperately           

        #TODO: FOURIER IN OPTION 1, also <30 if statement clean
        #TODO: Option 1 and two can probably be run from the same lines of code, with a small if statement for the uotient section
        if self.plot_option == 1:
            print('Plot option 1 selected for ' + plot_name)

            #Need to dervie values that will be used to define the plot parameters
            #Also need to produce the intensity array for on/off, and then remove the baseline

            #Specifiy the plotting window base on the velocity of the spectral feature
            vel_plot_lower = vel_source + self.min_plot_vel
            vel_plot_upper = vel_source + self.max_plot_vel

            intensity_array = spectra_array_on * intensity_scale

            #If the velocity is to be standardised, do that here
            if self.standardise_velocity == 'True':
                #Calculate the discrete cosine transform coefficients using scipy
                #Because of how velocity is calculated, need to flip the intensity array
                intensity_array = np.flip(intensity_array)

                #Now calulcated the dct coefficients
                #These can then be used to resample the spectra at the specified points
                #Type 3 is used because it was the only one I could get to work correctly.
                    #The documentation for the fft package was a bit disapointing, if anyone has a better implementation
                    #to offer please let me know
                dct_intensity = scipy.fft.dct(intensity_array, type=3)
                
                #Caclulating this for all 32 or 64 thousand values takes a long time.
                #TODO: Test if shortening the intensity array still produces the correct outputs.
                #Instead, shorten the velocity and count arrays.
                #Shorten them to slightly longer than they will be by adding a fraction of the total plotting window distance
                safety_factor = 0.4 * (vel_plot_upper - vel_plot_lower)

                #Create two empty arrays to append values to
                shortened_velocity = []
                shortened_counts = []
                for i in range(len(fourier_velocity)):
                    if (fourier_velocity[i] < vel_plot_upper + safety_factor) and (fourier_velocity[i] > vel_plot_lower - safety_factor):
                        #If the velocity is within the (extended with safety factor) range, then append the values.
                        shortened_velocity.append(fourier_velocity[i])
                        shortened_counts.append(fourier_counts[i])

                #Convert the shortened arrays into numpy arrays
                shortened_velocity = np.array(shortened_velocity)
                shortened_counts = np.array(shortened_counts)

                #Now add the intensity calculated from each cosine term.
                #First create an array of zeroes to add the values to.
                intensity_from_dct = np.zeros(len(shortened_counts))

                for i in range(len(dct_intensity)):
                    #The first bracketted term is the normalised amplitude of the ith cosine term
                    #The term inside the cosine is taken from the formula provided in scipy. (This formula is the one for type 2, but did not work with type 2 and did with 3. Again, a bit janky, but it works)
                    intensity_from_dct += (dct_intensity[i] / len(dct_intensity))*np.cos(shortened_counts * np.pi* (2*i +1) / (2*len(fourier_counts)))
                #Now procede as normal wiht the function, but with the intensity and velocity replaced by the dct ones
                intensity_array = intensity_from_dct
                velocity_array = shortened_velocity

            #Calculate the polynomial coefficients for the baseline 
            baseline_coefficients = self.fit_baseline(intensity_array = intensity_array, velocity_array=velocity_array, exclusion_min = vel_source_low, exclusion_max = vel_source_high, centre_velocity = vel_source, inclusion_size = inclusion_size, polynomial_order = polynomial_order)

            #Use the coefficients to remove the baseline using remove_baseline(). Skip this step if remove_baseline bool is False
            #rb for removed baseline
            if self.remove_baseline_bool == 'True':
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

            #Calculate the RMS noise of the spectra
            RMS = self.calculate_rms(intensity_array = intensity_array_rb, velocity_array = velocity_array, centre_velocity = vel_source, exclusion_max = vel_source_high, exclusion_min = vel_source_low, inclusion_size = inclusion_size)

            #If the spectra is to be normalised, do so.
            #This is done after RMS calculation so as to give RMS in Jy is scaling is provided
            if self.normalise == 'True':
                intensity_array_rb = intensity_array_rb / largest_flux_val
                flux_min = -0.2 #Scale the spectra window to fit the spectra properly.
                flux_max = 1.2

            #Write the file out to the chosen data format
            if format == 'spt':
                self.write_spt(file_name = plot_name, min_vel = vel_plot_lower, max_vel = vel_plot_upper, data_point_count = data_point_count, min_flux = flux_min, max_flux = flux_max, velocity_array = velocity_array, intensity_array = intensity_array_rb)
            elif format == 'csv':
                self.write_csv(file_name = plot_name, min_vel = vel_plot_lower, max_vel = vel_plot_upper, velocity_array = velocity_array, intensity_array = intensity_array_rb, epoch = epoch)
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
                #divide the ith value from the on data by the ith value from the off data, then subtract 1 (quotient formula)
                quotient_array[i] = spectra_array_on[i] / spectra_array_off[i] - 1 

            #Need to convert this ratio of on and off into a flux density.
            #This is done by calibrating the telescope
            intensity_array = quotient_array * intensity_scale
            
            #If the velocity is to be standardised, do that here
            if self.standardise_velocity == 'True':
                #Calculate the discrete cosine transform coefficients using scipy
                #Because of how velocity is calculated, need to flip the intensity array
                intensity_array = np.flip(intensity_array)

                #Now calulcated the dct coefficients
                #These can then be used to resample the spectra at the specified points
                #Type 3 is used because it was the only one I could get to work correctly.
                    #The documentation for the fft package was a bit disapointing, if anyone has a better implementation
                    #to offer please let me know
                dct_intensity = scipy.fft.dct(intensity_array, type=3)
                
                #Caclulating this for all 32 or 64 thousand values takes a long time.
                #TODO: Test if shortening the intensity array still produces the correct outputs.
                #Instead, shorten the velocity and count arrays.
                #Shorten them to slightly longer than they will be by adding a fraction of the total plotting window distance
                safety_factor = 0.2 * (vel_plot_upper - vel_plot_lower)

                #Create two empty arrays to append values to
                shortened_velocity = []
                shortened_counts = []
                for i in range(len(fourier_velocity)):
                    if (fourier_velocity[i] < vel_plot_upper + safety_factor) and (fourier_velocity[i] > vel_plot_lower - safety_factor):
                        #If the velocity is within the (extended with safety factor) range, then append the values.
                        shortened_velocity.append(fourier_velocity[i])
                        shortened_counts.append(fourier_counts[i])

                #Convert the shortened arrays into numpy arrays
                shortened_velocity = np.array(shortened_velocity)
                shortened_counts = np.array(shortened_counts)

                #Now add the intensity calculated from each cosine term.
                #First create an array of zeroes to add the values to.
                intensity_from_dct = np.zeros(len(shortened_counts))
                for i in range(len(dct_intensity)):
                    #The first bracketted term is the normalised amplitude of the ith cosine term
                    #The term inside the cosine is taken from the formula provided in scipy. (This formula is the one for type 2, but did not work with type 2 and did with 3. Again, a bit janky, but it works)
                    intensity_from_dct += (dct_intensity[i] / len(dct_intensity))*np.cos(shortened_counts * np.pi* (2*i +1) / (2*len(fourier_counts)))
                
                #Now procede as normal wiht the function, but with the intensity and velocity replaced by the dct ones
                intensity_array = intensity_from_dct
                velocity_array = shortened_velocity
            

            #Calculate the polynomial coefficients and remove the baseline using remove_baseline(), if it is to be removed
            #rb for removed baseline
            if self.remove_baseline_bool == 'True':
                baseline_coefficients = self.fit_baseline(intensity_array = intensity_array, velocity_array=velocity_array, exclusion_min = vel_source_low, exclusion_max = vel_source_high, centre_velocity = vel_source, inclusion_size = inclusion_size, polynomial_order = polynomial_order)
                intensity_array_rb = self.remove_baseline(intensity_array=intensity_array, velocity_array=velocity_array, polynomial_coefficients=baseline_coefficients)
            else:
                baseline_coefficients = 'Not removed'
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
            #This will likely need to be changed if Ceduna is used, as the noise will be much less than ~20 Jy
            if largest_flux_val > 30:
                flux_min = -0.2 * largest_flux_val
                flux_max = 1.2 * largest_flux_val
            else:
                flux_min = -30
                flux_max = 30

            #Caclulate the RMS noise for the spectra.
            RMS = self.calculate_rms(intensity_array = intensity_array_rb, velocity_array = velocity_array, centre_velocity = vel_source, exclusion_max = vel_source_high, exclusion_min = vel_source_low, inclusion_size = inclusion_size)
            
            #If the spectra is to be normalised, do so.
            #This is done after RMS calculation so as to give RMS in Jy if scaling is provided
            #TODO: Sometime this produces a peak of 1 like expected, sometimes it's 0.997~. Float errors should be smaller than that, not sure why this is occuring
            if self.normalise == 'True':
                intensity_array_rb = intensity_array_rb / largest_flux_val
                flux_min = -0.2 #Scale the spectra window to fit the spectra properly.
                flux_max = 1.2


            #Write the file out to the chosen data format
            if format == 'spt':
                self.write_spt(file_name = plot_name, min_vel = vel_plot_lower, max_vel = vel_plot_upper, data_point_count = data_point_count, min_flux = flux_min, max_flux = flux_max, velocity_array = velocity_array, intensity_array = intensity_array_rb)
            elif format == 'csv':
                self.write_csv(file_name = plot_name, min_vel = vel_plot_lower, max_vel = vel_plot_upper, velocity_array = velocity_array, intensity_array = intensity_array_rb, epoch = epoch)
            
            #Return the RMS and baseline coefficients so that they can be recorded in the output log
            return(RMS, baseline_coefficients)


    def plot_all(self):
        #Two different procedures. One for on / on-off plotting, and another for plotting using every other spectra as baseline (NOT DONE YET)
        
        #Both procedures will have an associated log file, which will store information on the rms, baseline and general parameters used
        log_file = open("log_" + str(self.name) + "_" + str(self.format) + ".txt", "x")
        log_file.write("Epoch: " + str(self.epoch) + "\n")
        log_file.write("Scan duration (including source change time): " + str(self.observation_length) + " seconds\n")
        log_file.write("Intensity scaling factors (for ifs a,b,c,d,e,f): " + str(self.intensity_scales))
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
                    #For trying to look at off spectra, can just set obs_order to ['test'] and pass in single spectra in question
                if scan_tag != 'off':
                    #Now that it is know that the scan needs to be plotted, the other information can be extracted from the name
                    obs_name = split_name[0]
                    IF = split_name[2] #Skipped indice here is for filler_text that is in positions 1,3,5
                    channel_number = int(split_name[4])
                    observatory = split_name[6]
                    num_channels = int(re.split('[p]', split_name[8])[0]) #split at the p in pt, and then just take the first part which is the number 
                    integration_time = int(re.split('[s]', split_name[9])[0]) #Same as above for s after time 
                    processing_format = split_name[10]

                    #TODO: Add error catcher if off scan picked is going to be insuitable (dif length)
                    #TODO: Add an efficiency optimisation to check if the if, channel and channel numbers are the same as the last scan
                    #   In this case, velocity array will be the same and it doesn't need to be recalculated


                    #Derive parameters that will be used in the creation of the velocity array and later in the spectra data extraction
                    #Can incorperate factor of two into both numchannel lines, as it is just due to twice as many channels as there are data points
                    #TODO: Currently bandwidth cannot be changed for 12m, but if in future it can, this can be updated so that the bandwidth can vary between IFs, rather than being fix like it is here
                    
                    #The following five lines of code were based on Guifre's code
                    #They calculate parameters useful for producing frequency arrays corresponding to the output intensity
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

                    #Because of an issue with the early proc files (at least fm025 and earlier), need to adjust the frequency taken from ff
                    if self.frequency_rounding == 'True':
                        #This would be a simple case of rounding, but Gabor decied to round down for the c and d IFs
                        if (IF == 'c') or (IF == 'd'):
                            centre_freq = math.floor(maser_freq / 1e6) * 1e6 #Divide by 1e6 so that the floor function can be used to round down, then bring back to Hz
                        else:
                            centre_freq = round(maser_freq / 1e6) * 1e6
                    else:
                        centre_freq = maser_freq
                    
                    #Repeat this process of checking and extraction for the intensity scale
                    if(IF_number >= len(self.intensity_scales)):
                        print('Error trying to determine the intensity scale for file: ' +str(file_name))
                        print('The IF number: ' + str(IF_number) + ' is out of bounds for the inensity scaling factor tuple: ' + str(self.intensity_scales))
                        print(self.termination_message)
                        quit()
                    
                    intensity_scale = float(self.intensity_scales[IF_number])
                
                    #Convert full frequency to velocity
                    #As ff is centred around 16 MHz, to move that to be centred on centre freq, need to subtract 16 and add centre freq
                    
                    #Fudge factor used to get velocity to line up if need be
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
                    #If it is not contained in the list, then the spectrum cannot be generated and the program will quit
                    if self.source_params.get(scan_tag) == None:
                        print('Name entered in -names is not contained in the dictionary of sources for tag: ' + scan_tag + 'from file: ' + file_name)
                        print(self.termination_message)
                        quit()
                    
                    source_parameters = self.source_params[scan_tag]

                    #Extract the information from the source parameters                   
                    #The format is (vel_centre, (vel_low_6GHz, vel_high_6GHz), (vel_low_12GHz, vel_high_12GHz), RA, Dec, full_name)
                    vel_centre = source_parameters[0]
                    RA = source_parameters[3]
                    dec = source_parameters[4]
                    full_name = source_parameters[5]

                    #There are two exclusion windows for each source, corresponding to the different frequency bands
                    #IF c and d are 6 GHz, e and f are 12
                    if (IF == 'c') or (IF == 'd'):
                        vel_low = source_parameters[1][0]
                        vel_high = source_parameters[1][1]
                    elif (IF == 'e') or (IF == 'f'):
                        vel_low = source_parameters[2][0]
                        vel_high = source_parameters[2][1]
                    else:
                        print('Velocity windows for IFs other than c,d,e and f are not yet implemented.')
                        print(self.termination_message)
                        quit()
                    


                    #Calculate the time of the scan
                    #To do this, multiply the scan number -1 by the length of a single scan, and add that to the starting epoch
                    #(Scan number -1 because the the first scan starts at the starting time, and the second one scan later ect)
                    #Also add on a centre time, which can be used to shift the time at which the velocity is calculated during each observation. Currently 0
                    centre_time = int(self.observation_length / 2) #Use half the observation length to place the time at the middle of the observation.
                    epoch_numpy = np.datetime64(self.epoch) + np.timedelta64(scan_number_int * self.observation_length + centre_time, 's')
                    

                    #Calculate the appropriate velocity correction and apply it to the velocity array
                    velocity_correction = self.velocity_correction(lat = lat, long = lon, alt = alt, epoch = epoch_numpy, RA = RA, dec=dec)
                    #print(velocity_correction)
                    #Apply the velocity correction
                    corrected_velocity_array = vel_freq_array - velocity_correction

                    #The following lines of code are based upon code written by Guifre, without much change
                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    #Extract the spectra information from the .bin file
                    file_temp_on = open(f"{file_path}{file_name}", "rb") #Open the .bin file
                    spectra_array = np.fromfile(file=file_temp_on, dtype="float32", count=Nfft) #Extract the data from the file
                    file_temp_on.close #Close the file
                    
                    #If an on-off plot is to be produced (plot option 2), then extract the scan that comes next in the list of files
                        #The next file shoudl always be the corresponding off scan.
                        #TODO: Add an error catcher if to ensure that the scan number is the next in line.
                    if self.plot_option == 2:
                        file_name_off = self.file_name_list[i+1]
                        file_path_off = self.path_list[i+1]
                        file_temp_off = open(f"{file_path_off}{file_name_off}", "rb") #Open the .bin file
                        off_spectra = np.fromfile(file=file_temp_off, dtype="float32", count=Nfft) #Extract the data from the file
                        file_temp_off.close #Close the file
                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    elif self.plot_option == 1:
                        #If plot option 1 is selected, then there is no need for an off spectra to be extracted
                        #To facilitate this, set off_spectra to None, so that this body of code doesn't need to be repeated or awkwardly defined
                        #Otherwise this whole section of code would need to be written twice, one with off and one without
                        off_spectra = None
                    #To be in this body of code requires plot option to be 1 or 2 (all in an if statement), so this is everything
                    
                    #TODO: Add a checker that ensures that the selected off scan is an off scan, and maybe even ensure that it
                        #should be associated with this scan. Can just raise a warning if not.


                    #Use the information above the construct the name for the plot
                    plot_name = str(full_name) + "_" + str(obs_name) + "_IF:" + str(IF) + "_" + str(maser_freq / 1e6) + "MHz_scan:" + str(scan_number) + "_" + str(observatory)

                    #If the plot is on-off, add that to the end of the name
                    if self.plot_option == 2:
                        plot_name = plot_name + str("_on_off")

                    
                    #TODO: Move all of the normal velocity stuff into the else case to save computation time. Not urgent, only small amount of code
                    if self.standardise_velocity == 'True':
                        #The following documentation is a little lengthy, but I found this section a bit confusing to construct, so I wanted
                            #To ensure that I could remember why I made these decissions later.
                        #The basic idea is to use a fourier transform (or in this case a descrete cosine transform) to construct the
                            #sum of functions required to construct the signal.
                            #We then use a standardised velocity array to sample from these. This array is contructed from an
                            #integer array that is shifted by a decimal value to compensate for the shifts used in further processing
                            #Convert from the integer array by multiplying by the velocity channel width and adding all of the terms
                            #that change the raw frequency input into the velocity-around-transition format
                        #The constant value is added so that once this array is converted to velocity units and centred it will
                            #be centred on the value 0, ie -2dv, -dv, 0, dv, 2dv rather than -2dv +delta, -dv +delta, 0 +delta ect
                        #This facilitates the easy comparison and averaging of different spectra, as a velocity range will give the same
                            #number of bins given the same dv
                        #Cannot just add rounded vel-cor, as that would sample at the same place, then not shift it to the correct amount
                        #The ammount that these bins need to be shifted by is the remainder of the velocity correction (LSRK) w.r.t delta v
                            #and the shift caused by the centre and transition frequency not being the same.
                        #This way the bins will be at the same location every observation, with a bin at velocity = 0.0

                        #The dv unit is given by the following equations
                        #df/freq = dv/c, BW/df = num_chan / 2, => dv = 2*c*BW/(num_chan*freq), for c=2.998*10^5 km/s
                        delta_vel = 2 * 2.998*1e5 * self.bandwidth / (num_channels * maser_freq)
                        
                        #Calculate the remainders that the bins need to be shifted by to have a 0 velocity value
                        #These are the non-integer components of the velocity and frequency corrections.
                        #This is different to just rounding the corrections, as the fourier transform is sampled by integer_counts, with no reference to velocity until later.
                        correction_remainder = (velocity_correction % delta_vel) / delta_vel#use the % python function to calculate the remainder                       
                        centre_maser_offset_bins = (centre_freq - maser_freq) / df #df is thw bin width in frequency space, calculated as part of Guifre's code
                        centre_maser_remainder = centre_maser_offset_bins%1
                        integer_counts = np.arange(0, Nfft, 1) + correction_remainder + centre_maser_remainder
                        
                        #Now create the velocity array that corresponds to the points sampled by this array
                        #Want the zero value to be centred, so first subtract half the array length
                        #Because of the difference in the centre frequency and the maser frequency, need to 
                            #subtract a small factor based on the number of bins (diff in freq / bin width)
                        
                        shifted_velocity = (integer_counts - int(len(integer_counts)/2) - centre_maser_offset_bins) * delta_vel - velocity_correction
                        

                        #All velocity arrays at this stage should be the same, but they are different due to small 
                            #errors introduced by computer float maths. Round to make consistent
                        #The level rounded to will be far smaller than the uncertainty, so it shouldn't be an issue.
                        shifted_velocity = np.around(shifted_velocity, 9)

                    else:
                        #If the velocity standardisation is not required, create empty variables to be passed in the write_plot argument
                        integer_counts = 'Velocity correction not used'
                        shifted_velocity = 'Velocity correction not used'
                    

                    #Can now use write_plot to turn this information into a spectra
                    #Write plot will return the rms and polynomial baseline coefficients
                    RMS, baseline_coefficients = self.plot_spectrum(velocity_array = corrected_velocity_array, spectra_array_on = spectra_array, spectra_array_off = off_spectra, fourier_counts=integer_counts, fourier_velocity=shifted_velocity, intensity_scale = intensity_scale, plot_name = plot_name, vel_source = vel_centre, vel_source_low = vel_low, vel_source_high = vel_high, epoch = epoch_numpy, format = self.format, inclusion_size = self.velocity_inclusion_size, polynomial_order = self.baseline_polynomial_order)
                    #RMS, baseline_coefficients = (0,0)
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