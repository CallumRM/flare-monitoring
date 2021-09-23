import numpy as np
import matplotlib.pyplot as plt
import os, argparse, re
import errno
import re
import sys
import PyAstronomy.pyasl
from astropy import units
import astropy.time
import astropy.coordinates

"""
A script that will take some parameters and a list of spectra .bin files, and then producea plot
Will either produce a simple on-off plot, where the on and off spectra are selected from the list
    or, will produce an spectra normalised by the average of the other spectra
    or, will just plot a single spectra.

"""
#TODO: Perl output, Do on-off or baseline for multiple inputs. Combine scans, should be easy, just need to convert vel before merging
#Will need to take in the start time of the observing epcoh in numpy YYYY-MM-DD'T'HH:MM:SS form, and then have scan length, use numpy datetime, and then use timedelta to add extra
#This way all data can be velocity corrected before before being plotted.
#TODO: Need to set intensity, will need to fit guassian and then use that to find conversion factor
    #Or use oother intensity calibration that Simmon talked about
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

        #Bandwidth argument
        #Bandwidth as an int in MHz
        parser.add_argument(
            "-bw",
            dest="bandwidth",
            help="observing bandwidth in MHz", 
            type=float, 
            default=32 #default set to 32 for the AuScope 12m antennas
        )

        #Number of channels argument
        #For formatting, see bandwidth above. Use -nc 'number of channels'
        #Must be an integer (there are an integer number of channels)
        parser.add_argument(
            "-nc",
            dest="num_channels",
            help="The number of channels recorded",
            type=int,
            default=32000 #default set to 32000 which is the number being used for the current flare monitoring setup
        )

        #List of spectra argument
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

        #Central frequency of the scan in MHz
        #The scans spectra by default have a 0 to 32 MHz scale.
        #To convert to velocity, the frequency corresponding to the 16 MHz position must be known
        parser.add_argument(
            "-cf",
            dest="centre_freq",
            default=6668 #Default set to centre freq used in fm tests for 6.7 GHz
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
        #Default set for 10 minutes or 600 seconds
        parser.add_argument(
            "-length",
            dest="observation_length",
            default=600
        )


        #Location of observatory
        #Either to initials of the included observaatories, or the lat, lon and altitude
        #Give coordinates in degrees, altitude in metres
        #Default set to Hobart 12 m
        #TODO: Allow DMS form to be given
        #Included observatories are:
        #   hb: Hobart 12 m
        #TODO: Add other utas facilities (ke, yg, ceduna ect)
        parser.add_argument(
            "-obs",
            dest="observatory",
            nargs="+",
            default='parkes'
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
            default=[-200,200] #defualt set to catch most all sources that will be observed. Will be very zoomed out
        )

        #Observation names
        #The names of the observations in a list (list of strings)
        #Default set for fm observations of 6 sources
        #If no names are provided, the default (which probably wont match) will be used
        #Name the off files to any string, but they need to be there as spacers
        #TODO IF there is a better way to do this, implement it
        #TODO Add something like group name, so that fm*** can be added to the end of each name.
        parser.add_argument(
            "-names",
            dest="observation_names",
            nargs="+",
            default=['G328.237-0.547', 'off', 'G328.809+0.633', 'off', 'G351.417+0.645', 'off', 'G323.740-0.263', 'off', 'G318.948-0.196', 'off', 'G9.621+0.196', 'off'] 
        )

        #Sky coordinates of the sources
        #Needs to be entered in degrees (RA, Dec)
        #Coordinates sourced from multibeam, converted using https://www.swift.psu.edu/secure/toop/convert.htm
        #TODO:This doesn't have split second accuracy, need to create internal converter and just enter HMS, D'" form
        parser.add_argument(
            "-coords",
            dest="coordinates",
            nargs="+",
            default=[[239.4917, -53.9897], 'off', [238.9542, -52.7183], 'off', [260.2208, -35.7836], 'off', [232.9375, -56.5139], 'off', [225.2292, -58.9811], 'off', [271.5625, -20.52567], 'off'] 
        )

        #Source velocities
        #A list of the recession velocities of the sources observed in Km/s
        #Default set for fm observations of 6 sources, velocities sourced from multibeam catalogues
        #For the G351.417 pair, use the velocity of +0.645, very close to 0.646
        #TODO: Finish this please
        parser.add_argument(
            "-velocities",
            dest="source_velocities",
            nargs="+",
            default=[-44.7, 'off', -44.5, 'off', -10.4, 'off', -50.5, 'off', -34.6, 'off', 1.3]
            
        )
        
        #Show spectrum
        #Boolean value. If true, the spectrum will be displayed.
        #If false, the spectrum will be saved as a .svg
        #Default false for ease of use with automation
        #TODO: Separate this and save plot, so that both can be done simultaneously
        parser.add_argument(
            "-show_spec",
            dest="show_spec",
            default=False
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
            default=(1,2) #Default to 1 and 1, as there will always be at least one spectrum
        )

        #End of input arguments
        #Now extract variables out of the parser for simplicity
        #TODO: Have a print out of the entered parameters
        #TODO: Make a converter from RA and Dec in HMS to Deg
        args =parser.parse_args() #args now is a list of all of the input variables
        self.bandwidth = int(args.bandwidth) * 1e6 #Convert bandwidth to Hz
        self.num_channels = int(args.num_channels)
        self.spec_list = args.spec_list
        self.centre_freq = int(args.centre_freq) * 1e6 #Convert centre freq to Hz
        self.epoch = str(args.epoch)
        self.observation_length = int(args.observation_length)
        self.plot_option = int(args.plot_option)
        self.observation_names = args.observation_names
        self.observation_coordinates = args.coordinates
        self.source_velocities = args.source_velocities
        self.show_spec = bool(args.show_spec)
        self.min_plot_vel = int(args.velocity_range[0]) #extract the minimum and maximum velocity to plot. Need int, as it reads as string
        self.max_plot_vel = int(args.velocity_range[1])
        self.single_spec_position = int(args.single_spec_position) - 1 #Subract 1 as python lists start at 0
        self.on_position = int(args.on_off_positions[0]) - 1 #Get the on (first) value from the tuple, subract 1 as ptyhon lists start at 0
        self.off_position = int(args.on_off_positions[1]) - 1 #Get the off (second) value from the tuple, subract 1 as ptyhon lists start at 0

        #Extract the observatory location data
        #If a name has been provided, use the associated coords
        #Location data from (J. Lovell et.al 2012, The AuScope geodetic VLBI array)
        #TODO: FInish these for other observatories
        #Converted into decimal degrees
        if str(args.observatory) == 'hb':
            self.lat = -42.80557
            self.long = 147.4381
            self.alt = 40.967
        elif args.observatory == 'parkes': #Parkes 64m (For vel calibration tests using ATNF online tool, sourced from Parkes user guide)
            self.lat = -32.99839
            self.long = 148.26352
            self.alt = 415
        #If it is not a name, and instead the three pieces of coord data have been provided, use them
        elif len(args.observatory) == 3:
            self.lat = args.observatory[0]
            self.long = args.observatory[1]
            self.alt = args.observatory[2]
        #If the input is not specified correctly, quit the program
        #(If velocity doesn't matter, and hence no observatory is required, just use hb)
        else:
            print('Observatory data entered incorrectly using -obs. \n Either enter lat lon alt or the abreviated observatory name.\n Script terminating.')
            quit()


        #TODO: Add error catchers for all of the other values (ie plot option, spec_pos out of bounds)
        #If the spec_list is empty, there is no data to be plotted, and the program quits
        if len(self.spec_list) < 1:
            print('Must have at least one spectrum passed using -spec')
            print('Script terminating')
            quit



        #First, check if the input .bin spectra files exist, and create a list of file names and a list of the paths to these spectra
        #The splitting of name and path may not be necessary, but I don't know enough about opening files to make a call.
        self.file_name_list = []
        self.path_list = []
        for i in range(len(self.spec_list)):
            if os.path.exists(self.spec_list[i]):
                path, file_name = os.path.split(self.spec_list[i]) #split the path and file name
                self.file_name_list.append(file_name)
                self.path_list.append("{}/".format(os.path.abspath(path)))
            else:
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), self.spec_list[i])


        #Now calculate some variables and lists to be used for the plot
        #TODO: Work out exactly what these variables are, so as to better document the code. Currently this is ripped
            #from Guifre's code
        self.SR = int(2 * self.bandwidth)
        self.df = float(self.SR / self.num_channels) #Differential freq unit
        self.Nfft = int((self.num_channels / 2) + 1)
        self.jf = np.arange(self.Nfft, dtype=int)
        self.ff = self.df * self.jf

        #Convert full frequency to velocity
        #As ff is centred around 16 MHz, to move that to be centred on centre freq, need to subtract 16 and add centre freq
        #TODO: Edit this so that it works for any ff centre
        self.ff_cen = self.ff + (self.centre_freq - 16*1e6)

        #Now convert the frequency to velocity using the frequency of the transition
        #This velocty does not factor in the Earth's motion realtive to LSR
        #That will be done for each plot individually, as the correction will change over the course of the observation
        #TODO: Update this to workk with other freq transitions
        self.freq_maser = 6668.5*1e6
        self.general_vel = 2.998*1e8 * (1 - self.ff_cen / self.freq_maser) / 1000 #Use doppler eq, divide by 1000 for Km/s


        #Make a list of the spectra data from the .bin files
        self.spec_data_list = []
        for i in range(len(self.spec_list)): #For each file name passed in
            file_temp = open(f"{self.path_list[i]}{self.file_name_list[i]}", "rb") #Open the .bin file
            temp_data = np.fromfile(file=file_temp, dtype="float32", count=self.Nfft) #Extract the data from the file
            self.spec_data_list.append(temp_data) #Put array of data in the list
            file_temp.close #Close the file
    
    #Method to calculate the velocity correction given an observatory location, epoch and Ra and Dec
    #lat long in degrees
    #epoch numpy datetime64 object
    #RA string in **h**m**.*s format
    #Dec string in ***d**m**.*s
    #Will return a float of the velocity correction
    #TODO: Set up a HMS to Degree converter, so that everything can be entered in the standard coords
    #TODO: Work out why there is a slight disagreement between this and ATNF

    #Function from https://gist.github.com/StuartLittlefair/5aaf476c5d7b52d20aa9544cfaa936a1 on github
    def velcorr(self, time, skycoord, location=None):
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


    def velocity_correction(self, lat, long, alt, epoch, RA, Dec):
        #First step is to generate the Heliocentric velocity using PyAstronomy
        #Convert numpy datetime64 to astropy, so it can be converted to Julian time
        epoch_astro = astropy.time.Time(val = epoch)
        #print(epoch_astro)
        epoch_julian = epoch_astro.jd
        
        #Use PyAstronomy to generate the heliocentric velocity
        #The velocity given appears to be the negative of the los velocity
        #TODO: Check that this is correct
        helio_velocity = -PyAstronomy.pyasl.helcorr(obs_long = long, obs_lat = lat, obs_alt = alt, ra2000 =RA, dec2000 = Dec, jd = epoch_julian)[0]
        #print(helio_velocity)

        #Next step is to convert from heliocentric to LSR
        #This can be done with a simple calculation, but requires galatic coordinates
        #To covnert RA and Dec to l and b, use AstroPy
        sky_coords = astropy.coordinates.SkyCoord(RA*units.deg, Dec*units.deg, frame = 'icrs')
        sky_coords = astropy.coordinates.SkyCoord('18h06m14.67s', '-20d31m32.4s', frame = 'icrs')


        #Test out new function from github person
        #First need to use astropy location object,
        obs_location = astropy.coordinates.EarthLocation.from_geodetic(long, lat, height=alt)
        git_vel = self.velcorr(time=epoch_astro, skycoord = sky_coords, location = obs_location)
        #This velocity comes with km / s after, need to split it out
        bar_vel, space_filler , space_filler , space_filler = re.split('[ ]', str(git_vel))
        bary_velocity = -float(bar_vel)
        print(bary_velocity)

        #Extract l and b (in radians for use in calculation below)
        galactic_l_dms = sky_coords.galactic.l
        galactic_b_dms = sky_coords.galactic.b

        #print(galactic_l_dms)
        #print(galactic_b_dms)

        #Convert from degree minute second to decimal degree then to radians
        #Space filler to catch blank space after second mark
        d_l, m_l, s_l, space_filler = re.split('[dms]', str(galactic_l_dms))
        d_b, m_b, s_b, space_filler = re.split('[dms]', str(galactic_b_dms))

        """
        print(d_l)
        print(m_l)
        print(s_l)
        print(d_b)
        print(m_b)
        print(s_b)
        """
        #Add the degree components and covert them to radians
        galactic_l = (float(d_l) + float(m_l)/60 + float(s_l)/(60*60)) * np.pi / 180
        galactic_b = (float(d_b) + float(m_b)/60 + float(s_b)/(60*60)) * np.pi / 180
        galactic_b = 9.621 * np.pi / 180
        galactic_l = 0.196 * np.pi / 180
        print(galactic_l * 180 / np.pi)
        print(galactic_b * 180 / np.pi)

        #With these, we can calculate the LSR velocity using
        #Sourced from https://www.atnf.csiro.au/people/Tobias.Westmeier/tools_hihelpers.php#restframes
        LSR_velocity = bary_velocity - (9 * np.cos(galactic_l) * np.cos(galactic_b) + 12 * np.sin(galactic_l) * np.cos(galactic_b) + 7 * np.sin(galactic_b))

        return LSR_velocity
    

    def plot_all(self, start_position = 1):
        #Will take in a start position, and then produce .spt for each on off pair
        #Will add time of previous scans to the epoch to properly correct veloity
        #This is the starting time.
        epoch_numpy = np.datetime64(self.epoch)

        #Convert from scna number to python list indicee
        i = start_position - 1
        #Use -1, as the off scans will be i + 1, and the last scan will be an off scan.
        #TODO Have option to cycle through if more than 12 sources are provided
        while i+1 < len(self.spec_data_list):
            #print(i)
            #Calculate the velocity correction
            RA_temp = self.observation_coordinates[i][0]
            Dec_temp = self.observation_coordinates[i][1]
            #print(RA_temp)
            #print(Dec_temp)

            vel_cor_temp = self.velocity_correction(self.lat, self.long, self.alt, epoch_numpy, RA_temp, Dec_temp)
            print(str(self.observation_names[i]) + " " + str(vel_cor_temp))

            #Write the .spt plot
            #self.write_plot(on_plot_position = i, off_plot_position = i+1, velocity_correction = vel_cor_temp)

            #Set the epoch to two scan times forwards (as an on and off have occured between previous and current)
            #on_off_time_delta = np.timedelta64(2 * self.observation_length, 's')
            #epoch_numpy = epoch_numpy + on_off_time_delta

            #Interate i to th next on scan
            i = i + 2


plot_maser = maser_plot()
plot_maser.plot_all()
