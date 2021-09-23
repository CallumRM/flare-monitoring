import numpy as np
import matplotlib.pyplot as plt
import os, argparse, re
import errno
import sys

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
            default=32000 #default set to 32000 which is the number being used for the current flair monitoring setup
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
            default='hb'
        )

        #Velocity range
        #The range of velocities to generate the plot for. First value is lower, second is upper limit
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
            default=['G328.237-0.547', 'off']
            #need to re-write this, some chars not accepted   
        )
        #Source velocities
        #A list of the recession velocities of the sources observed in Km/s
        #Default set for fm observations of 6 sources, velocities sourced from multibeam catalogues
        #TODO: Finish this please
        parser.add_argument(
            "-velocities",
            dest="source_velocities",
            nargs="+",
            #Default set for fm observations of 6 sources
            
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
        #Converted into decimal degrees
        if args.observatory == 'hb':
            self.lat = -42.80557
            self.long = 147.4381
            self.alt = 40.967
        #TODO: FInish these for other observatories
        """
        elif args.observatory == 'ke':
            self.lat =
            self.long = 
            self.alt = 
        """
        #If it is not a name, and instead the three pieces of coord data have been provided, use them
        elif len(args.observatory) == 3:
            self.lat = args.observatory[0]
            self.long = args.observatory[1]
            self.alt = args.observatory[2]
        #If the input is not specified correctly, quit the program
        #(If velocity doesn't matter, and hence no observatory is required, just use hb)
        else:
            print('Observatory data entered incorrectly using -obs. \n Either enter lat lon alt or the abreviated observatory name.\n Script terminating.')
            quit


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
    #epoch in string fits format'YYYY-MM-DDTHH:MM:SS'
    #RA string in **h**m**.*s format
    #Dec string in ***d**m**.*s
    #Will return a float of the velocity correction
    #TODO: Set up a HMS to Degree converter, so that everything can be entered in the standard coords
    def velocity_correction(lat, long, alt, epoch, RA, Dec):
        #First step is to generate the Heliocentric velocity using PyAstronomy
        
        return LSR_vel
    
    #Method to produce a plot of a spctra
    #If the plot position are left unspecified, they will be set to the input spec position.
    #For the plot_all method, they will be set as required.
    #Will wrtie out a file
    def write_plot(self, single_plot_position = None, on_plot_position = None, off_plot_position = None):
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
        if self.plot_option == 1:
            #Single spectrum plot
            print("Plot option 1 selected, writing out single spectrum")
            print("Spectrum file position in list is: " + str(single_plot_position + 1))

            plot = open(self.observation_names[single_plot_position] + ".spt", "x")
            #First, inster the name and axis label components
            plot.write(" heading \n " + self.observation_names[single_plot_position] +"\n xlabel \n Velocity w.r.t. LSR (km/s) \n ylabel \n Flux (Jy) \n ")

            #Calculate the plotting window
            """
            #Calculate the velocity 
            #TODO: Add function to move time forwards across observations, current system will only be reliable for first 1 or 2 obs
            time_obv = np.datetime64(self.epoch)
            #Convert the epoch into julian time for use with barycorrpy
            julian_time = time_obs.jd
            """


            """
            #Now specifiy the plotting window
            vel_source = source_velocities[single_spec_pos]
            vel_lower = vel_source + self.min_plot_vel
            vel_
            """
            #TODO:Currently have it set up to enter the range manually, update this 
            adjusted_vel = self.general_vel #Need to fix this
            vel_lower = self.min_plot_vel
            vel_upper = self.max_plot_vel

            #Need an input for the plot y scale, currently hard set
            flux_min = -1
            flux_max = 20

            #,spt form needs data \n 'number of points'
            #Run through the list of velocities, if they are within the range, add to the count
            #TODO: Come up with a more efficient way to do this (maybe calculate the number of points base on vel)
            data_point_count = 0
            for i in range(len(adjusted_vel)):
                if (adjusted_vel[i]>vel_lower and adjusted_vel[i]<vel_upper):
                    data_point_count += 1
            
            #Now write these sections into the file
            plot.write("baxis \n " + str(vel_lower) + " \n " + str(vel_upper) + " \n taxis \n " + str(data_point_count) + " \n 1 \n laxis \n " + str(flux_min) + " \n " + str(flux_max) + " \n data \n " + str(data_point_count) + " \n ")
            
            print(single_plot_position)
            print(self.spec_data_list)
            print(len(self.spec_data_list[single_plot_position]))
            print(self.spec_data_list[single_plot_position])
            

            #Now the data needs to be written into the file
            for i in range(len(adjusted_vel)):
                if (adjusted_vel[i]>vel_lower and adjusted_vel[i]<vel_upper):
                    print(i)
                    plot.write(str(adjusted_vel[i]) + "    " + str(self.spec_data_list[single_plot_position][i]) + " \n ")

            #Now add an end comment and close the file
            plot.write("end")
            plot.close()
            #Should be everything

        if self.plot_option == 2:
            print("Plot option 2 selected, on off plot" + str(on_plot_position) + " " +str(off_plot_position))
            plot = open(self.observation_names[single_spec_pos] + ".spt", "x")
            #First, insert the name and axis label components
            plot.write(" heading \n " + self.observation_names[on_plot_position] +" On Off \n xlabel \n Velocity w.r.t. LSR (km/s) \n ylabel \n Flux (Jy) \n ")

            #Calculate the plotting window
            """
            #Calculate the velocity 
            #TODO: Add function to move time forwards across observations, current system will only be reliable for first 1 or 2 obs
            time_obv = np.datetime64(self.epoch)
            #Convert the epoch into julian time for use with barycorrpy
            julian_time = time_obs.jd
            """


            """
            #Now specifiy the plotting window
            vel_source = source_velocities[single_spec_pos]
            vel_lower = vel_source + self.min_plot_vel
            vel_
            """
            #TODO:Currently have it set up to enter the range manually, update this 
            adjusted_vel = self.general_vel #Need to fix this
            vel_lower = self.min_plot_vel
            vel_upper = self.max_plot_vel

            #Need an input for the plot y scale, currently hard set
            flux_min = -1
            flux_max = 20

            #,spt form needs data \n 'number of points'
            #Run through the list of velocities, if they are within the range, add to the count
            #TODO: Come up with a more efficient way to do this (maybe calculate the number of points base on vel)
            data_point_count = 0
            for i in range(len(adjusted_vel)):
                if (adjusted_vel[i]>vel_lower and adjusted_vel[i]<vel_upper):
                    data_point_count += 1
            
            #Now write these sections into the file
            plot.write("baxis \n " + str(vel_lower) + " \n " + str(vel_upper) + " \n taxis \n " + str(data_point_count) + " \n 1 \n laxis \n " + str(flux_min) + " \n " + str(flux_max) + " \n data \n " + str(data_point_count) + " \n ")
            
            #Calculate the quotient between the on and off plots,
            print(self.spec_data_list[on_plot_position])
            print(self.spec_data_list[off_plot_position])
            quotient_array = np.zeros(len(self.spec_data_list[on_plot_position])) #Use on position for size as it must be within list length anyway
            print(quotient_array)
            for i in range(len(quotient_array)):
                #divide the ith value from the on data by the ith value from the on data, then subtract 1 (quotient formula)
                quotient_array[i] = self.spec_data_list[on_plot_position][i] / self.spec_data_list[off_plot_position][i] - 1 

            print(adjusted_vel)
            print(quotient_array)
            #Now the data needs to be written into the file
            for i in range(len(adjusted_vel)):
                if (adjusted_vel[i]>vel_lower and adjusted_vel[i]<vel_upper):
                    plot.write(str(adjusted_vel[i]) + "    " + str(quotient_array[i]) + " \n ")

            #Now add an end comment and close the file
            plot.write("end")
            plot.close()
            #Should be everything
    #def plot_all(self):
        #Need to be able to combine plots (double feature, position + (len/2 + pos))

#TODO Work out what I'm doing with the names
"""
    #We now have a list of file names and a corresponding list of file paths.
    #If no observation name was provided, set it to the file name
    if observation_name == None:
        if self.plot_option == 2:
            #If it is an on off plot, use the on file
            observation_name = file_name_list[on_position]
        else:
            #If using 1 or 3, use the single source file name
            observation_name = file_name_list[single_spec_position]
"""

#Need to convert this to splot output.
"""
    #Now plot the data, using the selected option
    if plot_option == 1:
        #Single spectrum plot
        print("Plot option 1 selected, plotting single spectrum")
        print("Spectrum file position in list is: " + str(single_spec_position + 1))

        plt.figure()
        plt.plot(vel, spec_data_list[single_spec_position], color = "blue")
        plt.xlabel("Recession Velocity Km/s")
        plt.title(observation_name)
        plt.xlim(min_plot_vel, max_plot_vel)  #Look at centre 400 Km/s, should have everything in here

        if show_spec == True:
            plt.show()
        else:
            plt.savefig(fname=(observation_name + ".svg"), bbox_inches="tight")

    elif plot_option == 2:
        #On-off plot
        print("Plot option 2 selected, plotting the quotient of the on and off file")
        print("On file position in list is: " + str(on_position + 1))
        print("Off file position in list is: " + str(off_position + 1))
        #Need to generate the quotient of the on and off data
        #Make a new array to hold the quotient output
        quotient_array = np.zeros(len(spec_data_list[on_position])) #Use on position for size as it must be within list length anyway
        for i in range(len(quotient_array)):
            #divide the ith value from the on data by the ith value from the on data, then subtract 1 (quotient formula)
            quotient_array[i] = spec_data_list[on_position][i] / spec_data_list[off_position][i] - 1 

        plt.figure()

        plt.plot(vel, quotient_array, color = "blue")
        plt.xlabel("Recession Velocity Km/s")
        plt.title(observation_name)
        plt.xlim(min_plot_vel, max_plot_vel)  #Look at centre 400 Km/s, should have everything in here

        if show_spec == True:
            plt.show()
        else:
            plt.savefig(fname=(observation_name + ".svg"), bbox_inches="tight")

    elif plot_option == 3:
        #Averaged baseline plot
        print("Plot option 3 selected, plotting the quotient of the selected spectrum and the avergae of the other spectra")
        print("Spectrum file position in list is: " + str(single_spec_position + 1))

        #Make a new array to hold the average values of the other spectra
        average_array = np.zeros(len(spec_data_list[on_position])) #Use on position for size as it will be within the list
        #For each element in the array, take the average value from the spectra other than the selected spectrum
        for i in range(len(average_array)):
            for j in range(len(spec_data_list)):
                if j != single_spec_position: #All but the selected position
                    average_array[i] += spec_data_list[j][i] #Add the ith element of the jth spectra to the ith position in the average
        #Now we have the sum of all the other spectra, to get average divide every element in the array by the number of spectra used
        average_array = average_array / (len(spec_data_list) - 1) #Subtract 1 as the list contains the selected spectrum which was not avergae over


        #Make a new array to hold the quotient output
        quotient_array = np.zeros(len(spec_data_list[on_position])) #Use on position for size as it must be within list length anyway
        for i in range(len(quotient_array)):
            #divide the ith value from the on data by the ith value from the on data, then subtract 1 (quotient formula)
            quotient_array[i] = spec_data_list[on_position][i] / average_array[i] - 1 

        #Plot the spectrum
        plt.figure()

        plt.plot(vel, quotient_array, color = "blue")
        plt.xlabel("Recession Velocity Km/s")
        plt.title(observation_name)
        plt.xlim(min_plot_vel, max_plot_vel) #Look at centre 400 Km/s, should have everything in here

        if show_spec == True:
            plt.show()
        else:
            plt.savefig(fname=(observation_name + ".svg"), bbox_inches="tight")
"""

plot_maser = maser_plot()
plot_maser.write_plot()