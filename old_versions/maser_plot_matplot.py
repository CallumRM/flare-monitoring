import numpy as np
#from pysctrack import handler, core
import matplotlib.pyplot as plt
import os, argparse, re
import errno
import sys

"""
A script that will take some parameters and a list of spectra .bin files, and then producea plot
TODO: DOUBLE CHECK THIS IS TRUE
Will either produce a simple on-off plot, where the on and off spectra are selected from the list
    or, will produce an spectra normalised by the average of the other spectra
    or, will just plot a single spectra.

"""

#Use argparse to accept the input arguments
#Will need the bandwidth in MHz, number of channels, plotting option variables, and the spectra files.

#Parser will hold the argument data, will be extracted into other variables below.
parser = argparse.ArgumentParser()

#Bandwidth argument
#To pass in the bandwidth, type -bw 'bandwidth' after the script name, replacing 'bandwidth'
    #with the bandwidth in MHz
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
    default=32000 #default set to 32000 which is the number being used for the flair monitoring setup
)

#List of spectra argument
#Will accept a single spectrum or a list of spectra
#If a single spectrum is to be plotted, only one spectrum needs to be passed in
#For a simple on-off, both the on and off spectra must be passed in
#For the on averaged-off, at least two spectra need to be passed in (for sensible result, need at least ~5)
#TODO: Check that just one argument works (might not with list[0], instead wants list)
parser.add_argument(
    "-spec",
    dest="spec_list",
    nargs="+",
    default=[] #This argument obviously needs to be passed in, and if it isn't then an error will be returned
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
    default=(1,1) #Default to 1 and 1, as there will always be at least one spectrum
)

#End of input arguments
#Now extract variables out of the parser for simplicity
args =parser.parse_args() #args now is a list of all of the input variables
bandwidth = args.bandwidth * 1e6 #Convert bandwidth to Hz
num_channels = args.num_channels
spec_list = args.spec_list 
plot_option = args.plot_option
single_spec_position = args.single_spec_position - 1 #Subract 1 as ptyhon lists start at 0
on_position = args.on_off_positions[0] - 1 #Get the on (first) value from the tuple, subract 1 as ptyhon lists start at 0
off_position = args.on_off_positions[1] - 1 #Get the off (second) value from the tuple, subract 1 as ptyhon lists start at 0


#TODO: Add error catchers for all of the other values (ie plot option, spec_pos out of bounds)
#If the spec_list is empty, there is no data to be plotted, and the program quits
if len(spec_list) < 1:
    print('Must have at least one spectrum passed using -spec')
    print('Script terminating')
    quit



#First, check if the input .bin spectra files exist, and create a list of file names and a list of the paths to these spectra
#The splitting of name and path may not be necessary, but I don't know enough about opening files to make a call.
file_name_list = []
path_list = []
for i in range(len(spec_list)):
    if os.path.exists(spec_list[i]):
        path, file_name = os.path.split(spec_list[i]) #split the path and file name
        file_name_list.append(file_name)
        path_list.append("{}/".format(os.path.abspath(path)))
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), spec_list[i])

#We now have a list of file names and a corresponding list of file paths.

#Now calculate some variables and lists to be used for the plot
#TODO: Work out exactl what these variables are, so as to better document the code. Currently this is largely ripped
SR = int(2 * bandwidth)
df = float(SR / num_channels) #Differential freq unit
Nfft = int((num_channels / 2) + 1)
jf = np.arange(Nfft, dtype=int)
ff = df * jf # 
#TODO: Make x as velocity rather than frequency

#Make a list of the spectra data from the .bin files
spec_data_list = []
for i in range(len(spec_list)): #For each file name passed in
    file_temp = open(f"{path_list[i]}{file_name_list[i]}", "rb") #Open the .bin file
    temp_data = np.fromfile(file=file_temp, dtype="float32", count=Nfft) #Extract the data from the file
    spec_data_list.append(temp_data) #Put array in the list
    file_temp.close #Close the file


#Now plot the data, using the selected option
if plot_option == 1:
    #Single spectrum plot
    print("Plot option 1 selected, plotting single spectrum")

    plt.figure()
    plt.plot(ff, spec_data_list[single_spec_position], color = "blue")
    plt.xlabel("Frequency")
    plt.title(str(file_name_list[single_spec_position]))

    plt.show()

elif plot_option == 2:
    #On-off plot
    print("Plot option 2 selected, plotting the quotient of the on and off file")
    print("On file position in list is: " + str(on_position + 1))
    print("Off file position in list is: " + str(off_position + 1))
    #Need to generate the quotient of the on and off data
    #Make a new array to hold the quotient output
    quotient_array = np.zeros(len(spec_data_list[on_position])) #Use on position for size as it must be within list length anyway
    for i in range(len(quotient_array)):
        quotient_array[i] = spec_data_list[on_position][i] / spec_data_list[off_position][i] - 1 #divide the ith value from the on data by the ith value from the on data, then subtract 1 (quotient formula)
    
    plt.figure()

    plt.plot(ff, quotient_array, color = "blue")
    plt.xlabel("Frequency")
    plt.title("On off for " + str(file_name_list[on_position]) + " and " + str(file_name_list[off_position]))

    plt.show()


