import subprocess
import argparse

#Author: Callum Macdonald
#Direct any questions to either Callum.Macdonald@utas.edu.au or CallumRMacdonald@hotmail.com
#If another author contributes new methods to the script, please add your name and contact here
"""
#Program designed to process all IFs and channels using SWSpec
#Inifile will have channels and name completely replaced
#Everything else in the file should be set as required
#The line numbers of the inifile are used for text replacement, so it should have the format that is specified at the bottom of the file.

#WARNING: At the completion of this script, it will delete all files in the directory of the form inifile.temp*
#Ensure that all permanent inifiles have more sensible names than .temp

#This program is optimised for processing involving more scans than there are cores, and if there are fewer scan than cores, it will
	#only processes the number of scans simultaneously. If this usage is required, then it might need to be re-written. 
	#Sorry for not doing this myself, but it was not pertinent to my work
#For example, if 2 scans are to be processed by 8 cores (there will be 32 total processes because of 2 * 4 IFs * 4 channels), 
	#it will only processes two items at once
"""

#TODO: Add other parameters that Guifre may want to adjust in other work (eg VDIX). Can be done manually currently.

class processAll:
    def __init__(self):
        #Collect information about what to process and how to do it using argparse
        parser = argparse.ArgumentParser()

        #Name of the inifile for the c and d IFs
        parser.add_argument(
            "-inicd",
            dest="inifile_name_cd",
            help="Name of the inifile for the c and d IFs", 
            type=str, 
            default='inifile.cd' #default set to base inifile name for c and d
        )

        #Name of the inifile for the e and f IFs
        parser.add_argument(
            "-inief",
            dest="inifile_name_ef",
            help="Name of the inifile for the e and f IFs", 
            type=str, 
            default='inifile.ef' #default set to base inifile name for e and f
        )

        #File name format
        #Option to allow different file name conventions to be used
        #Currently implemented options are 'fm' for Callum's flare monitoring, and 'sc' for Guifre's default format
        parser.add_argument(
            "-nf",
            dest="file_name_format",
            help="Format to use when creating file names from the given parameters. 'fm', 'sc', ...", 
            type=str, 
            default='fm' #default set to 'fm' because Callum wrote the script
        )

        #Observation name
        #Name of the observation being processed, for use in locating the data and writing the file name
        #Default set to quit program, as propper input must be given
        parser.add_argument(
            "-name",
            dest="observation_name",
            help="Name of the observation", 
            type=str, 
            default='not_set' #default set to quit the program
        )

        #Processing cores
        #The number of processes to run simultaneously
        #Increasing this number will decrease the processing time, up to the maximum cores of the machine the processing is occuring on
        parser.add_argument(
            "-cores",
            dest="cores",
            help="Number of processes to run simultaneously", 
            type=int, 
            default=2 #default set to 2, so as to not copletely slow down flexbuff machines. 4 seems to work fine, but will use more comuting power
        )

        #Station
        #The station the observations were made at
        parser.add_argument(
            "-st",
            dest="station",
            help="Code for the station the observations were made at", 
            type=str, 
            default='hb' #default set to Hb, the station that the flexbuff corresponds to
        )

        #Path
        #The path to the folder in which the raw data is mounted
        parser.add_argument(
            "-path",
            dest="path",
            help="The path to the folder with the raw observations mounted into it", 
            type=str, 
            default='not_set' #if default value is left, will be set o ../observation_name+station/ which will work for fm observations
        )

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

        #First scan
        #The number of the first scan to be processed.
        #All scans between this one and the last (inclusive) will be processed
        #TODO: Might be good to have an option to use a list of scans instead of a range
        parser.add_argument(
            "-first",
            dest="first_scan",
            help="The number of the first scan to be processed", 
            type=int, 
            default=1 #default set to 1 for use in the fm observations
        )

        #Last scan
        #The number of the last scan to be processed.
        #All scans between the first scan and this one(inclusive) will be processed
        parser.add_argument(
            "-last",
            dest="last_scan",
            help="The number of the last scan to be processed", 
            type=int, 
            default=24 #default set to 24 for use in the fm observations
        )



        #Extract the argeparse inputs into global variables
        args =parser.parse_args()
        self.inifile_name_cd = str(args.inifile_name_cd)
        self.inifile_name_ef = str(args.inifile_name_ef)
        self.file_name_format = str(args.file_name_format)
        self.observation_name = str(args.observation_name)
        self.path = str(args.path)
        self.station = str(args.station)
        self.IFs = args.IFs
        self.channels = args.channels
        self.first_scan = int(args.first_scan)
        self.last_scan = int(args.last_scan)
        self.cores = int(args.cores)

        #Error catchers
        if self.observation_name == 'not_set':
            print('Observation name not provided\n Program terminating')
            quit()
        
        if self.path == 'not_set':
            print('Path not provided\n Using ../' + str(self.observation_name) + str(self.station) + "/ for the path")
            self.path = '../' + str(self.observation_name) + str(self.station) +'/'

        if self.first_scan > self.last_scan:
            print('First scan number is larger than last scan number\n Program terminating')
            quit()


    def edit_inifile(self, inifile_base, channel_num, IF, scan_num):
        """
        Function that will take in an open base inifile andsome parameters
        Returns a list of edited lines (read/write lines)
        The name format is decided by the file_name_format -nf input
        """
        #Open the text as a list of lines, so that sspecific lines can be replaced as needed
        inifile_edited = inifile_base.readlines()

        #Replace the channel numbers in lines 10 and 11 (python counts from 0, sop list pos 9 and 10)
        inifile_edited[9] = 'UseFile1Channel   = ' + str(channel_num) + '\n'
        inifile_edited[10] = 'UseFile2Channel   = ' + str(channel_num) + '\n'

        #Replace the name of the output file
        if self.file_name_format == 'fm':
            #flare monitoring format
            inifile_edited[28] = 'BaseFilename1 = ' + self.observation_name + '_IF:' + IF + '_ch:' + str(channel_num) + '_st:' + self.station + '_' + str(scan_num) + '_%fftpoints%pt_%integrtime%s_VDIX\n'
        elif self.file_name_format == 'sc':
            #Guifre's format
            inifile_edited[28] = 'BaseFilename1 = ' + self.observation_name + '_' + self.station + '_VDIX_n' + IF + str(scan_num) + '_%fftpoints%pt_%integrtime%s_ch' + str(channel_num) + '\n'
        else:
            print('Invalid file name format\n Program terminating')
            quit()

        return inifile_edited

    def proc_one(self, channel_num, IF, scan_num):
        """
        Function that will process a single scan.
        Will not return a value, but will use SWSpec to write out the three default 
        outputs to files in the directory the script was called in
        """

        #Edit the temporary inifile to have the correct parameters
        #Open the inifile corresponding to the IF being processed
        if IF == 'c' or IF == 'd':
            inifile_base = open(self.inifile_name_cd)
        elif IF == 'e' or IF == 'f':
            inifile_base = open(self.inifile_name_ef)
        else:
            print('Invalid IF selected while attempting to open the inifile\n Program terminating')
            quit()
        
        #Use edit_inifile to create a temporary inifile
        inifile_edited_lines = self.edit_inifile(inifile_base, channel_num, IF, scan_num)

        #Create and open the new temporary inifile, then write out the edited lines
        #With this name it will be unique, allowing simultaneous processesing
        inifile_temp = open('inifile.temp' + str(IF) + str(channel_num) + str(scan_num), 'w')
        inifile_temp.writelines(inifile_edited_lines)
        inifile_temp.close()

        #Run SWSpectrometer using the new inifile and the data specified by the scan number
        #First, for clarity, make a string defining the path and file for the raw data
        #Needs to be in the format of '/path/fm014_hb_no0004_c'
        path_to_file = str(self.path) + str(self.observation_name) + "_" + str(self.station) + "_no" + str(scan_num) + "_" + str(IF)
        
        #The following line of code is based upon Code originally written by Guifre
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        process = subprocess.Popen(['swspectrometer','inifile.temp' + str(IF) + str(channel_num) + str(scan_num), path_to_file])
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        return process

    
    def proc_all(self):
        '''
        Will cycle through all IFs, and for each IF through all channels, and for each channel all scans.
        For each scan, process_one will be called with the IF, channel and scan number
        '''
        #Cycle through all IFs to be processed
        for IF in self.IFs:
            #Cycle through all channels to be processed (for each IF)
            for channel in self.channels:
                #Set the scan position tracker to the first scan
                current_scan_possition = self.first_scan
                #Continue to process scans until the scan position is equal to the last scan
                while current_scan_possition <= self.last_scan:
                    #Creat a list to store the simultaneous processes in
                    processes = []
                    #Use a for loop to set up 'cores' number of simultaneous processes
                    for core in range(self.cores):
                        #Check that the current scan position has not exceeded the final scan (only an issue if total number of scans is not divisible by cores)
                        if current_scan_possition <= self.last_scan:
                            #Convert the current scan position to the #### format used by telescopes and SWSpec
                            #TODO: Use fstring to do this much quicker
                            if current_scan_possition < 10:
                                scan_number = '000' + str(current_scan_possition)
                            elif current_scan_possition < 100:
                                scan_number = '00' + str(current_scan_possition)
                            elif current_scan_possition < 1000:
                                scan_number = '0' + str(current_scan_possition)
                            else:
                                scan_number = str(current_scan_possition)
                            
                            #Use proc_one to processes the scan corresponding to the current parameters
                            print('Now processesing IF: ' + str(IF) + ', channel: ' + str(channel) + ', scan number: ' + str(scan_number))
                            current_process = self.proc_one(channel, IF, scan_number)
                            processes.append(current_process)

                            #Move the counter to the next scan
                            current_scan_possition += 1
                        #If current_scan_position exceeds last scan pos during the loop, then the program will not be operating at maximum efficiency
                        #This little message alerts the user when this occurs
                        else:
                            print('Cores is not an integer factor of the total scans\n There are currently cores idle, consider changing number of cores in future processing runs')


                    #Now wait for these processes to complete, at which stage the cores will be freed and the while loop can continue
                    for proc in processes:
                        proc.wait()

        #Loops done, now clean up inifiles
        clean_inifiles = subprocess.Popen(['rm inifile.temp*'], shell=True)
        clean_inifiles.wait()
                        

processor = processAll()
processor.proc_all()


#Below is the text content of the inifiles in the correct format for this script
#The FFTpoints will need to be set depending on what resolution is required
#(For fm observations, ef FFTpoints will be approximately half that of cd for ~6 and ~12 GHz)
"""
[Spectrometer]
FFTpoints = 32000

SourceFormat      = VDIX
SourceChannels    = 4
BitsPerSample     = 2
BandwidthHz       = 32e6
SourceSkipSeconds = 0

UseFile1Channel   = 1
UseFile2Channel   = 1

ChannelOrderIncreasing = True

PlotProgress = No
DoCrossPolarization = No
MaxSourceBufferMB = 256
FFTIntegrationTimeSec = 540

NumCPUCores = 4
WindowType = Hanning
FFToverlapFactor = 2

SinkFormat = Binary

PCalOffsetHz = 10e3
ExtractPCal = No

BaseFilename1 = none
BaseFilename2 = none
"""
