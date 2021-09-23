import os, argparse, re
import subprocess
import sys
import astropy.coordinates
import numpy as np
import pandas as pd
import scipy.fft
import matplotlib.pyplot as plt

#test = subprocess.Popen(['rm *.vex'], shell=True)
#test.wait()
'''
channel_num = 2
ini_edit = open('inifile.Hb').readlines()
ini_edit[9] = 'Test\n'
ini_edit[10] = 'Test_line_two\n'

ini_output = open('inifile.Hb' + str(channel_num), 'w')
ini_output.writelines(ini_edit)
ini_output.close()
'''
#sky_coords_1 = astropy.coordinates.SkyCoord('18h06m14.67s', '-80d31m32.4s', frame = 'icrs')
#sky_coords_2 = astropy.coordinates.SkyCoord('19h06m14.67s', '-80d31m32.4s', frame = 'icrs')

#Works
#dif = sky_coords_1.spherical_offsets_to(sky_coords_2)

#dif = sky_coords_1.separation(sky_coords_2)

#print(dif)

#test_name = re.split('[_:_:_:_____]', 'fm024_IF:c_ch:2_st:hb_0003_64000pt_540s_VDIX_swspec')
#print(test_name)


#CSV tests and fourier


parser = argparse.ArgumentParser()
parser.add_argument(
    "-csv",
    dest="csv",
    nargs="+",
    default=[] #This argument obviously needs to be passed in, and if it isn't then an error will be returned
    )
args = parser.parse_args()

csv_files = args.csv


#print(csv_file)

#csv_data = np.genfromtxt(csv_file, delimiter=',', dtype=None)

#print(csv_data[1][0])
#print(csv_data[1][1])
velocities=[]
intensities=[]
for i in range(len(csv_files)):
    csv_pandas_df = pd.read_csv(csv_files[i], sep=',', skiprows=1)
    velocities.append(csv_pandas_df.loc[:,'vel(km/s)'].to_numpy())
    intensities.append(csv_pandas_df.loc[:,'intensity(Jy)'].to_numpy())


average_intensity = np.zeros(len(intensities[0]))
for i in range(len(intensities)):
    average_intensity += intensities[i]

average_intensity = average_intensity / len(intensities)

plt.figure()

plt.plot(velocities[0], average_intensity, color='red')
#plt.plot(velocities[0], intensities[0], color='green')
plt.show()


