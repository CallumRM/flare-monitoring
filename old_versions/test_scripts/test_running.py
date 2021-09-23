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

csv_file = args.csv[0]


#print(csv_file)

#csv_data = np.genfromtxt(csv_file, delimiter=',', dtype=None)

#print(csv_data[1][0])
#print(csv_data[1][1])

csv_pandas_df = pd.read_csv(csv_file, sep=',', skiprows=1)

#print(csv_pandas_df)
velocity = csv_pandas_df.loc[:,'vel(km/s)'].to_numpy()
intensity = csv_pandas_df.loc[:,'intensity(Jy)'].to_numpy()

fft_intensity = scipy.fft.fft(intensity)
#print(fft_intensity)

ifft_intensity = scipy.fft.ifft(fft_intensity)

#print(len(intensity))
#print(ifft_intensity[1])
sin_from_vel = np.sin(2*np.pi * velocity / 10)

fft_sin = scipy.fft.fft(sin_from_vel)
dct_sin = scipy.fft.dct(sin_from_vel)

#print(fft_sin[30])
#print(dct_sin)

#test_array = np.zeros(len(sin_from_vel), dtype='complex_')
#test_array[29] = fft_sin[29]
#test_array[30] = fft_sin[30]
#test_array[31] = fft_sin[31]

#print(velocity[30])
#print(len(velocity))
vel_new = np.arange(velocity[-1], velocity[0], 0.01)
vel_count_shift = np.arange(0, len(velocity), 1) + 0.5
vel_count = np.arange(0.5, len(velocity)+0.5, 1)
#print(2*np.pi*velocity[30]*vel_new)
#print(fft_sin[30])
sin_from_fft = 2*abs(fft_sin[30])*np.cos(vel_count * 2*np.pi* 30 / (len(velocity))) / len(velocity)
sin_from_fft_shift = 2*abs(fft_sin[30])*np.cos(vel_count_shift * 2*np.pi* 30 / (len(velocity))) / len(velocity)
array=[sin_from_fft, sin_from_fft, sin_from_fft]
sum_of_things = np.zeros(len(sin_from_fft))
for i in range(len(array)):
    sum_of_things += array[i]
#sin_from_fft = scipy.fft.ifft(test_array)

#Do the dct sum
sin_from_dct = np.zeros(len(vel_count))
for i in range(len(dct_sin)):
    sin_from_dct += dct_sin[i]*np.cos(vel_count * np.pi* i / (len(velocity))) / len(velocity)

sin_from_dct_shift = np.zeros(len(vel_count))
for i in range(len(dct_sin)):
    sin_from_dct_shift += dct_sin[i]*np.cos(vel_count_shift * np.pi* i / (len(velocity))) / len(velocity)


range_1 = list(range(50,550,1))
range_2 = list(range(51,551,1))

"""
plt.figure()

plt.plot(vel_count, sin_from_vel, color = 'blue')
#plt.plot(vel_count, sin_from_fft, color = 'red')
plt.plot(vel_count, sin_from_dct, color = 'purple')
plt.plot(vel_count_shift, sin_from_dct_shift, color = 'green')
#plt.xlim(-15,-10)

plt.show()
"""

#Try this but for actual data
intensity_flip = np.flip(intensity)
dct_intensity = scipy.fft.dct(intensity, type=3)
dct_intensity_flip = scipy.fft.dct(intensity_flip, type=3)

vel_count_shift = np.arange(0, len(velocity), 1) + 0.5
vel_count = np.arange(0, len(velocity), 1)

#Set up a veloicty shift array based on how much final velocity will need to be shifted.
#Need to add the value of the that will put the term closest to zero to zero divided by the step size
#([0,1,2,..]+vsmall/vstep)*vstep = [0,1,2...]*vstep + vsmall
smallest_positive_vel = 1000000 #Large so that first encountered is smaller, then everything else will be compared to that
for i in range(len(velocity)):
    if velocity[i] < smallest_positive_vel and velocity[i] > 0:
        smallest_positive_vel = velocity[i]

#print(smallest_positive_vel)
#print(velocity[1]-velocity[2])
#print(velocity[-2]-velocity[-1])

velocity_step = velocity[0] - velocity[1]
vel_count_shift = np.arange(0, len(velocity), 1) + (smallest_positive_vel / velocity_step) + 150
vel_shift = vel_count_shift * velocity_step + velocity[-1]

#This above stuff doesn't work, try velocity-smallest, then subtract and divide to get 0 to 666. That way, reverse gives perfect 0.
#print((velocity - velocity[-1]) / velocity_step)

intensity_from_dct = np.zeros(len(vel_count))
for i in range(len(dct_intensity)):
    intensity_from_dct += dct_intensity[i]*np.cos(vel_count * np.pi* i / (len(vel_count))) / len(dct_intensity)

intensity_from_dct_formula = np.zeros(len(vel_count))
for i in range(len(dct_intensity)):
    intensity_from_dct_formula += dct_intensity[i]*np.cos(vel_count * np.pi* (2*i +1) / (2*len(vel_count))) / len(dct_intensity)

intensity_from_dct_formula_shift = np.zeros(len(vel_shift))
for i in range(len(dct_intensity)):
    intensity_from_dct_formula_shift += dct_intensity[i]*np.cos(vel_count_shift * np.pi* (2*i +1) / (2*len(vel_count))) / len(dct_intensity)

intensity_from_dct_formula_shift_flip = np.zeros(len(vel_shift))
for i in range(len(dct_intensity)):
    intensity_from_dct_formula_shift_flip += dct_intensity_flip[i]*np.cos(vel_count_shift * np.pi* (2*i +1) / (2*len(vel_count))) / len(dct_intensity)

plot_vel = np.flip(velocity)
#vel_shift = np.flip(vel_shift)

plt.figure()

plt.plot(velocity, intensity, color = 'blue')
#plt.plot(vel_count, intensity_from_dct, color = 'red')
#plt.plot(plot_vel, intensity_from_dct_formula, color = 'green')
#plt.plot(vel_shift, intensity_from_dct_formula_shift, color = 'red')
#plt.plot(vel_shift, intensity_from_dct_formula_shift_flip, color = 'green')

plt.show()





"""
csv_pandas_df = pd.read_csv(csv_file, sep=' ', header=None)

print(csv_pandas_df)

"""






