## Simple Visual Plot

## Working on making the visual plots look simpler and all at once, not each individual orbit 

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Polygon, Wedge
import numpy as np
import pandas as pd
import datetime
from ovationpyme import ovation_utilities
import aacgmv2
from ovationpyme.ovation_prime import FluxEstimator

TARGET = 'AT2021sdu'
OBSID = '4202470102'
DATE = '12-02-2021'


plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'monospace'
plt.rcParams['legend.fontsize'] = 7

## Setting up the cusp region ##
RE = 1
REkm = 6371
plot_limit = RE * 3
#circ_RE = 2*np.pi
#cuspLAT = 75 #degrees in lat
#cuspMLT = 12 #magnetic local time in hours
#cuspLATwidth = 5
#cuspMLTwidth = 1.5
#cusp_long = ((cuspMLT % 24) / 24) * 360
#lat_rad = np.radians(cuspLAT)
#lon_rad = np.radians(cusp_long)
x_cusp_gsm = (0.3+0.1) / 2
y_cusp_gsm = (-0.3+0.1) / 2
#x_cusp_gsm = RE * np.cos(lat_rad) * np.cos(lon_rad)
#y_cusp_gsm = RE * np.cos(lat_rad) * np.sin(lon_rad)
cuspLATwidthRE = 0.1 - (-0.3)
cuspMLTwidthRE = 0.3 - 0.1
#cuspLATwidthRE = (cuspLATwidth / 360) * circ_RE
#cuspMLTwidthRE = ((cuspMLTwidth * 15) / 360) * circ_RE


## Getting ISS and NICER pointing data ##
full_data = pd.read_csv(f'C:/Users/aduja/NICER/{OBSID}_FullDataValues.csv')
iss_x_gsm = full_data['ISS X_GSM'] / REkm
iss_y_gsm = full_data['ISS Y_GSM'] / REkm
nicer_x_point = full_data['NICER X_Point']# [::35] # every 20th vector plotted
nicer_y_point = full_data['NICER Y_Point'] #[::35]
unix_times = full_data['time_Unix']
full_data['datetime'] = pd.to_datetime(unix_times, unit='s')

#iss_x_20th = iss_x_gsm[::35] # the iss location every 20th nicer point vector
#iss_y_20th = iss_y_gsm[::35]

## Calculating the auroral flux for the diffuse aurora and energy flux
flux_estimator = FluxEstimator(atype='diff', energy_or_number='energy')
flux_values = []
for dt in full_data['datetime']:
    mlat, mlt, flux = flux_estimator.get_flux_for_time(dt, hemi='N')
    flux_values.append(flux)

full_data['flux'] = flux_values

#print(type(full_data['flux'][0]))
#print(full_data['flux'][0].shape)
# ^^ This is giving a 2D array that includes the magnetic local time index and the magnetic latitude index
# need to convert to 1D in order to test the plot and make sure everything works. 

## Plotting to Visualize the 2D array
# plt.figure(figsize=(10, 8))
# plt.imshow(flux_values[0], cmap='viridis', aspect='auto')
# plt.colorbar(label='Flux')
# plt.title(f'Target = {TARGET} OBSID = {OBSID} {DATE} \n Auroral Flux Heatmap at First Timestamp')
# plt.xlabel('Magnetic Local Time Index')
# plt.ylabel('Magnetic Latitude Index')
# plt.show()


## Calc the mean of each 2D array in flux_values
flux_means = [np.nanmean(flux) for flux in flux_values]
full_data['flux_mean'] = flux_means

## Practicing the auroral flux plot ##
# plt.figure(figsize=(10, 5))
# plt.plot(full_data['datetime'], full_data['flux_mean'], label='Mean Auroral Flux', color='orange')
# plt.xlabel('Time')
# plt.ylabel('Mean Flux (mW/m^2)')
# plt.title(f'Target = {TARGET} OBSID = {OBSID} {DATE} \n Mean Auroral Flux')
# plt.legend()
# plt.show()

flux_threshold = 0.5
mlat_range = np.linspace(50, 80, 80)
mlt_range = np.linspace(0, 24, 96) 

def mlt_mlat_to_geo(mlt, mlat, date):
    #convert mlt to mag long
    magnetic_longitude = (mlt / 24.0)*360.0-180.0
    #convert mag lat, lon to geo coords at 100 km altitude
    geo_lat, geo_lon, _ = aacgmv2.convert_latlon(mlat, magnetic_longitude, 100, date, a2g=True)
    x = geo_lon
    y = geo_lat
    return x, y

fig, ax = plt.subplots()
date = datetime.datetime.strptime(DATE, "%m-%d-%Y")

for flux_array in flux_values:
    boundary_points = []
    for i, mlat in enumerate(mlat_range):
        for j, mlt in enumerate(mlt_range):
            if flux_array[i, j] >= flux_threshold:
                x, y = mlt_mlat_to_geo(mlt, mlat, date)
                boundary_points.append((x, y))
    
    if boundary_points:
        boundary_polygon = Polygon(boundary_points, closed=True, edgecolor='red', fill=None)
        ax.add_patch(boundary_polygon)

ax.plot(iss_x_gsm, iss_y_gsm, color='black', linestyle='-', marker='None', label='ISS Orbit', linewidth=0.7)
ax.set_xlim(-plot_limit, plot_limit)
ax.set_ylim(-plot_limit, plot_limit)
ax.legend()

plt.show()

## Breaking the points for separate orbits
x_diff = np.abs(np.diff(iss_x_gsm))
threshold = 150 / REkm # This is the threshold where the next orbit starts ... 150 km 
# ^^^ Will need to check this threshold for other obs... but for 4202470102 & 103, this satisfies
#x_diff = np.insert(x_diff,0,0) # inserts a 0 at the beginning to align with the data length
#iss_x_gsm[x_diff > threshold] = np.nan
#iss_y_gsm[x_diff > threshold] = np.nan
orbit_start = np.where(x_diff > threshold)[0]+1

orbit_indices = [0] +list(orbit_start) + [len(iss_x_gsm)]

for idx in orbit_start:
    iss_x_gsm[idx] = np.nan
    iss_y_gsm[idx] = np.nan

## Plotting
fix, ax = plt.subplots()

e_outline = plt.Circle((0,0), RE, color='blue', fill=False)
ax.add_artist(e_outline)

cusp = Ellipse((x_cusp_gsm,y_cusp_gsm), width=cuspMLTwidthRE, height=cuspLATwidthRE, angle=0, edgecolor='black', facecolor='brown', label='Approx. Static Cusp Location', zorder=3)
ax.add_patch(cusp)
ax.set_aspect('equal')

ax.plot(iss_x_gsm, iss_y_gsm, color = 'black', linestyle='-', marker='None', label='ISS Orbit', linewidth=0.7)
#ax.quiver(iss_x_20th, iss_y_20th, nicer_x_point, nicer_y_point, color='orange', scale=1, scale_units='xy', angles='xy', linewidth=0.5)

line_length = plot_limit*2.5  #set line length for vectors to end at plot edges

##plot 10 NICER pointing directions per orbit
for start, end in zip(orbit_indices[:-1], orbit_indices[1:]):
    num_points = end - start
    step = max(num_points // 7, 1) #7 pointing vectors
    selected_indices = range(start, end, step)[:7]  #select up to 7 indices

    #plot NICER pointing directions as lines
    for i in selected_indices:
        x, y, dx, dy = iss_x_gsm[i], iss_y_gsm[i], nicer_x_point[i], nicer_y_point[i]
        ax.plot([x, x + line_length * dx], [y, y + line_length * dy], color='orange', linestyle='-', linewidth=0.5)



night_side = Wedge(center=(0,0), r=RE, theta1=90, theta2=270, color='black', alpha=0.2)
ax.add_patch(night_side)

ax.set_xlim(-plot_limit,plot_limit)
ax.set_ylim(-plot_limit,plot_limit)
ax.invert_xaxis()
ax.invert_yaxis()

major_ticks = [-3,-2,-1,0,1,2,3]
ax.set_xticks(major_ticks)
ax.set_yticks(major_ticks)
minor_ticks = (np.arange(-plot_limit, plot_limit + 0.1, 0.2))
ax.tick_params(axis='both', which='both', direction='in', top=True, right=True,labelsize=8)
ax.set_xticks(minor_ticks, minor=True)
ax.set_yticks(minor_ticks, minor=True)


ax.set_xlabel('GSM X (RE)', fontsize=9)
ax.set_ylabel('GSM Y (RE)', fontsize=9)
ax.set_title(f'Target = {TARGET} OBSID = {OBSID} \n {DATE}', fontsize=9)
ax.legend(loc='best') #upper right
plt.show()


