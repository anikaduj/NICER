## Calculating the Auroral FLux for Obs: 4202470102 & 4202470103

import os
import sys
sys.path.append('C:/Users/aduja/NICER/OvationPyme/ovationpyme')
from ovationpyme import ovation_prime
import matplotlib.pyplot as pp
from matplotlib.patches import Ellipse, Wedge
import numpy as np
import pandas as pd
import datetime
from ovationpyme import ovation_utilities
from ovationpyme.ovation_prime import FluxEstimator, AverageEnergyEstimator
from ovationpyme.visual_test_ovation_prime import draw_weighted_flux

fig_directory = r"C:\Users\aduja\Documents\Research"

## This is going to generate a plot for each hour

def generate_daily_flux_plots(start_datetime, end_datetime):
    duration = end_datetime - start_datetime
    total_hours = duration.total_seconds() / 3600 # converting seconds to hours

    #generate plots for each hour between start and end times
    for hour in range(int(total_hours)):
        current_dt = start_datetime + datetime.timedelta(hours=hour)
        f = draw_weighted_flux(current_dt, atype='diff', jtype='energy')
        #mlatgridN, mltgridN, fluxgridN
        filename = f"flux_output_{current_dt.strftime('%Y%m%d_%H%M%S')}.png"
        full_path = f"{fig_directory}\\{filename}"
        f.savefig(full_path)
        pp.close(f)

        #saving to csv
        #df = pd.DataFrame({
        #    'Magnetic Latitude': mlatgridN.flatten(),
        #    'Magnetic Local Time': mltgridN.flatten(),
        #    'Flux': fluxgridN.flatten()
        #})
        #csv_filename: f"flux_data_{current_dt.strftime('%Y%m%d_%H%M%S')}.csv"
        #csv_filename = "flux_data_{}.csv".format(current_dt.strftime('%Y%m%d_%H%M%S'))
        #csv_full_path = f"{fig_directory}\\{csv_filename}"
        #df.to_csv(csv_full_path, index=False)


# start and end times for both observations
dt1s = datetime.datetime(2021, 12, 21, 19, 37, 23)
dt1e = datetime.datetime(2021, 12, 22, 3, 50, 23) # need to add an extra hour (really ends at 3:50:23)
dt2s = datetime.datetime(2021, 12, 2, 4, 45, 4)
dt2e = datetime.datetime(2021, 12, 3, 2, 34, 24) # need to add an extra hour (really ends at 2:34:24)

generate_daily_flux_plots(dt1s,dt1e)
generate_daily_flux_plots(dt2s,dt2e)



## Now generate one plot for entire time for each obs

# def generate_flux_plots_period(start_datetime, end_datetime):
#     duration = end_datetime - start_datetime
#     total_hours = int(duration.total_seconds() /3600) #convert sec to hours

#     all_mlats, all_mlts, all_fluxes = [],[],[]

#     #compute flux for each hour and collect data
#     for hour in range(total_hours + 1): #use +1 to include the final hour
#         current_dt = start_datetime + datetime.timedelta(hours=hour)
#         fig, flux, mlat, mlt = draw_weighted_flux(current_dt,atype='diff', jtype='energy')
#         all_fluxes.append(flux)
#         all_mlats.append(mlat)
#         all_mlts.append(mlt)

#     # #convert to numpy array
#     all_mlats = np.array(all_mlats)
#     all_fluxes = np.array(all_fluxes)
#     all_mlts = np.array(all_mlts)

#     mean_flux = np.mean(all_fluxes, axis=0) #average fluxes over time

#     fig, ax = plt.subplots()
#     c = ax.pcolormesh(all_mlts[0], all_mlats[0], mean_flux, shading='auto') 
#     fig.colorbar(c, ax=ax)
#     ax.set_title(f"Auroral Flux from {start_datetime} to {end_datetime}")
#     ax.set_xlabel('Magnetic Local Time')
#     ax.set_ylabel('Magnetic Latitude')
#     return fig

# fig1 = generate_flux_plots_period(dt1s, dt1e)
# filename = f"auroral_flux_20211221.png"
# full_path = os.path.join(fig_directory,filename)
# fig1.savefig(full_path)
# plt.close(fig1)

# fig2 = generate_flux_plots_period(dt2s, dt2e)
# filename = f"auroral_flux_20211202.png"
# full_path = os.path.join(fig_directory,filename)
# fig2.savefig(full_path)
# plt.close(fig2)


## Generating just one plot at the starting time 19:37:23 for obs 4202470103

# fig = draw_weighted_flux(dt1s)
# filename = f"4202470103_flux_start{dt1s.strftime('%Y%m%d_%H%M%S')}.png"
# full_path = os.path.join(fig_directory, filename)
# fig.savefig(full_path)
# plt.show()


    