## Trying to combine to one auroral flux plot for each obs

#polar plot to represent auroral flux that mimics the view of the Earth's pole

import matplotlib.pyplot as pp
import numpy as np
import datetime
from ovationpyme.ovation_prime import FluxEstimator
from geospacepy import satplottools, special_datetime


def plot_auroral_flux(dt1s, dt1e, atype='diff', jtype='energy'):
    flux_estimator = FluxEstimator(atype, jtype)
    
    time_range = [dt1s + datetime.timedelta(minutes=15 * i) for i in range(int((dt1e - dt1s).total_seconds() / 900))]
    ##^^ This will list times from dt1s to dt1e at 15-min intervals

    mlat_d, mlt_d, flux_d = [],[],[]

    for current_time in time_range:
        mlatgridN, mlatgridN, fluxgridN = flux_estimator.get_flux_for_time(current_time, hemi='N')
        mlat_d.append(mlatgridN)
        mlt_d.append(mlatgridN)
        flux_d.append(fluxgridN)

    mlat_d = np.array(mlat_d)
    mlt_d = np.array(mlt_d)
    flux_d = np.array(flux_d)

    mlat, mlt, flux = mlat_d.mean(axis=0), mlt_d.mean(axis=0), flux_d.mean(axis=0)

    #Convert MLT and MLAT to cartesian coords for polar plot
    X, Y = satplottools.latlt2cart(mlat.flatten(), mlt.flatten(), 'N')
    X = X.reshape(mlat.shape)
    Y = Y.reshape(mlt.shape)

    fig, ax = pp.subplots(subplot_kw={'projection': 'polar'}, figsize=(8, 8))

    # Convert MLT to radians for polar plot
    theta = np.radians(mlt * 15 - 90)  # MLT to degrees and shift so that MLT=0 is at top
    r = 90 - mlat  # Radial coordinate is 90 - MLAT to place the pole at the center
    Theta, R = np.meshgrid(theta, r)

    c = ax.pcolormesh(Theta, R, flux, shading='auto', cmap='viridis')
    pp.colorbar(c, label='Auroral Energy Flux')

    ax.set_title("Northern Hemisphere Auroral Flux")
    pp.show()

    # pp.figure(figsize=(12,6))
    # pp.pcolormesh(mlt_d[0], mlat_d[0], flux_d.mean(axis=0), shading='auto', cmap='viridis')
    # pp.colorbar(label='Auroral Energy Flux')
    # pp.xlabel('Magnetic Local Time (MLT)')
    # pp.ylabel('Magnetic Latitude (MLAT)')
    # pp.title('Northern Hemisphere Auroral Flux')
    # pp.show()

dt1s = datetime.datetime(2021, 12, 21, 19, 37, 23)
dt1e = datetime.datetime(2021, 12, 22, 3, 50, 23)


plot_auroral_flux(dt1s,dt1e)