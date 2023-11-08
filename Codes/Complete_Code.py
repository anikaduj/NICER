## This code contains all 3 codes I've been working on
## Plots spectrum (JP-NB), looks through cusp (NICER_POINTING), looks through AO (THROUGH_AO_WITH_B)


## Things to fix
# - for the AO code... need to find out how to combine all the data from the cusp code into a single file straight from this code
# - install heasoftpy and xspec

# - matplotlib notebook and inline for the jp-nb code (I commented these out... we'll see if it changes the output at all)


## Starting with the JP-NB code to plot the spectrum of the observation


print('------------OBSID SPECTRA STARTING------------')


from matplotlib.pyplot import *
import glob
import os
import scipy as sp
from scipy import stats
from scipy import optimize
from scipy.stats import norm
from scipy.optimize import curve_fit
import math

import numpy as np
import matplotlib.pyplot as plt

import heasoftpy as hsp ## Need to figure out how to install heasoftpy
## heasoftpy is installed in my Linux... but not in this pyspedas part
import xspec ## Need to figure out how to install xspec... should be easy once heasoftpy is installed

from astropy.table import Table
from astropy.time import Time
from astropy.io import fits
import astropy.units as u
from astropy import modeling
from astropy.modeling import fitting

## Set up NICER obsid directory
cwd = os.getcwd()

nicerdatadir = os.path.join(os.environ['HOME'])
nicerobsID = '5671040102' ## change this for different directories / also update the 'outdir' a few lines below
obsdir = os.path.join(nicerdatadir, nicerobsID)

# place cleaned output in a separate directory
outdir = os.path.join(os.environ['HOME'],'5671040102','nicer_output/'+nicerobsID+'_out')
# if outdir doesn't exit, create it
if not os.path.exists(outdir):
    os.makedirs(outdir)
    print(f'Created {outdir}')
    # copy the mkf file from the input directory to the outdir

## Setting mkf file to the new output directory
mkf = os.path.join(obsdir,'auxil',f'ni{nicerobsID}.mkf')
if os.path.exists(mkf):
    # see if mkf is gzipped
    cmd = f'cp {mkf} {outdir}/.' 
    stat = os.system(cmd)
    mkf = os.path.join(outdir, os.path.split(mkf)[1])
    print(f'Setting mkf file to {mkf}')
elif os.path.exists(mkf+'.gz'):
    # try to copy gzipped mkf
    cmd = f'cp {mkf}.gz {outdir}/.'
    print(cmd)
    os.system(cmd)
    mkf = os.path.join(outdir, os.path.split(mkf)[1])
    print(f'Setting mkf file to {mkf}')
    # cmd = f'gunzip -f {mkf}.gz' #* Unzip gz file
    # print(cmd)
    # stat = os.system(cmd)
cmd = f'chmod u+w {mkf}*' #* Changes the mode of path to the passed numeric mode
print(cmd)
stat = os.system(cmd)

## Good place to test to see the lengths
mkf_test = fits.open('/home/adujakovich/5671040102/auxil/ni5671040102.mkf')
mkf_data = mkf_test[1].data
mkf_hdr = mkf_test[1].header
#mkf_data.columns
#print(mkf_data.columns)
#print(mkf_test[1].data['TIME'][18336])
#print(mkf_test[1].data['SAT_LAT'][0]) 
#print(mkf_test[1].data['SAT_LON'][0]) 
#print(mkf_test[1].data['POSITION'][0]) 
#print(mkf_test[1].data['PNTUNIT'][0]) 
#print(mkf_test[1].data['POSITION']) 
#len(mkf_test[1].data['TIME'])


#%matplotlib inline
#plot(mkf_test[1].data['TIME'],mkf_test[1].data['TIME'])

## Create the nicerl2 task
tstart = Time.now()
print(f'Start at: {tstart.iso[:19]}')
nicerl2 = hsp.HSPTask('nicerl2')

nicerl2.clobber = "yes"
# nicerl2.cldir = outdir
# nicerl2.mkffile = mkf
# add the KP values to the mkf file during nicerl2 processing
nicerl2.geomag_path = "https://heasarc.gsfc.nasa.gov/FTP/caldb/data/gen/pcf/geomag/"
nicerl2.geomag_columns = "kp_noaa.fits(KP)"

resl2 = nicerl2(indir=nicerobsID, noprompt=True, cldir = outdir, mkffile = mkf)

tend = Time.now()
print(f'End at: {tend.iso[:19]}')
print(f'nicerl2 took: {(tend.mjd-tstart.mjd)*86400} seconds')


if resl2.returncode !=0:
    print('\n')
    for o in resl2.output[:]:
        print(o)

## Extract products from cleaned event file
## Can filter the times here with 'timefile=times' and updating the 'times' .txt file
times = '/home/adujakovich/Sept6NotAOTimes.txt'
clevt = f'{outdir}/ni{nicerobsID}_0mpu7_cl.evt' 
phafile = f'{outdir}/ni{nicerobsID}_0mpu7_cl.pha'
lcfile = f'{outdir}/ni{nicerobsID}_0mpu7_cl.lc' 

res = hsp.extractor(filename=clevt, phafile=phafile, clobber='yes', binlc=10.0,fitsbinlc=lcfile, 
                    eventsout='NONE', imgfile='NONE', regionfile='NONE', timefile='NONE', tcol='TIME',
                   ecol='PI', xcolf='RAWX', xcolh='RAWX',ycolf='RAWY', ycolh='RAWY',stokes='NONE')

print(res)

## Analyzing the NICER Spectra
# get the on-axis rmf
res = hsp.quzcif(mission='nicer', instrument='xti',detector='-',
             filter='-', date='-', time='-',expr='-',codename='MATRIX')
rmf = [x.split()[0] for x in res.output if 'nixtiref'  in x][0]

# get the on-axis arf
res = hsp.quzcif(mission='nicer', instrument='xti',detector='-',
             filter='-', date='-', time='-',expr='-',codename='SPECRESP')
arf = [x.split()[0] for x in res.output if 'nixtiaveonaxis'  in x][0]

xspec.AllData.clear()
spec = xspec.Spectrum(phafile)
spec.response = rmf
spec.response.arf = arf
spec.ignore('0.0-0.3, 10.0-**') # I can ignore and notice specific channels here

target = fits.open(spec.fileName)[1].header['OBJECT']
target

#exposure_time = spec.exposure
#print(f"Exposure Time: {exposure_time} seconds")

model = xspec.Model('wabs*agauss') # Change bknpow to ga for a guassina model
xspec.Fit.perform()
xspec.Fit.show()


## -------------------Plotting the Spectra ----------------------------------

#* wabs is a model used for xspec - photoelectric absorption
#%matplotlib ipympl # Need to find out the matplotlib notebook and inline for python
#%matplotlib inline
# ^^^^^ taking out these two because it's not supported in vscode... we'll see if it changes the plot at all
xspec.Plot.device='/null'
xspec.Plot.xAxis='keV' 
xspec.Plot('lda')
cr=xspec.Plot.y() # cr is the y values of the xspec
crerr = xspec.Plot.yErr() # crerr is the y value error arrays
en = xspec.Plot.x() # en are the log values of the energy in keV
enwid = xspec.Plot.xErr() # enwild is the x value error arrays
mop = xspec.Plot.model() # mop are the model values from the wabs model conducted in the previous cell
target = fits.open(spec.fileName)[1].header['OBJECT'] # just labeling the target to assign it to the title

#\n Time Filter: \n All Cusp Times +300s

fig = figure(figsize=[8,6])
ylabel('Cts/s/keV', fontsize=12)
xlabel('Energy (keV)', fontsize=12)
title('Target = '+target+' OBSID = '+nicerobsID+' \n Geomag Storm on 2022-09-06 \n Time Filter: None \n Exposure Time: 14950s', fontsize=12) 
yscale('log') 
xscale('log')
xlim([0.35,2.5])
ylim([0.09,8])
minorticks_off()
xticks([0.3, 0.4, 0.5, 0.6, 0.7, 1])
#xtickslabels([10e-3, 10e-5, 10e-6, 10e1])
yticks([0.1, 1, 1.6655534278604096])
plt.text(0.57500002,1.80, "Peak Amplitude (0.565)", color='blue', fontsize=10)
#set_xticklabels([0.3, 0.4, 0.5, 0.6, 1])
errorbar(en, cr, xerr=enwid, yerr=crerr, fmt='k.', alpha=0.2)
#plot(en, mop,'r')


#z = [0.56500002]
#q = [1.605352701552202]

#plot(z,q, marker="o", markersize = 4, markeredgecolor = "black", markerfacecolor = "red")

## Going through data to plot marker for amplitude and ion detection
xdata = np.asarray(en[10:])
ydata = np.asarray(cr[10:])
ydata1 = np.asarray(cr)
xdata1 = np.asarray(en)

ymax = max(ydata)
#print(max(ydata))
xpos = np.where(ydata == ymax)
xmax = xdata[xpos]
#print(xmax)


#for i in range(len(ydata)):
#    if ydata[i] ==1.9570098887972835:
 #       print(i)

# i = 24

#print(ydata)
#print(xdata)
#print(ydata_orig[24])
#print(xdata[14])

# amp is at (0.54500002,1.734622401433956)

print('------------JP-NB SPECTRA COMPLETED------------')







## Now take the observation and check if it looks through the cusp (NICER_POINTING)


print('------------CUSP POINTING STARTING------------')


from astropy.io import fits
from datetime import datetime
from datetime import timedelta
import geopack
import geopack.geopack as gp
import numpy as np
import pytz
import copy
from dateutil import parser
import math
import pandas as pd
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate

##### Import of data
# Data location is just a folder on my computer
# RootPath = 'C:\Users\aduja\NICER\home\adujakovich\Orbits\'
# DataName = 'ni4202470103'

RootPath = 'C:/Users/aduja/NICER/home/adujakovich/Orbits/'
DataName = 'ni4202470103'
DataCuspSE = np.array([]).reshape(0,4) #Initialize variable to list cusp observation times
DataNonCusp = np.array([]).reshape(0,4)
#print(DataNonCusp)
#print(DataCuspSE)
# fits.open is used to read the fits files
orb_fits = fits.open(RootPath+DataName+'.orb')
att_fits = fits.open(RootPath+DataName+'.att')
# the header contains relevant variabl information
orb_header = orb_fits[1].header
att_header = att_fits[1].header
# the data file contains the relevant data
orb_data = orb_fits[1].data
att_data = att_fits[1].data
# If you are searching for a variable, this is one way you can find it
#cols = att_data.columns
#print(cols)
# If you are searching for info on the header, this is how you print it
#print(repr(orb_header))
# apple
# Getting the orbit position data
Times_Pos = orb_data['TIME']
#print(Times_Pos[0])
X = orb_data['X']
Y = orb_data['Y']
Z = orb_data['Z']
Times_Q = att_data['TIME']
#print(Times_Q[0])

Quats = att_data['QPARAM']
#print(Quats)
# The TIME data is set as mission elapsed time since the start
# This needs to be calibrated to UTC time
NICERStart = datetime(year=2014,month=1,day=1,hour=0,minute=0,second=0)
#print(NICERStart)
NICERStart = datetime.timestamp(NICERStart)
#print(NICERStart)

LeapSec = 2
Times_Pos = Times_Pos + NICERStart - LeapSec
#print(Times_Pos[0])
Times_Q = Times_Q + NICERStart - LeapSec
#print(Times_Q[0])


##### Data manipulation
# Unfortunatley the quaternion data is higher resolution than the position data
# This means that we need to interpolate one into the other to match data points
Q1_Func = interpolate.interp1d(Times_Q, Quats[:,0], kind='cubic')
Q2_Func = interpolate.interp1d(Times_Q, Quats[:,1], kind='cubic')
Q3_Func = interpolate.interp1d(Times_Q, Quats[:,2], kind='cubic')
Q4_Func = interpolate.interp1d(Times_Q, Quats[:,3], kind='cubic')

# Apply the times of the position data to the quaternion function to generate
# quaternions at every time of the position data
Q1_at_P = Q1_Func(Times_Pos)
Q2_at_P = Q2_Func(Times_Pos)
Q3_at_P = Q3_Func(Times_Pos)
Q4_at_P = Q4_Func(Times_Pos)

# Initialize some new variables for the calculations
RX = np.zeros((len(Times_Pos ),3))
RY = np.zeros((len(Times_Pos ),3))
RZ = np.zeros((len(Times_Pos ),3))
RA = np.zeros((len(Times_Pos ),1))
Dec = np.zeros((len(Times_Pos ),1))
Roll = np.zeros((len(Times_Pos ),1))

## Added Aug 4 2023 ##
# Initialize variables for orbit data
OrbitNumber = 1
OrbitIter = None
OrbitData = {}

# Initialize CurrTime and PrevTime
CurrTime = Times_Pos[0]
PrevTime = CurrTime
OrbitIterData = []
## ##

# Now convert the quaternions into angles
for i in range(0,len(Times_Pos)):
    q1 = Q1_at_P[i] #x
    q2 = Q2_at_P[i] #y
    q3 = Q3_at_P[i] #z
    q4 = Q4_at_P[i] #w This is the scalar

    # This is the calculation of the rotation matrix. This rotation matrix rotates
    # The coordinate system (like ECI) into the coordinates pointing of the subject.
    # Therefore, multiplying the rotation matrix by [[1],[0],[0]] gives the pointing
    # of the X axis, multiply by [[0],[1],[0]] gives pointing of Y axis.... etc
    # because each row of the rotation matrix is the unit vector by which the
    # the axis points
    RX[i,0] = +q1*q1 - q2*q2 - q3*q3 + q4*q4
    RX[i,1] = 2*(q1*q2 + q3*q4)
    RX[i,2] = 2*(q1*q3 - q2*q4)
    RY[i,0] = 2*(q1*q2 - q3*q4)
    RY[i,1] = -q1*q1 + q2*q2 - q3*q3 + q4*q4
    RY[i,2] = 2*(q2*q3 + q1*q4)
    RZ[i,0] = 2*(q1*q3 + q2*q4)
    RZ[i,1] = 2*(q2*q3 - q1*q4)
    RZ[i,2] = -q1*q1 - q2*q2 + q3*q3 + q4*q4
    # Calculation of the pointing angles from quaternions
    # This comes from code by Fred Eckert and Aspire - specifically Quaternion.cpp
    #####
    # Aspire - High-performance, lightweight IP-based messaging middleware.
    # *   Copyright (C) 2010-2011 Ken Center, Fred Eckert & Rod Green
    #####
    #THIS IS ROLL
    denom = q1*q3 - q2*q4 #Prep denominator
    psi = math.atan2(-(q2*q3 + q1*q4),denom) #Calculate arc tan with quadrant dependency
    psi = psi + 0
    # These if statements bring the roll into the positive rotation by wrapping 2pi
    if (psi < -np.pi):
        n = (np.ceil(-psi / (2*np.pi)))
        Roll[i] = psi+(n*2*np.pi)
    elif (psi >= np.pi):
        n = (psi / (2*np.pi));
        Roll[i] = psi-(n*2*np.pi)
    else:
        Roll[i] = psi
    #####
    # THIS IS RIGHT ASCENSION
    denom = q3*q1 + q2*q4 #Prep denominator
    phi = math.atan2(q3*q2 - q1*q4,denom) #Calculate arc tan with quadrant dependency
    # These if statements bring the RA into the positive rotation by wrapping 1pi
    if (phi < 0):
        n = (np.ceil(-phi / (2*np.pi)))
        RA[i] = phi+(n*2*np.pi)
    elif (phi >= 2*np.pi):
        n = (phi / (2*np.pi));
        RA[i] = phi-(n*2*np.pi)
    else:
        RA[i] = phi
    #####
    # THIS IS DECLINATION
    acos = q3*q3 + q4*q4 - q2*q2 - q1*q1 #This is actually an element of a directional cosine matrix
    if (acos > 1):
        acos = 1
    if (acos < -1):
        acos = -1
    theta = math.acos(acos) #get the angle
    # These if statements bring the Declination to be defined from the equator
    if (theta >=0):
        Dec[i] = np.pi/2 - theta
    else:
        Dec[i] = -np.pi/2 - theta

# Convert the RA and Dec to a point in ECI
x_eci = np.cos(RA)*np.cos(Dec)
y_eci = np.sin(RA)*np.cos(Dec)
z_eci = np.sin(Dec)
# Summarize it in a column
RaDec = np.c_[x_eci,y_eci,z_eci]
# Initialize some variable matricies
P_gsm = np.zeros((len(Times_Pos),3))
RaDec_gsm = np.zeros((len(RaDec[:,0] ),3))

## Added Aug 4 2023 ##
# Set the reference time for the next pass through the loop
PrevTime = CurrTime
CurrTime = Times_Pos[i] # Update CurrTime with the current time
## ##


###### Conversion to GSM and splitting into orbits
# Iterate through the matricies to convert to GSM with geopack.
PrevTime = Times_Pos[0] #Initiate the reference time
OrbitData = {} #Initiate the dictionary
OrbitNumber = 1 #Call out the first orbit
OrbitIter = [] #Make a blank matrix

#print(len(Times_Pos))

for i in range(len(Times_Pos)):
    gp.recalc(Times_Pos[i]) # set the time in geopack to get proper location of ECI
    tmp1,tmp2,tmp3 = gp.geigeo(X[i],Y[i],Z[i], 1) # the 1 means go from GEI to GEO
    # P_gsm[i,0],P_gsm[i,1],P_gsm[i,2] = gp.geogsm(tmp1,tmp2,tmp3, 1) # from GEO to GSM
    P_gsmX,P_gsmY,P_gsmZ = gp.geogsm(tmp1,tmp2,tmp3, 1) # from GEO to GSM
    # P_gse[i,0],P_gse[i,1],P_gse[i,2] = gp.gsmgse(P_gsm[i,0],P_gsm[i,1],P_gsm[i,2], 1)

    tmp1,tmp2,tmp3 = gp.geigeo(RaDec[i,0],RaDec[i,1],RaDec[i,2], 1)
    # RaDec_gsm[i,0],RaDec_gsm[i,1],RaDec_gsm[i,2] = gp.geogsm(tmp1,tmp2,tmp3, 1)
    RaDec_gsmX,RaDec_gsmY,RaDec_gsmZ = gp.geogsm(tmp1,tmp2,tmp3, 1)
    # RaDec_gse[i,0],RaDec_gse[i,1],RaDec_gse[i,2] = gp.gsmgse(RaDec_gsm[i,0],RaDec_gsm[i,1],RaDec_gsm[i,2], 1)

    CurrTime = Times_Pos[i] #Set the current time
    if CurrTime-PrevTime<=10: #Compare the current time to the previous time
        # Orbit data is in 10 second time steps,
        # so a time jump of less than 10 seconds
        # means these are sequential data points.
        # Append this data to the matrix of this orbit
        appenddata = [CurrTime,P_gsmX/1000,P_gsmY/1000,P_gsmZ/1000,RaDec_gsmX,RaDec_gsmY,RaDec_gsmZ]
        OrbitIterData.append(appenddata)
        #if OrbitIter is None:
        #    OrbitIter = np.array[(appenddata)]
        #else:
        #    OrbitIter = np.vstack((OrbitIter, appenddata))

    if CurrTime-PrevTime>600 or i == len(Times_Pos)-1:
        # If the time step is greater than 10 seconds...
        # here say its greater than 10 minutes, then its
        # a new orbit. Or, if its the last entry.
        # Convert to an array
        #OrbitIter = np.asarray(OrbitIter) # Commented this out on aug 4 2023
        # Define its dictionary name
        var_name = "OrbitNumber%d" % OrbitNumber
        # Put it in the dictionary
        OrbitData[var_name] = np.array(OrbitIterData)
        # Reset the iterative matrix
        #OrbitIter = None # Used to be [] but changed to None on aug 4 2023
        OrbitIterData = []
        # Increment the orbit number
        OrbitNumber = OrbitNumber+1
        # Append the data to the new blank matrix
        appenddata = [CurrTime,P_gsmX/1000,P_gsmY/1000,P_gsmZ/1000,RaDec_gsmX,RaDec_gsmY,RaDec_gsmZ]
        OrbitIter.append(appenddata)
    # Set the reference time for the next pass through the loop
    PrevTime = CurrTime


##### Plotting
##### Iterate over the number of orbits
#for r in range(1,len(OrbitData)+1):
##### Iterate over just one orbit (for testing)

for r in range(1,2):

    # Grab just the one orbit of data
    TheDataToPlot = OrbitData["OrbitNumber"+str(r)]
    np.savetxt("/Users/aduja/OneDrive/Documents/School - Research/"+"Orbit"+str(r)+"_DataValues.csv",TheDataToPlot,delimiter=',',newline='\n', header='time_UTC,ISS X_GSM,ISS Y_GSM,ISS Z_GSM,NICER X_Point,NICER Y_Point,NICER Z_Point', fmt='%.10f')    
    # ^^^ This has the correct orbit values...
    # TheDataToPlot is Time, ISS X, ISS Y, ISS Z, NICER X, NICER Y, NICER Z
    #print(TheDataToPlot)
    #print(np.sqrt((TheDataToPlot[1,1]-TheDataToPlot[2,1])**2+(TheDataToPlot[1,2]-TheDataToPlot[2,2])**2+(TheDataToPlot[1,3]-TheDataToPlot[2,3])**2))
    # apple
    #Step through the orbit positions to find the times NICER views the cusp
    scale = 10# [km] the size steps we project the view
    limit = 1000# [steps] maximum number of steps to iterate through
    # cuspgroundradius = 500# [km] the radius of the cusp on earth surface
    DMSPalt = 6371 #km #850+6371
    cuspMLT = 12
    cuspLAT = 75
    cuspMLTwidth = 1.5 #Hours of half length
    cuspLATwidth = 4 #degrees of half width
    x_cusp_gsm = DMSPalt*np.cos(cuspLAT*np.pi/180)
    y_cusp_gsm = 0
    z_cusp_gsm = DMSPalt*np.sin(cuspLAT*np.pi/180)
    cusptime = [] # initiate the flag matrix
    sw_params = {}
    sw_params["p_dyn"] = 2
    sw_params["by_gsm"] = -3
    sw_params["bz_gsm"] = -2
    param = [sw_params["p_dyn"], -30, sw_params["by_gsm"], sw_params["bz_gsm"], 0, 0, 0, 0, 0, 0]
    rho = 1.0+850/6371 # altitude in re for cusp position #1.0+850/6371
    rho = 1
    # rho = 1.2
    # rho = 1.4
    # rho = 1.6
    # Dont change these ======================
    # Constants from Tsyganenko and Russell JGR 1999
    phi_c0 = 0.24
    alpha1 = 0.1287
    alpha2 = 0.0314
    # iterate through the length of the orbit data
    for j in range(0,len(TheDataToPlot[:,0])):
        # print('next')
        dView = 1 #initial step
        # grab the appropriate data to iterate with
        Time = TheDataToPlot[j,0]
        X = TheDataToPlot[j,1]
        Y = TheDataToPlot[j,2]
        Z = TheDataToPlot[j,3]
        ps = gp.recalc(Time) #ps is the dipole tilt in radians
        # phi_1 = phi_c0 - (alpha1*ps+alpha2*ps**2)
        # x1 = np.sqrt(rho)
        # x2 = np.sqrt(rho+np.sin(phi_1)**(-2)-1)
        # phi_c = np.arcsin(x1/x2)+ps
        # z_cusp_gsm = rho*np.sin(np.pi/2-phi_c)*6371
        # x_cusp_gsm = rho*np.cos(np.pi/2-phi_c)*6371
        # y_cusp_gsm = 0
        # cusplatitude = np.arctan(abs(z_cusp_gsm)/abs(x_cusp_gsm))#*180/np.pi
        #print(x_cusp_gsm)
        #print(z_cusp_gsm)
        # print(cusplatitude)
        # apple

        #print(z_cusp_gsm)

        # while loop steps through the observation direction
        while True:
            step = scale*dView #How far to step
            # the pointing is a unit vector (magnitude 1)
            # so extending it by a uniform value projects
            # the XYZ position along the vector
            X = (X+TheDataToPlot[j,4]*step)
            Y = (Y+TheDataToPlot[j,5]*step)
            Z = (Z+TheDataToPlot[j,6]*step)
            posMLT = (np.arctan(Y/X)*180/np.pi)*(1/15)+12
            posLAT = np.arctan(Z/X)*180/np.pi
            # rdist = np.sqrt(X**2 + Y**2 + Z**2)
            # if rdist > DMSPalt:
            #     Xcusp, Ycusp, Zcusp, _, _, _  = gp.trace(X, Y, Z, dir=1, rlim=10, r0=0.9999, parmod=param,
            #                               exname='t96', inname='igrf', maxloop=10000)

            # trace the posiiton to the surface of the earth
            # print('starting trace')
            # print(rho)
            # print(rdist)
            # print(Xsurf)
            # print(Ysurf)
            # print(Zsurf)
            # print(np.sqrt(Xsurf**2 + Ysurf**2 + Zsurf**2))

            # apple
            # print('traced')
            # Determine if the position is inside the cusp region
            Distance2CuspMLT = abs(posMLT-cuspMLT)
            Distance2CuspLAT = abs(posLAT-cuspLAT)

            # MLTang = np.arctan(Ysurf/Xsurf)*180/np.pi
            # print(MLTang)

            if Distance2CuspMLT <= cuspMLTwidth:
                # print('si')
                if Distance2CuspLAT<=cuspLATwidth:
                    print('yes')
                    cusptime = np.append(cusptime,1)
                    break

            # print(np.sqrt(Distance2CuspX**2+Distance2CuspY**2+Distance2CuspZ*2))
            # if np.sqrt(Distance2CuspX**2+Distance2CuspY**2+Distance2CuspZ*2)< cuspgroundradius:
                # print('yes')
                # cusptime = np.append(cusptime,1)
                # break
            # should the position never be inside the cusp region
            # at the end of the distance, then say that it never
            # looked through the cusp, beak the while loop.
            if step >=scale*limit:
                #print('no')
                cusptime = np.append(cusptime,0)
                break
            dView = dView + 1

    # Calculate the proper view angle
    # Start = first ISS position of observation
    Start = np.asarray([TheDataToPlot[0,1],TheDataToPlot[0,2],TheDataToPlot[0,3]])
    # Half = middle time of observation
    Half = int(len(TheDataToPlot[:,0])/2)
    # Mid = middle ISS position of observatino using Half
    Mid  = np.asarray([TheDataToPlot[Half,1],TheDataToPlot[Half,2],TheDataToPlot[Half,3]])
    # End = last ISS position of observation
    End = np.asarray([TheDataToPlot[-1,1],TheDataToPlot[-1,2],TheDataToPlot[-1,3]])
    First = Mid-Start #Vector from middle to start
    Second =End-Start #Vector from end to start
    Direction = np.cross(First,Second) #Cross product of these vectors

    Direction=Direction / np.sqrt(np.sum(Direction**2)) #normalize
    # Calculate the azimuth and elevation angles of the normal vector to the orbit plane
    a=np.degrees(np.arctan(Direction[1]/Direction[0]))
    e=np.degrees(np.arctan(np.sqrt((Direction[0]**2)+(Direction[1]**2))/Direction[2]))
    # print(a)
    # print(e)
    fig = plt.figure()
    ax =  plt.axes(projection = '3d')
    ax.set_title('NICER Pointing GSM')
    ax.set_xlabel('GSM X (km)')
    ax.set_ylabel('GSM Y (km)')
    ax.set_zlabel('GSM Z (km)')
    ax.set_xticks([-6000, -3000, 0, 3000, 6000])
    ax.set_yticks([-6000, -3000, 0, 3000, 6000])
    ax.set_zticks([-6000, -3000, 0, 3000, 6000])
    ax.set_xlim(-7000.0, 7000.0) #was originally 7000 for each limit
    ax.set_ylim(-7000.0, 7000.0)
    ax.set_zlim(-7000.0, 7000.0)
    ax.view_init(e-90, a) #Dictates what angle you want the plot shown first before you can rotate it
    #ax.plot([0,7000], [0,0], [0,0], color="y") # plot the Earth-sun lin
    ax.plot3D(TheDataToPlot[:,1], TheDataToPlot[:,2], TheDataToPlot[:,3], color = 'r',label='ISS orbit') #The +500 just makes sure that it's not too close to the Earth so you can actually see it
    #print(x_cusp_gsm, y_cusp_gsm, z_cusp_gsm)
    #ax.scatter([x_cusp_gsm], [y_cusp_gsm], [z_cusp_gsm+500], color = 'b',label='CUSP POSITION',s=100)
#Plot an Earth sphere (thank you StackOverflow)
    plotradius = 4000 #4000
    ud, vd = np.mgrid[-np.pi/2:np.pi/2:50j, 0:np.pi:30j]
    Sxd = np.cos(ud)*np.sin(vd)*plotradius
    Syd = np.sin(ud)*np.sin(vd)*plotradius
    Szd = np.cos(vd)*plotradius
    ax.plot_surface(Sxd, Syd, Szd, color="#fffff4",linewidth=0.0)
    un, vn = np.mgrid[np.pi/2:3*np.pi/2:50j, 0:np.pi:30j]
    Sxn = np.cos(un)*np.sin(vn)*plotradius
    Syn = np.sin(un)*np.sin(vn)*plotradius
    Szn = np.cos(vn)*plotradius
    ax.plot_surface(Sxn, Syn, Szn, color="k",linewidth=0.0)

#Plot the surface of the cusp region
    cuspu, cuspv = np.mgrid[-cuspMLTwidth*(15/1)*(np.pi/180):cuspMLTwidth*(15/1)*(np.pi/180):20j, (90-cuspLAT-cuspLATwidth)*(np.pi/180):(90-cuspLAT+cuspLATwidth)*(np.pi/180):10j]
    Sxd = np.cos(cuspu)*np.sin(cuspv)*DMSPalt
    Syd = np.sin(cuspu)*np.sin(cuspv)*DMSPalt
    Szd = np.cos(cuspv)*DMSPalt
    ax.plot_surface(Sxd, Syd, Szd, color='red',linewidth=0.0)



# Plot the first vector so you get the label
    ax.quiver(TheDataToPlot[0,1],TheDataToPlot[0,2],TheDataToPlot[0,3],TheDataToPlot[0,4],TheDataToPlot[0,5],TheDataToPlot[0,6],length=2000,color="k",label='NICER Pointing Vector')
    #ax.quiver(TheDataToPlot[0,1],TheDataToPlot[0,2],TheDataToPlot[0,3],TheDataToPlot[0,4],TheDataToPlot[0,5],TheDataToPlot[0,6],length=2000,color="g")
    #ax.quiver(TheDataToPlot[-1,1],TheDataToPlot[-1,2],TheDataToPlot[-1,3],TheDataToPlot[0,4],TheDataToPlot[0,5],TheDataToPlot[0,6],length=2000,color="r")



 #For loop to go over the pointing vectors
    for i in range(1,len(TheDataToPlot[:,0])):
        if i % 5 == 0: # use only every n'th value. It gets confusing with too many
            if cusptime[i] == 0:
                ax.quiver(TheDataToPlot[i,1],TheDataToPlot[i,2],TheDataToPlot[i,3],TheDataToPlot[i,4],TheDataToPlot[i,5],TheDataToPlot[i,6],length=3000,color="k")



                #first_noncusptime_index = np.where(cusptime == 0)[0][0]
                #second_noncusptime_index = np.where(cusptime == 0)[0][1]
                #third_noncusptime_index = np.where(cusptime == 0)[0][2]
                #fourth_noncusptime_index = np.where(cusptime == 0)[0][3]
                #last_noncusptime_index = np.where(cusptime == 0)[0][-1]


                #first_noncusptime_value = TheDataToPlot[first_noncusptime_index, 0]
                #second_noncusptime_value = TheDataToPlot[second_noncusptime_index, 0]
                #third_noncusptime_value = TheDataToPlot[third_noncusptime_index, 0]
                #fourth_noncusptime_value = TheDataToPlot[fourth_noncusptime_index, 0]
                #last_noncusptime_value = TheDataToPlot[last_noncusptime_index, 0]


                #print("First noncusptime value:", first_noncusptime_value)
                #print("Second noncusptime value:", second_noncusptime_value)
                #print("Third noncusptime value:", third_noncusptime_value)
                #print("Fourth noncusptime value:", fourth_noncusptime_value)
                #print("Last noncusptime value:", last_noncusptime_value)



            if cusptime[i] == 1:
                ax.quiver(TheDataToPlot[i,1],TheDataToPlot[i,2],TheDataToPlot[i,3],TheDataToPlot[i,4],TheDataToPlot[i,5],TheDataToPlot[i,6],length=3000,color="y")
                #ax.plot(TheDataToPlot[i,1],TheDataToPlot[i,2],TheDataToPlot[i,3],TheDataToPlot[i,4],TheDataToPlot[i,5],TheDataToPlot[i,6], 'bo', markersize=20)
                #^^ This just created a giant blue dot at the top of the north pole

                #print("Yellow line data:")
                #print("Time:", TheDataToPlot[i, 0])
                #print("ISS X:", TheDataToPlot[i, 1])
                #print("ISS Y:", TheDataToPlot[i, 2])
                #print("ISS Z:", TheDataToPlot[i, 3])
                #print("NICER X:", TheDataToPlot[i, 4])
                #print("NICER Y:", TheDataToPlot[i, 5])
                #print("NICER Z:", TheDataToPlot[i, 6])
                #print("---")



                #first_cusptime_index = np.where(cusptime == 1)[0][0]
                #second_cusptime_index = np.where(cusptime == 1)[0][1]
                #third_cusptime_index = np.where(cusptime == 1)[0][2]
                #fourth_cusptime_index = np.where(cusptime == 1)[0][3]
                #last_cusptime_index = np.where(cusptime == 1)[0][-1]



                #first_cusptime_value = TheDataToPlot[first_cusptime_index, 0]
                #second_cusptime_value = TheDataToPlot[second_cusptime_index, 0]
                #third_cusptime_value = TheDataToPlot[third_cusptime_index, 0]
                #fourth_cusptime_value = TheDataToPlot[fourth_cusptime_index, 0]
                #last_cusptime_value = TheDataToPlot[last_cusptime_index, 0]

                #first_noncusptime_index_aftercusp = np.where(cusptime == 0)[0][np.where(np.where(cusptime == 0)[0] > last_cusptime_index)[0][0]]
                #second_noncusptime_index_aftercusp = np.where(cusptime == 0)[0][np.where(np.where(cusptime == 0)[0] > first_noncusptime_index_aftercusp)[0][0]]
                #third_noncusptime_index_aftercusp = np.where(cusptime == 0)[0][np.where(np.where(cusptime == 0)[0] > second_noncusptime_index_aftercusp)[0][0]]
                #fourth_noncusptime_index_aftercusp = np.where(cusptime == 0)[0][np.where(np.where(cusptime == 0)[0] > third_noncusptime_index_aftercusp)[0][0]]
                #fifth_noncusptime_index_aftercusp = np.where(cusptime == 0)[0][np.where(np.where(cusptime == 0)[0] > fourth_noncusptime_index_aftercusp)[0][0]]

                #indexes_before_cusp = np.arange(first_cusptime_index - 10, first_cusptime_index)
                #indexes_after_cusp = np.arange(last_cusptime_index + 1, last_cusptime_index + 11)

                #values_before_cusp = TheDataToPlot[indexes_before_cusp, 0]
                #values_after_cusp = TheDataToPlot[indexes_after_cusp, 0]


                #print("Fourth cusptime value:", fourth_cusptime_value)
                #print("Second cusptime value:", second_cusptime_value)
                #print("Third cusptime value:", third_cusptime_value)
                #print("Fourth cusptime value:", fourth_cusptime_value)
                #print("Last cusptime value:", last_cusptime_value)

                #ax.plot(TheDataToPlot[last_cusptime_index, 1], TheDataToPlot[last_cusptime_index, 2], TheDataToPlot[last_cusptime_index, 3], 'go', markersize=5)
                #ax.plot(TheDataToPlot[first_cusptime_index, 1], TheDataToPlot[first_cusptime_index, 2], TheDataToPlot[first_cusptime_index, 3], 'bo', markersize=5)
                #ax.text(TheDataToPlot[last_cusptime_index, 1], TheDataToPlot[last_cusptime_index, 2], TheDataToPlot[last_cusptime_index, 3]-1300,
                #    'tstop', color='g', fontsize=8) #'t\u2080'
                #ax.text(TheDataToPlot[first_cusptime_index, 1], TheDataToPlot[first_cusptime_index, 2], TheDataToPlot[first_cusptime_index, 3]-1300,
                #    'tstart', color='b', fontsize=8) #'t\u0066'
                                
                #ax.plot(TheDataToPlot[last_noncusptime_index, 1], TheDataToPlot[last_noncusptime_index, 2], TheDataToPlot[last_noncusptime_index, 3], 'go', markersize=5)
                #ax.plot(TheDataToPlot[fifth_noncusptime_index_aftercusp, 1], TheDataToPlot[fifth_noncusptime_index_aftercusp, 2], TheDataToPlot[fifth_noncusptime_index_aftercusp, 3], 'bo', markersize=5)
                #ax.text(TheDataToPlot[last_noncusptime_index, 1], TheDataToPlot[last_noncusptime_index, 2], TheDataToPlot[last_noncusptime_index, 3]-1800,
                #    'tstop', color='g', fontsize=8) #'t\u2080'
                #ax.text(TheDataToPlot[fifth_noncusptime_index_aftercusp, 1], TheDataToPlot[fifth_noncusptime_index_aftercusp, 2], TheDataToPlot[fifth_noncusptime_index_aftercusp, 3]-1800,
                #    'tstart', color='b', fontsize=8) #'t\u0066'


    mask = np.where(np.diff(cusptime) != 0)
    #print(mask)
    # UTCcusptime = TheDataToPlot[:,0][mask]
    UTCcusptime = np.reshape(TheDataToPlot[:,0][mask], (-1, 2))
    #print(UTCcusptime)
    METcusptime = UTCcusptime - NICERStart + LeapSec
    UTCtime = TheDataToPlot[:,0]
    #print(UTCtime[0])
    METtime = UTCtime - NICERStart + LeapSec
    #print(METtime[0])
    CuspStartEnd = np.concatenate((UTCcusptime, METcusptime),axis=1)
    # print(CuspStartEnd)
    DataCuspSE = np.vstack([DataCuspSE,CuspStartEnd])
    #METUTC = np.concatenate((UTCtime,METtime), axis=1)
    METandUTC = np.vstack([UTCtime,METtime]).T



    #first_noncusptime_value_aftercusp = TheDataToPlot[first_noncusptime_index_aftercusp, 0]
    #fifth_noncusptime_value_aftercusp = TheDataToPlot[fifth_noncusptime_index_aftercusp, 0]
    #print("First noncusp time value after the cusp:", first_noncusptime_value_aftercusp)
    #first_noncusptime_value_aftercusp_MET = first_noncusptime_value_aftercusp - NICERStart + LeapSec
    #fifth_noncusptime_value_aftercusp_MET = fifth_noncusptime_value_aftercusp - NICERStart + LeapSec
    #print("First noncusp time value in MET after the cusp:", first_noncusptime_value_aftercusp_MET)
    #last_noncusptime_value_MET = last_noncusptime_value - NICERStart + LeapSec
    #print("Last noncusp time value in MET:", last_noncusptime_value_MET)

    #UTCnoncusptime = np.array([[fifth_noncusptime_value_aftercusp, last_noncusptime_value]])
    #METnoncusptime = np.array([[fifth_noncusptime_value_aftercusp_MET, last_noncusptime_value_MET]])
    #NonCuspStartEnd = np.concatenate((UTCnoncusptime, METnoncusptime),axis=1)
    #DataNonCusp = np.vstack([DataNonCusp,NonCuspStartEnd])


    ax.legend()
    plt.show()
    fig.savefig('/Users/aduja/OneDrive/Documents/School - Research/OrbitResults'+str(r)+'_GSM.png')
    plt.close()

    print('Saving data')
    #np.savetxt("/Users/aduja/OneDrive/Documents/School - Research/OrbitResults"+str(r)+"_Data.csv",OrbitData["OrbitNumber1"],delimiter=',',newline='\n',header='time,ISS X_GSM,ISS Y_GSM,ISS Z_GSM,NICER X_Point,NICER Y_Point,NICER Z_Point', fmt='%.10f')
    #np.savetxt("/Users/aduja/OneDrive/Documents/School - Research/"+"Orbit"+str(r)+"_DataValues.csv",TheDataToPlot,delimiter=',',newline='\n', header='time,ISS X_GSM,ISS Y_GSM,ISS Z_GSM,NICER X_Point,NICER Y_Point,NICER Z_Point', fmt='%.10f')    


    np.savetxt("/Users/aduja/OneDrive/Documents/School - Research/OrbitResults"+str(r)+"_CuspSE.csv",CuspStartEnd,delimiter=',',newline='\n',header='StartUTC,EndUTC,StartMET,EndMET', fmt='%.10f')
    np.savetxt("/Users/aduja/OneDrive/Documents/School - Research/OrbitResults"+str(r)+"_UTC_MET.csv",METandUTC,delimiter=',',newline='\n',header='Time-UTC,Time-MET', fmt='%.10f')



np.savetxt("/Users/aduja/OneDrive/Documents/School - Research/OrbitResults"+str(r)+"_CuspSE.csv",DataCuspSE,delimiter=',',newline='\n',header='StartUTC,EndUTC,StartMET,EndMET', fmt='%.10f')
np.savetxt("/Users/aduja/OneDrive/Documents/School - Research/OrbitResults"+str(r)+"_NonCusp.csv",DataNonCusp,delimiter=',',newline='\n',header='StartUTC,EndUTC,StartMET,EndMET', fmt='%.10f')

print("Doin Science Yo!")



print('------------CUSP POINTING COMPLETED------------')






## Finally using the AO code to check if the observation looks through the AO


print('------------AO POINTING STARTING------------')


import math
import numpy as np
from numpy import sqrt
import csv

re = 6371. # Earth radii in km

RootPath = 'C:/Users/aduja/OneDrive/Documents/School - Research/'
DataName = 'AllOrbitResults_Data'
filepath = RootPath+DataName+'.csv'
Nicer_Pointing = []
Nicer_Position = []

with open(filepath, 'r') as csvfile:
    csv_reader = csv.reader(csvfile)
    data = list(csv_reader)
    ## Save the Nicer Pointing and Nicer Position entries to an array
for row in data[1:]:
    Nicer_X_point = float(row[4])  
    Nicer_Y_point = float(row[5])  
    Nicer_Z_point = float(row[6])
    Nicer_X_pos = float(row[1])  
    Nicer_Y_pos = float(row[2])  
    Nicer_Z_pos = float(row[3])

    Nicer_Pointing.append([Nicer_X_point, Nicer_Y_point, Nicer_Z_point])
    Nicer_Position.append([Nicer_X_pos, Nicer_Y_pos, Nicer_Z_pos])
    

r_position = Nicer_Position # defining the position of NICER to the r variable
s_pointing = Nicer_Pointing # defining pointing of Nicer to s variable

#r = np.array([209.1986317224,-4046.1352032042,5450.5767400454])
#s = np.array([-0.1620193516,-0.9748193518,-0.1532219337])

r = np.array(r_position)
s = np.array(s_pointing)
#print(len(r))

r_re = r/re # position of NICER in units of RE
#print(len(r_re))
r_mag = np.linalg.norm(r_re, axis=1)

lat_low = 60. #deg min for auroral region
lat_high = 75. # deg max auroral region
lat_tmp = np.abs(np.degrees(np.arctan(r_re[:,2]/np.sqrt(r_re[:,0]**2+r_re[:,1]**2))))
#print('Latitude', lat_tmp)
#print(len(lat_tmp))

dr = 250./re # step size in km
rmax = 2.5 # RE max distance for ray trace
status = 0 # if not auroral viewing return 0.  If aurora -> 1

#print('Initial before loop',lat_tmp[0],r_mag[0],dr)
status_list = []

# step through view direction
while (r_mag < rmax).any() and (status == 0):
    # step view forward
    r_re = r_re + s * dr
    r_mag = np.linalg.norm(r_re, axis=1)
    lat_tmp = np.abs(np.degrees(np.arctan(r_re[:,2]/np.sqrt(r_re[:,0]**2+r_re[:,1]**2))))
    #print(lat_tmp,r_mag,r_re)

    LatMagPos = np.column_stack((lat_tmp, r_mag, r_re))
    boundaries = ((lat_tmp < lat_high) & (lat_tmp > lat_low)).astype(int) # .astype(int) means that True = 1 and False = 0
    status_list.append(boundaries)

    if boundaries.any(): # this is saying if those conditions are met for any of the values, then enter a 1
        status = 1

status_list = np.array(status_list)
status_listT = status_list.T
StatusLatMagPos = np.column_stack((status_listT, LatMagPos))
######### Old code in place of *boundaries = [through] status = 1   ################
#    if (lat_tmp < lat_high).any() and (lat_tmp > lat_low).any():
#        status = 1
#        print(status)
####################################################################################
#np.savetxt("/Users/aduja/OneDrive/Documents/School - Research/"+"Results"+"_LatMagPos.csv",LatMagPos,delimiter=',',newline='\n',header='lat_tmp,r_mag,r_re', comments='')
#np.savetxt("/Users/aduja/OneDrive/Documents/School - Research/"+"Results"+"_status.csv",status_listT,delimiter=',',newline='\n',header='status', comments='')
np.savetxt("/Users/aduja/OneDrive/Documents/School - Research/"+"Results"+"_AuroraView.csv",StatusLatMagPos,delimiter=',',newline='\n',header='Through AO?, Lat, R_Mag, NICER X Position, NICER Y Position, NICER Z Position', comments='')
np.savetxt("/Users/aduja/OneDrive/Documents/School - Research/"+"Results"+"_S++.csv",s,delimiter=',',newline='\n',header='NICER X Pointing, NICER Y Pointing, NICER Z Pointing', comments='')

    

print('------------AO POINTING COMPLETED------------')