import pyFAI, fabio
import os
import numpy as np
import glob
from info import *

runstart, runend = 0, 300 
for run in np.arange(runstart, runend+1):
    if run == 121 or not os.path.isdir(data_dir+"run%d"%run):
        print "No data for run %d!"%run
        continue 
    elif len(glob.glob(data_dir+"converted/run%d_evt*_data.npy"%run)) > 0:
        print "Run %d already converted! Skipping..."%run
        continue

    ## find the last event
    i_evt = 1
    while os.path.isfile(
            data_dir+"run%d/run_%d_evt_%d_%s.tif"%(run,run,i_evt,pads[0])):
        i_evt += 1
    i_evt -= 1

    I = []
    print "Loading data for run %d..." % run
    for i_pad, pad in enumerate(pads):
        ai = pyFAI.load("calibration/pad%d.poni" % i_pad)
        img = fabio.open(
          data_dir+"run%d/run_%d_evt_%d_%s.tif"%(run,run,i_evt,pad)).data
        if i_pad in [0,1]:
            img = img.T
        elif i_pad in [2]:
            img = np.flipud(img)
        elif i_pad in [3,4]:
            img = np.fliplr(img)
        I_i, tth, phi = ai.integrate2d(img, npt_rad=npt_rad, npt_azim=npt_azim, \
                   radial_range=(20., 110.), azimuth_range=(-180., 180.),\
                   unit='2th_deg')
        if I == []:
            I = I_i
        else:
            I += I_i
    I = np.fliplr(I.T)
    np.save(data_dir+"converted/run%d_evt%d_data.npy"%(run, i_evt), I)
    print "data saved"

phi_grid, tth_grid = np.meshgrid(phi, tth)
np.savez(data_dir+"converted/coordinates.npz", \
         phi=phi, tth=tth, phi_grid=phi_grid, tth_grid=tth_grid)
print "coordinates saved"
