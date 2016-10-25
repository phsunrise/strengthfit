import pyFAI, fabio
import matplotlib.pyplot as plt
import numpy as np
import os, sys
sys.path.insert(0, "/Users/phsun/LL20/information/")
import info

#run = int(sys.argv[1])
pads = ['Cspad.0', 'Cspad.1', 'Cspad.2', 'Cspad2x2.1']

## diamond lines (tth values) for 10 keV
diamond_lines = {
        '(111)': 35.0385, '(220)': 58.8882, '(311)': 70.3985,\
        '(400)': 88.0849, '(331)': 98.5003}

for run in range(174, 200):
    if run == 121 or not os.path.isdir(info.data_dir+"run%d"%run):
        continue
    I = None
    for i_pad, pad in enumerate(pads):
        ai = pyFAI.load("calibration/pad%d.poni" % i_pad)
        ## find the last event
        i_evt = 1
        while os.path.isfile(
                info.data_dir+"run%d/run_%d_evt_%d_%s.tif"%(run,run,i_evt,pad)):
            i_evt += 1
        i_evt -= 1
        img = fabio.open(
          info.data_dir+"run%d/run_%d_evt_%d_%s.tif"%(run,run,i_evt,pad)).data
        if i_pad in [0,1]:
            img = np.rot90(img, 3) 
        elif i_pad in [2]:
            img = np.rot90(img, 2) 
        elif i_pad in [3,4]:
            pass
        I_i, tth, chi = ai.integrate2d(img, npt_rad=500, npt_azim=500, \
                   radial_range=(20., 110.), azimuth_range=(-180., 180.),\
                   unit='2th_deg')
        if I == None:
            I = I_i
        else:
            I += I_i
        
    fig = plt.figure(figsize=(15, 20))
    ax1 = fig.add_subplot(2,1,1)
    I[I>0.3*np.max(I)] = 0.5*np.max(I)
    ax1.imshow(I.T, origin='lower', aspect='auto',\
               extent=[min(chi), max(chi), min(tth), max(tth)])
    ## draw unshocked diamond lines
    for line, line_tth in diamond_lines.iteritems():
        ax1.axhline(y=line_tth, ls='--', color='r')
        ax1.text(x=185, y=line_tth, s=line, color='r')

    ax1.set_xlabel("chi (deg)")
    ax1.set_ylabel("2 theta (deg)")
    ax1.set_title("Run %d, Event %d" % (run, i_evt))

    #try:
    #    i_evt -= int(sys.argv[2])
    #except IndexError:
    #    i_evt -= 1
    if i_evt == 1: # if only one event exists, plot run 74, evt 1 instead
        run2 = 74
        i_evt2 = 1
    else:
        run2 = run
        i_evt2 = i_evt-1
    I2 = None
    for i_pad, pad in enumerate(pads):
        ai = pyFAI.load("calibration/pad%d.poni" % i_pad)
        img = fabio.open(
          info.data_dir+"run%d/run_%d_evt_%d_%s.tif"%(run2,run2,i_evt2,pad)).data
        if i_pad in [0,1]:
            img = img.T
        elif i_pad in [2]:
            img = np.flipud(img)
        elif i_pad in [3,4]:
            img = np.fliplr(img)
        I_i, tth, chi = ai.integrate2d(img, npt_rad=500, npt_azim=500, \
                   radial_range=(20., 110.), azimuth_range=(-180., 180.),\
                   unit='2th_deg')
        if I2 == None:
            I2 = I_i
        else:
            I2 += I_i
     
    ax2 = fig.add_subplot(2,1,2)
    ax2.imshow(np.fliplr((I2/np.max(I2)*np.max(I)).T), origin='lower', aspect='auto',\
               extent=[min(chi), max(chi), min(tth), max(tth)])
    ## draw unshocked diamond lines
    for line, line_tth in diamond_lines.iteritems():
        ax2.axhline(y=line_tth, ls='--', color='r')
        ax2.text(x=185, y=line_tth, s=line, color='r')
    ax2.set_xlabel(r"$\phi$ (deg)")
    ax2.set_ylabel(r"$2\theta$ (deg)")
    ax2.set_title("Run %d, Event %d" % (run2, i_evt2))
    #if os.path.isfile("plots/run%d.pdf"%run):
    #    os.system("mv plots/run%d.pdf plots/run%d_1.pdf"%(run,run))
    fig.savefig("plots/run%d.pdf"%run)
    plt.close()
    print "Done run %d!" % run
    sys.exit(0)
