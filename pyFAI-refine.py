'''
This script allows user to refine the calibration.
To run the script:
    python pyFAI-refine.py [calibration image] -c [calibrant] -f [.poni file to be refined] 
The rest should be straightforward
'''
import matplotlib
matplotlib.use('Qt4Agg')
import pyFAI, fabio
import pyFAI.calibrant
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import getopt
import sys, os

vals_name = ['dist', 'poni1', 'poni2', 'rot1', 'rot2', 'rot3']
vals_dict = {name:i for i, name in enumerate(vals_name)}

data_file = sys.argv[1]
data = fabio.open(data_file).data

npoints = 500
npoints_azim = 360
calibrant = "" 

opts, args = getopt.getopt(sys.argv[2:], "c:f:n:")
for o, a in opts:
    if o == '-c':
        calibrant = a
    elif o == '-f':
        calib_file = a
    elif o == '-n':
        npoints = int(a)

## import poni file
try:
    ai = pyFAI.load(calib_file)
except NameError:
    calib_file = raw_input("Please specify calibration file:\n")
    ai = pyFAI.load(calib_file)

## get list of tth values for calibrant
calibrant = pyFAI.calibrant.ALL_CALIBRANTS(calibrant)
calibrant.set_wavelength(ai.get_wavelength())
tth_list = np.array(calibrant.get_2th())/np.pi*180.

fig = plt.figure(figsize=(15,15))
ax1 = fig.add_subplot(2,1,1) # 2d plot
ax2 = fig.add_subplot(2,1,2, sharex=ax1) # 1d lineout
plt.subplots_adjust(hspace=0, bottom=0.2)

intensity, tth, chi = ai.integrate2d(data, npt_rad=npoints, \
                        npt_azim=npoints_azim, unit="2th_deg") 
tth_low = 1.5*tth[0]-0.5*tth[1]
tth_high = 1.5*tth[-1]-0.5*tth[-2]
chi_low = 1.5*chi[0]-0.5*chi[1]
chi_high = 1.5*chi[-1]-0.5*chi[-2]
im1 = ax1.imshow(intensity, origin='lower', aspect='auto', \
           extent=[tth_low, tth_high, chi_high, chi_low])
line2, = ax2.plot(tth, np.sum(intensity, axis=0))

## plot expected peaks
for tth in tth_list:
    ax1.axvline(x=tth, ls='--', c='k')
    ax2.axvline(x=tth, ls='--', c='k')

ax1.set_xlim(tth_low, tth_high)

## get current values
vals = [ai.get_dist()*1000., ai.get_poni1()*1000., ai.get_poni2()*1000., \
        ai.get_rot1()/np.pi*180., ai.get_rot2()/np.pi*180., ai.get_rot3()/np.pi*180.]
## print out values
print "Current values:"
print "dist=%f mm, poni1=%f mm, poni2=%f mm" % (vals[0], vals[1], vals[2])
print "rot1=%f deg, rot2=%f deg, rot3=%f deg" % (vals[3], vals[4], vals[5])

## set up sliders
axcolor = 'lightgoldenrodyellow'
vals_low = [] # to store the limits of slider
vals_high = []
ax_sliders = []
s_sliders = []
for i in xrange(6):
    ax_sliders.append(plt.axes([0.1, 0.175-0.025*i, 0.8, 0.015], \
                        axisbg=axcolor))
    if i == 0:
        vals_low.append(vals[i] * 0.95)
        vals_high.append(vals[i] * 1.05)
    elif i in [1, 2, 3, 4, 5]:
        vals_low.append(vals[i] - 2.)
        vals_high.append(vals[i] + 2.)
    s_sliders.append(Slider(ax_sliders[i], vals_name[i], \
                     vals_low[i], vals_high[i], valinit=vals[i]))

def update(val):
    ## update values
    for i in xrange(6):
        vals[i] = s_sliders[i].val
    ai.set_dist(vals[0]/1000.)
    ai.set_poni1(vals[1]/1000.)
    ai.set_poni2(vals[2]/1000.)
    ai.set_rot1(vals[3]/180.*np.pi)
    ai.set_rot2(vals[4]/180.*np.pi)
    ai.set_rot3(vals[5]/180.*np.pi)

    ## update plot
    intensity, tth, chi = ai.integrate2d(data, npt_rad=npoints, \
                            npt_azim=npoints_azim, unit="2th_deg") 
    tth_low = 1.5*tth[0]-0.5*tth[1]
    tth_high = 1.5*tth[-1]-0.5*tth[-2]
    chi_low = 1.5*chi[0]-0.5*chi[1]
    chi_high = 1.5*chi[-1]-0.5*chi[-2]
    im1.set_data(intensity)
    im1.set_extent([tth_low, tth_high, chi_high, chi_low])
    line2.set_xdata(tth)
    line2.set_ydata(np.sum(intensity, axis=0))
    ax1.set_xlim(tth_low, tth_high)
## link sliders
for i in xrange(6):
    s_sliders[i].on_changed(update)

## reset button
ax_reset = plt.axes([0.1, 0.01, 0.1, 0.02])
b_reset = Button(ax_reset, 'Reset', color=axcolor, hovercolor='0.975')
def reset(event):
    for i in xrange(6):
        s_sliders[i].reset()
b_reset.on_clicked(reset)

## save button
ax_save = plt.axes([0.45, 0.01, 0.1, 0.02])
b_save = Button(ax_save, 'Save', color=axcolor, hovercolor='0.975')
def save(event):
    os.system("cp %s %s_backup"%(calib_file, calib_file))
                    # backup calibration
    ai.save(calib_file)
    print "saved to file %s, original file backed up in %s_backup" % (\
                    calib_file, calib_file)
b_save.on_clicked(save)

## exit button
ax_exit = plt.axes([0.8, 0.01, 0.1, 0.02])
b_exit = Button(ax_exit, 'Exit', color=axcolor, hovercolor='0.975')
def exit(event):
    plt.close()
b_exit.on_clicked(exit)

plt.show()
