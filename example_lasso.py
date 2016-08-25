import numpy as np
import sys, os
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path

low, high, n = 5., 10., 5
pix = np.linspace(low, high, n, endpoint=False) + (high-low)/2./n
xx, yy = np.meshgrid(pix,pix)
pix = np.vstack((xx.ravel(), yy.ravel())).T

fig = plt.figure(figsize=(10,5))
ax1 = fig.add_subplot(1,2,1)
ax1.set_title("lasso")
im1 = ax1.imshow(xx, origin='lower', aspect='auto', \
        interpolation='nearest', extent=[low,high,low,high])

array = np.zeros((n, n))
ax2 = fig.add_subplot(1,2,2, sharex=ax1, sharey=ax1)
mask = ax2.imshow(array, origin='lower', aspect='auto',\
        vmax=1, interpolation='nearest', extent=[low,high,low,high])

ax2.set_title("mask")
ax2.set_xlim([low-1, high+1])
ax2.set_ylim([low-1, high+1])

#def updateArray(array, indices):
#    #lin = np.arange(array.size)
#    newArray = array.flatten()
#    #newArray[lin[indices]] = 1
#    newArray[indices] = 1
#    return newArray.reshape(array.shape)

def onselect(verts):
    global array, pix
    p = Path(verts)
    ind = p.contains_points(pix, radius=0)
    #array = updateArray(array, ind)
    array.ravel()[ind] = 1
    xx1 = np.copy(xx)
    xx1[array.astype(bool)] = 0.
    im1.set_data(xx1)
    mask.set_data(array)
    fig.canvas.draw()

def keypress(event):
    print "pressed", event.key
    sys.stdout.flush()

    if event.key == 'y':
        im1.set_data(yy)
        plt.draw()
    elif event.key == 'x':
        im1.set_data(xx)
        plt.draw()
    elif event.key == 'q':
        plt.close()
        sys.exit(0)

lasso = LassoSelector(ax1, onselect)
fig.canvas.mpl_connect('key_press_event', keypress)
plt.show()
