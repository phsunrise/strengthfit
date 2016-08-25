from PIL import Image
import sys, os
import numpy as np
import scipy.misc
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    print "must provide run number!"
    sys.exit(0)

pads = ['Cspad.0', 'Cspad.1', 'Cspad.2', \
        'Cspad2x2.1', 'Cspad2x2.2', 'Cspad2x2.4']
run = int(sys.argv[1])
directory = "/Users/phsun/LL20_data/run%d/" % run
i = 1
im_sum = []
while os.path.isfile(directory + "run_%d_evt_%d_Cspad.0.tif"%(run, i)):
    for i_pad, pad in enumerate(pads):
        im = Image.open(directory + "run_%d_evt_%d_%s.tif"%(run, i, pad))
        try:
            im_sum[i_pad] += np.array(im)
        except IndexError:
            im_sum.append(np.array(im))

    print "done event %d" % i
    i += 1

i -= 1
for i_pad, pad in enumerate(pads):
    scipy.misc.imsave("run_%d_ave_%s.tif"%(run, pad), im_sum[i_pad]*1./i) 
