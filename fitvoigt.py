import matplotlib
matplotlib.use("Qt4Agg") # for mac osx
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os, sys, getopt, glob, re
import pickle
import info, voigt

mode = 'TEST' # default mode
try:
    opts, args = getopt.getopt(sys.argv[1:], "r:tmf",\
            ['dt', 'dm'])
except getopt.GetOptError as err:
    print str(err)
for o, a in opts:
    if o == '-r':
        run = int(a)
    elif o == '-f':
        mode = 'FIT'
    elif o == '-m':
        mode = 'MASK' # mask mode, will let user select mask region 
    elif o == '-t':
        mode = 'TEST'  # test mode, will let user select e values
    elif o == '--dt':
        fname = info.data_dir+"parameters/run%d.pickle"%run
        if os.path.isfile(fname):
            val = raw_input("Delete %s? Enter y to confirm: "%fname)
            if val == 'y':
                os.system("rm %s"%fname)
                print "File %s deleted!" % fname
            else:
                print "not deleted"
        else:
            print "File %s does not exist!" % fname
        sys.exit(0)
    elif o == '--dm':
        fname = info.data_dir+"parameters/run%d_mask.npy"%run
        if os.path.isfile(fname):
            val = raw_input("Delete %s? Enter y to confirm: "%fname)
            if val == 'y':
                os.system("rm %s"%fname)
                print "File %s deleted!" % fname
            else:
                print "not deleted"
        else:
            print "File %s does not exist!" % fname
        sys.exit(0)

## import parameters
from parameters import *

## read data
try:
    fname = glob.glob(info.data_dir+"converted/run%d_evt*_data.npy"%run)[0]
    I = np.load(fname)
    i_evt = int(re.search('(?<=evt)\d+', fname).group(0))
    maxI = np.max(I)
    print "Data loaded"
except IndexError:
    print "data for run %d does not exist" % run
    sys.exit(0)
coord = np.load(info.data_dir+"converted/coordinates.npz")
tth = coord['tth']
phi = coord['phi']
tth_grid = coord['tth_grid']
phi_grid = coord['phi_grid']

## try to read parameters from file
try:
    with open(info.data_dir+"parameters/run%d.pickle"%run, 'r') as f:
        params = pickle.load(f)
        exx = params['exx']
        ezz = params['ezz']
        exxlow = params['exxlow']
        exxhigh = params['exxhigh']
        ezzlow = params['ezzlow']
        ezzhigh = params['ezzhigh']
        cutoff_level = params['cutoff']
        print "Loaded config file"
except IOError:
    print "Parameter file does not exist; running TEST mode..."
    mode = 'TEST'
    exx = 0.
    ezz = 0.
    cutoff_level = 0.05

try:
    mask_orig = np.load(info.data_dir+"parameters/run%d_mask.npy"%run)
    print "Loaded mask file"
except IOError:
    if mode == 'FIT':
        print "Mask file does not exist; running MASK mode..."
        mode = 'MASK'
    mask_orig = []

## TEST mode
if mode == 'TEST':
    print "TEST mode"
    fig = plt.figure(figsize=(12,8))
    ax1 = fig.add_subplot(1,1,1)

    cutoff = cutoff_level*maxI
    I1 = np.copy(I)
    I1[I1<cutoff] = maxI
    im1 = ax1.imshow(I1, origin='lower', aspect='auto',\
               extent=[min(phi), max(phi), min(tth), max(tth)])
    divider1 = make_axes_locatable(ax1)
    cax1 = divider1.append_axes("right", size="10%", pad=0.05)
    lines = [] 
    for h, k, l in hkl_list: 
        d0 = a0 / np.sqrt(h**2+k**2+l**2)
        _phi_array, tth_voigt = voigt.voigt(chi_deg, exx, exx, ezz,\
                d0, wavelength, phi_array=phi/180.*np.pi, tol=1.e-12)
        lines.append(ax1.plot(phi, tth_voigt/np.pi*180., \
                color='r', ls='--'))
    ax1.set_xlim(min(phi), max(phi))
    ax1.set_xlabel(r"$\phi$ (deg)")
    ax1.set_ylabel(r"$2\theta$ (deg)")
    ax1.set_title("Run %d, TEST mode"%run)
    fig.colorbar(im1, cax=cax1)
    #fig.savefig(info.plot_dir+"run%d_sat_%.2f_ex_%.3f_ez_%.3f.pdf"%(run,\
    #                                cutoff_level, exx, ezz))
    plt.ion()
    plt.show()

    print "Enter new parameters, or q to save and quit, or h for help"
    while True:
        print "  Current values:"
        print "\tcutoff level = %f" % cutoff_level,
        print "exx = %f" % exx,
        print "ezz = %f" % ezz
        val = raw_input("")
        if val == 'q':
            plt.close()
            fname = info.data_dir+"parameters/run%d.pickle"%run
            with open(fname, 'w') as f:
                pickle.dump({'exx':exx, 'ezz':ezz, \
                        'exxlow':exx*0.8, 'exxhigh':exx*1.2,\
                        'ezzlow':ezz*0.8, 'ezzhigh':ezz*1.2,\
                        'cutoff':cutoff_level}, f)
                print "Saved to file %s" % fname
            sys.exit(0)

        elif val == 'h':
            print "Parameter format:"
            print "\"c %f\" to change cutoff level (as ratio of max intensity)"
            print "\"x %f\" to change exx" 
            print "\"z %f\" to change ezz" 
            continue
        
        elif len(val) >= 3 and val[0] == 'c':
            try:
                _, cutoff_level = val.split()
                cutoff_level = float(cutoff_level)
                cutoff = cutoff_level*maxI
                I1 = np.copy(I)
                I1[I1<cutoff] = maxI
                im1.set_data(I1)
                plt.draw()
                continue
            except ValueError:
                print "Invalid input!"
                continue

        elif len(val) >= 3 and val[0] == 'x':
            try:
                _, exx = val.split()
                exx = float(exx)
                for i_hkl, (h, k, l) in enumerate(hkl_list): 
                    d0 = a0 / np.sqrt(h**2+k**2+l**2)
                    _phi_array, tth_voigt = voigt.voigt(chi_deg, exx, exx, ezz,\
                            d0, wavelength, phi_array=phi/180.*np.pi, tol=1.e-12)
                    lines[i_hkl].pop(0).remove()
                    lines[i_hkl] = ax1.plot(phi, tth_voigt/np.pi*180., \
                            color='r', ls='--')
                continue
            except ValueError:
                print "Invalid input!"
                continue
               
        elif len(val) >= 3 and val[0] == 'z':
            try:
                _, ezz = val.split()
                ezz = float(ezz)
                for i_hkl, (h, k, l) in enumerate(hkl_list): 
                    d0 = a0 / np.sqrt(h**2+k**2+l**2)
                    _phi_array, tth_voigt = voigt.voigt(chi_deg, exx, exx, ezz,\
                            d0, wavelength, phi_array=phi/180.*np.pi, tol=1.e-12)
                    lines[i_hkl].pop(0).remove()
                    lines[i_hkl] = ax1.plot(phi, tth_voigt/np.pi*180., \
                            color='r', ls='--')
                continue
            except ValueError:
                print "Invalid input!"
                continue
                
        else:
            print "Invalid input!"
            continue

## MASK mode
elif mode == 'MASK':
    print "MASK mode"
    print "If the starting mask has unnecessary regions, quit and use --dm to remove the mask file"
    from matplotlib.widgets import LassoSelector
    from matplotlib.path import Path

    if mask_orig == []:
        mask_orig = np.zeros_like(I).astype(bool)
    mask = mask_orig
    nomask = False # flag to turn off mask in ax1 display

    ## first make the plots
    fig = plt.figure(figsize=(12,16))
    ax1 = fig.add_subplot(2,1,1)
    cutoff = cutoff_level*maxI
    I1 = np.copy(I)
    if not nomask:
        I1[mask] = 0.
    I1[I1<cutoff] = maxI
    im1 = ax1.imshow(I1, origin='lower', aspect='auto',\
               extent=[min(phi), max(phi), min(tth), max(tth)])
    divider1 = make_axes_locatable(ax1)
    cax1 = divider1.append_axes("right", size="10%", pad=0.05)
    for h, k, l in hkl_list: 
        d0 = a0 / np.sqrt(h**2+k**2+l**2)
        _phi_array, tth_voigt = voigt.voigt(chi_deg, exx, exx, ezz,\
                d0, wavelength, phi_array=phi/180.*np.pi, tol=1.e-12)
        ax1.plot(phi, tth_voigt/np.pi*180., color='r', ls='--')
    ax1.set_xlim(min(phi), max(phi))
    ax1.set_xlabel(r"$\phi$ (deg)")
    ax1.set_ylabel(r"$2\theta$ (deg)")
    fig.colorbar(im1, cax=cax1)

    ax2 = fig.add_subplot(2,1,2, sharex=ax1, sharey=ax1)
    im2 = ax2.imshow(mask, origin='lower', aspect='auto',\
           extent=[min(phi), max(phi), min(tth), max(tth)],\
           vmin=0, vmax=1, cmap='Greys')
    ax2.set_xlabel(r"$\phi$ (deg)")
    ax2.set_ylabel(r"$2\theta$ (deg)")

    ax1.set_title("Select lasso region. Press h for help (plot window needs to be active).")
    ax2.set_title("Mask")

    pixels = np.vstack((phi_grid.ravel(), tth_grid.ravel())).T
    regions = []
    def onselect(verts):
        global nomask, mask, I1, regions
        nomask = False # selecting lasso region will automatically turn on mask 
        lassopath = Path(verts)
        newregion = lassopath.contains_points(pixels, radius=0.).reshape(I1.shape)
        regions.append(newregion)
        mask = np.logical_or(mask, newregion)
        I1[mask] = maxI 
        im1.set_data(I1)
        im2.set_data(mask)
        plt.draw()

    def press(event):
        global nomask, mask, I1, regions
        print "pressed", event.key
        sys.stdout.flush()

        if event.key == 'h':
            print "Press d to delete the last region, r to remove/recover the mask, v to save, or q to quit"
            print "Plot window needs to be active when key is pressed"

        elif event.key == 'r':
            if nomask: # if currently not masked
                print "turning mask on"
                nomask = False
                I1[mask] = maxI
                im1.set_data(I1)
                plt.draw()
            else: # if currently masked
                print "turning mask off"
                nomask = True
                I1 = np.copy(I)
                I1[I1<cutoff] = maxI
                im1.set_data(I1)
                plt.draw()

        elif event.key == 'd':
            if regions == []:
                print "No region selected!"
            else:
                regions.pop()
                ## update image and mask
                mask = mask_orig
                for region in regions:
                    mask = np.logical_or(mask, region) 
                mask = mask.astype(bool)
                I1 = np.copy(I)
                I1[mask] = maxI
                I1[I1<cutoff] = maxI
                im1.set_data(I1)
                im2.set_data(mask)
                plt.draw()

        elif event.key == 'v':
            fout = info.data_dir+"parameters/run%d_mask.npy"%run
            np.save(fout, mask)
            print "Mask saved to %s" % fout

        elif event.key == 'q':
            sys.exit(0)

    lasso = LassoSelector(ax1, onselect)
    fig.canvas.mpl_connect('key_press_event', press)
    plt.show()


## FIT mode
elif mode == 'FIT':
    print "FIT mode"
    fig = plt.figure(figsize=(18,12))
    ax1 = fig.add_subplot(1,1,1)

    cutoff = cutoff_level*maxI
    I1 = np.copy(I)
    I1[mask_orig] = maxI
    I1[I1<cutoff] = maxI
    im1 = ax1.imshow(I1, origin='lower', aspect='auto',\
               extent=[min(phi), max(phi), min(tth), max(tth)])
    divider1 = make_axes_locatable(ax1)
    cax1 = divider1.append_axes("right", size="10%", pad=0.05)
    lines_low = []
    lines_high = []
    for h, k, l in hkl_list: 
        d0 = a0 / np.sqrt(h**2+k**2+l**2)
        _phi_array, tth_voigt = voigt.voigt(chi_deg, exxlow, exxlow, ezzlow,\
                d0, wavelength, phi_array=phi/180.*np.pi, tol=1.e-12)
        lines_low.append(ax1.plot(phi, tth_voigt/np.pi*180., \
                color='c', ls='--'))
        _phi_array, tth_voigt = voigt.voigt(chi_deg, exxhigh, exxhigh, ezzhigh,\
                d0, wavelength, phi_array=phi/180.*np.pi, tol=1.e-12)
        lines_high.append(ax1.plot(phi, tth_voigt/np.pi*180., \
                color='y', ls='--'))
    ax1.set_xlim(min(phi), max(phi))
    ax1.set_xlabel(r"$\phi$ (deg)")
    ax1.set_ylabel(r"$2\theta$ (deg)")
    ax1.set_title("Run %d, FIT mode"%run)
    fig.colorbar(im1, cax=cax1)
    #fig.savefig(info.plot_dir+"run%d_sat_%.2f_ex_%.3f_ez_%.3f.pdf"%(run,\
    #                                cutoff_level, exx, ezz))
    plt.ion()
    plt.show()

    print "Enter new parameters, or q to save and continue, or h for help"
    while True:
        print "  Current values:"
        print "\tcutoff level = %f" % cutoff_level,
        print "exx = %f" % exx,
        print "ezz = %f" % ezz
        print "\t exx range: %f to %f" % (exxlow, exxhigh)
        print "\t ezz range: %f to %f" % (ezzlow, ezzhigh)
        val = raw_input("")
        if val == 'exit':
            print "exiting..."
            sys.exit(0)

        elif val == 'q':
            plt.close()
            fname = info.data_dir+"parameters/run%d.pickle"%run
            with open(fname, 'w') as f:
                pickle.dump({'exx':exx, 'ezz':ezz, \
                        'exxlow':exxlow, 'exxhigh':exxhigh,\
                        'ezzlow':ezzlow, 'ezzhigh':ezzhigh,\
                        'cutoff':cutoff_level}, f)
                print "Saved to file %s" % fname
            break

        elif val == 'h':
            print "Parameter format:"
            print "\"c %f\" to change cutoff level (as ratio of max intensity)"
            print "\"x %f %f\" to change exx range" 
            print "\"z %f %f\" to change ezz range" 
            print "Enter \"exit\" to exit the program"
            continue
        
        elif len(val) >= 3 and val[0] == 'c':
            try:
                _, cutoff_level = val.split()
                cutoff_level = float(cutoff_level)
                cutoff = cutoff_level*maxI
                I1 = np.copy(I)
                I1[I1<cutoff] = maxI
                im1.set_data(I1)
                plt.draw()
                continue
            except ValueError:
                print "Invalid input!"
                continue

        elif len(val) >= 3 and val[0] == 'x':
            try:
                _, exxlow, exxhigh = val.split()
                _exxlownhigh = float(exxlow), float(exxhigh)
                exxlow, exxhigh = min(_exxlownhigh), max(_exxlownhigh) 
                for i_hkl, (h, k, l) in enumerate(hkl_list): 
                    d0 = a0 / np.sqrt(h**2+k**2+l**2)
                    lines_low[i_hkl].pop(0).remove()
                    lines_high[i_hkl].pop(0).remove()
                    _phi_array, tth_voigt = voigt.voigt(chi_deg, exxlow, exxlow, ezzlow,\
                            d0, wavelength, phi_array=phi/180.*np.pi, tol=1.e-12)
                    lines_low[i_hkl] = ax1.plot(phi, tth_voigt/np.pi*180., \
                            color='c', ls='--')
                    _phi_array, tth_voigt = voigt.voigt(chi_deg, exxhigh, exxhigh, ezzhigh,\
                            d0, wavelength, phi_array=phi/180.*np.pi, tol=1.e-12)
                    lines_high[i_hkl] = ax1.plot(phi, tth_voigt/np.pi*180., \
                            color='y', ls='--')
                continue
            except ValueError:
                print "Invalid input!"
                continue
               
        elif len(val) >= 3 and val[0] == 'z':
            try:
                _, ezzlow, ezzhigh = val.split()
                _ezzlownhigh = float(ezzlow), float(ezzhigh)
                ezzlow, ezzhigh = min(_ezzlownhigh), max(_ezzlownhigh) 
                for i_hkl, (h, k, l) in enumerate(hkl_list): 
                    d0 = a0 / np.sqrt(h**2+k**2+l**2)
                    lines_low[i_hkl].pop(0).remove()
                    lines_high[i_hkl].pop(0).remove()
                    _phi_array, tth_voigt = voigt.voigt(chi_deg, exxlow, exxlow, ezzlow,\
                            d0, wavelength, phi_array=phi/180.*np.pi, tol=1.e-12)
                    lines_low[i_hkl] = ax1.plot(phi, tth_voigt/np.pi*180., \
                            color='c', ls='--')
                    _phi_array, tth_voigt = voigt.voigt(chi_deg, exxhigh, exxhigh, ezzhigh,\
                            d0, wavelength, phi_array=phi/180.*np.pi, tol=1.e-12)
                    lines_high[i_hkl] = ax1.plot(phi, tth_voigt/np.pi*180., \
                            color='y', ls='--')
                continue
            except ValueError:
                print "Invalid input!"
                continue
                
        else:
            print "Invalid input!"
            continue

    plt.ioff()
    plt.close()
    
    fig = plt.figure(figsize=(9, 12))
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)
    cutoff = cutoff_level*maxI
    I1 = np.copy(I)
    I1[I1<cutoff] = 0.
    I1[mask_orig] = 0.
    im1 = ax1.imshow(I1, origin='lower', aspect='auto',\
               extent=[min(phi), max(phi), min(tth), max(tth)])
    divider1 = make_axes_locatable(ax1)
    cax1 = divider1.append_axes("right", size="10%", pad=0.05)
    ### draw unshocked diamond lines
    #for line, line_tth in diamond_lines.iteritems():
    #    ax1.axhline(y=line_tth, ls='--', color='r')
    #    ax1.text(x=185, y=line_tth, s=line, color='r')

    while True:
        val = raw_input("Number of steps (x and z): ")
        try:
            exx_n, ezz_n = val.split()
            exx_n = int(exx_n)
            ezz_n = int(ezz_n)
            break
        except ValueError:
            print "Invalid input! Try again"
            continue

    exx_array = np.linspace(exxlow, exxhigh, exx_n)
    ezz_array = np.linspace(ezzlow, ezzhigh, ezz_n)
    weightedsum = np.zeros((len(exx_array), len(ezz_array)))
    optim_res = [float('inf'), 0., 0.] # array to store the optimum result
    for i_exx, exx in enumerate(exx_array):
        for i_ezz, ezz in enumerate(ezz_array):
            dist_grid = []
            ## get all three lines, calculate weighted squared distance
            for h, k, l in hkl_list: 
                d0 = a0 / np.sqrt(h**2+k**2+l**2)
                _phi_array, tth_voigt = voigt.voigt(chi_deg, exx, exx, ezz, \
                        d0, wavelength, phi_array=phi/180.*np.pi, tol=1.e-12)
                dist_grid.append(np.abs(tth_grid/180.*np.pi-np.tile(tth_voigt, [npt_rad,1])))
            ## for each point, find the minimal distance from the three tth_voigt values
            dist_grid = np.array(dist_grid)
            dist_grid = np.amin(dist_grid, axis=0)
            #im1 = ax1.imshow(dist_grid, origin='lower')
            #fig.colorbar(im1, cax=cax1)
            #fig.savefig("/Users/phsun/Dropbox/LL20_plots/run%d.pdf"%run)
            #plt.close()
            #sys.exit(0)
            weightedsum[i_ezz, i_exx] = np.sum(dist_grid**2 * I1)
            if weightedsum[i_ezz, i_exx] < optim_res[0]:
                optim_res = [weightedsum[i_ezz, i_exx], exx, ezz]
            
        print "done exx =", exx

    ## now save to parameter file and plot the optimum  
    print "Best fit: exx=%.4f, ezz=%.4f, weightedsum=%.1f" % (optim_res[1], 
                                                    optim_res[2], optim_res[0])
    fname = info.data_dir+"parameters/run%d.pickle"%run
    with open(fname, 'w') as f:
        pickle.dump({'exx':optim_res[1], 'ezz':optim_res[2],\
                'exxlow':optim_res[1]*0.8, 'exxhigh':optim_res[1]*1.2,\
                'ezzlow':optim_res[2]*0.8, 'ezzhigh':optim_res[2]*1.2,\
                'cutoff':cutoff_level}, f)
        print "Saved parameters to %s" % fname
    for h, k, l in hkl_list: 
        d0 = a0 / np.sqrt(h**2+k**2+l**2)
        _phi_array, tth_voigt = voigt.voigt(chi_deg, \
                optim_res[1], optim_res[1], optim_res[2], \
                d0, wavelength, phi_array=phi/180.*np.pi, tol=1.e-12)

        ax1.plot(phi, tth_voigt/np.pi*180., color='r', ls='--')

    ax1.set_xlabel(r"$\phi$ (deg)")
    ax1.set_ylabel(r"$2\theta$ (deg)")
    ax1.set_xlim(-180., 180.)
    ax1.set_ylim(20., 110.)
    ax1.set_title("Run %d, Evt %d, exx=%.4f, ezz=%.4f" % \
                    (run, i_evt, optim_res[1], optim_res[2]))
    fig.colorbar(im1, cax=cax1)

    im2 = ax2.imshow(weightedsum, origin='lower', aspect='auto', \
               extent=[exxlow, exxhigh, ezzlow, ezzhigh],\
               interpolation='nearest')
    divider2 = make_axes_locatable(ax2)
    cax2 = divider2.append_axes("right", size="10%", pad=0.05)
    ax2.set_xlabel(r"$\epsilon_{xx}$")
    ax2.set_ylabel(r"$\epsilon_{zz}$")
    fig.colorbar(im2, cax=cax2)

    fig.savefig(info.plot_dir+"run%d_cutoff_%.3f.pdf"%(run, cutoff_level))
    plt.show()
    plt.close()

    with open(info.plot_dir+
          "run%d_cutoff_%.3f.npz"%(run, cutoff_level), 'w') as outfile:
        np.savez(outfile,\
          weightedsum=weightedsum, exx_array=exx_array, ezz_array=ezz_array)
