# strengthfit
Right now, this module has implemented data fitting to
voigt function (fitvoigt.py)

**NOTE**: the "Qt4Agg" backend is being used in fitvoigt.py
because it is the only one I found that works fine with my 
office desktop (Mac Os X). Change it if necessary.

**FIXIT**: for some exx/ezz values, the voigt calculation 
module voigt.py cannot find a root and will raise 
RuntimeError, thereby ending the program. This error will 
be caught in future revisions.

## Installation
After cloning the repository, change the directory setting
in info.py before running any program.

* data\_dir: should point to the "LL20 data" folder in Dropbox
* plot\_dir: can be any folder where plots can be saved

## Voigt fitting
The fitvoigt.py script can run in three modes: TEST (default or -t), 
MASK (-m), and FIT (-f). To run it, provide the run number
and the mode option, *e.g.*,
'''
python fitvoigt.py -r 174 -m
'''

### TEST mode (-t)
This mode lets the user experiment with different values
for exx, ezz, and cutoff level (defined as ratio of maximum
intensity; all values below this level will be set to 0 and
shown as maximum intensity on the display).

The plotting window is detached in this mode to allow
command-line input. The input formats are:
1. x %f  
To change exx to %f

2. z %f  
To change ezz to %f

3. c %f  
To change cutoff level to %f (*e.g.*, if c 0.05 is entered,
 then all intensities less than 0.05\*max intensitiy will
 be set to 0 and shown as max intensity on display)

4. q  
To save and quit. The range for FIT mode will also be saved
as 0.8 and 1.2 of the saved exx/ezz values.

5. h  
To display help.

### MASK mode (-m)
This mode lets the user select areas of data that need to
be masked (so the intensity will be set to 0 when doing the
fit). It uses the LassoSelector in matplotlib.widgets.

The user is allowed to undo the area selections one by one,
if they are done within the session. If the mask to start
with (read from the run%d\_mask.npy) has unnecessary areas,
the user needs to quit and run 
```
python fitvoigt.py -r RUN# --dm
```
to delete the mask file.

After the data is plotted, the user can either select (on
the upper graph) areas that need to be masked, or press
one of the keys in the list below. **NOTE**: the plot window
**must be active** when the key is pressed.

* d  
To delete the last selected area.

* r  
To (temporarily) remove/recover the mask in display.

* v  
To save the mask.

* q 
To quit (without saving).

* h  
To display help.


### FIT mode (-f)
This mode lets the user specify the fitting range and then
does the fitting. The first part is very similar to TEST
mode, except that for exx and ezz, a range needs to be
specified. The format is
```
x %f %f
```
To change the range (exxlow, exxhigh), and similar for ezz.
 The two lines on display use (exxlow, ezzlow) and (exxhigh,
 ezzhigh), respectively.

When the user is satisfied with the range, "q" can be 
entered to start performing the fit. The user will be
promted to specify the number of steps (%d %d) for exx and
ezz, respectively. After the fit is done, the optimal exx
and ezz value will be selected and saved automatically 
(along with the default range 0.8 - 1.2 times the value), 
and the plot will be saved and displayed.

**NOTE**: It is recommended that fewer and larger
steps are taken at the beginning to get a rough fit, and
then the script can be run again to get a better fit.
