workbook = "/Users/phsun/Dropbox/LL20 shot log - Peihao.xlsx"
data_dir = "/Users/phsun/LL20_data/"
plot_dir = "/Users/phsun/Dropbox/LL20_plots/"

timeppix = {1: 0.00939001, 2: 0.00999682} 
        # time per pixel in ns for VISAR 1 and 2

pads = ['Cspad.0', 'Cspad.1', 'Cspad.2', 'Cspad2x2.1']
hkl_list = [(1,1,1), (2,2,0), (3,1,1), (4,0,0)]

## diamond lines (tth values) for 10 keV
diamond_lines = {
        '(111)': 35.0385, '(220)': 58.8882, '(311)': 70.3985,\
        '(400)': 88.0849, '(331)': 98.5003}
chi_deg = 30.
a0 = 3.56683 # diamond at 300K, in Angstrom
E_xray = 10. # keV
wavelength = 12.3984193 / E_xray 

# data dimension
npt_rad = 500
npt_azim = 500

