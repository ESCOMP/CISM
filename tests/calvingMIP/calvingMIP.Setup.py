#!/usr/bin/env python

# This script sets up experiments for the Calving Model Intercomparison (CalvingMIP).
# More information about CalvingMIP can be found on github:
# https://github.com/JRowanJordan/CalvingMIP/wiki.
# For now, the script sets up Experiments 1 to 5 of CalvingMIP Phase 1.

import sys, os
import shutil
import fileinput
import numpy as np
from netCDF4 import Dataset
from configparser import ConfigParser
#from optparse import OptionParser
import argparse
parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# small helper function so argparse will understand unsigned integers
def unsigned_int(x):
    x = int(x)
    if x < 1:
        raise argparse.ArgumentTypeError("This argument is an unsigned int type! Should be an integer greater than zero.")
    return x

#############
# Constants #
#############

R = 800.e3   # (m) size of the circular and Thule domains; x and y range over (-R, R)

# For the circular domain (Experiments 1 and 2)
Bc_circ = 900.    # (m) maximum bed topography elevation at the center of the domain
Bl_circ = -2000.  # (m) minimum bed topography elevation at the edge of the domain

# For the Thule domain (Experiments 3, 4 and 5)
Bc_thule = 900.    # (m) maximum bed topography elevation at the center of the domain
Bl_thule = -2000.  # (m) minimum bed topography elevation at the edge of the domain
Ba_thule = 1100


# For initialization
#WHL - Spinup uses which_ho_calving_front = 0, marine_margin = 5 (calving mask)
#      CalvingMIP runs use which_ho_calving_front = 1, marine_margin = 6
initThickness = 0.   # initial uniform ice thickness (m)
rcalve = 750000.     # radial distance of calving front from center (m)
restartfreqSpinup = 5000.    # frequency at which restart file is written (yr)

# Physical constants
# See here: https://github.com/JRowanJordan/CalvingMIP/wiki/Physical-constants-and-assumptions
# Notes on units:
# For the rate factor A:
# Protocol has A = 10^{-9} kPa^{-3} a^{-1}
# This translates to 10^{-18} Pa^{-3} a^{-1} in CISM units
# For basal sliding:
# Protocol has C = 0.001 m a^{-1} kPa^{-3} in the relation u_b = C * tau_b^m with m = 3
# CISM powerlaw sliding uses the relation tau^b = beta * u_b^{1/m}
# Thus, beta = (1/C)^{1/m} = 10^4 Pa m^{-1/3} a^{1/3}

g = 9.81        # gravitational acceleration (m s-2)
accum = 0.3     # Surface mass balance (m a-1)
rhoi = 917.     # Ice density (kg m-3)
rhoo = 1028.    # Sea water density (kg m-3)
A = 2.9377e-18  # Ice rate factor (Pa-3 a-1)
n = 3           # Flow law stress exponent
C = 10000.      # Basal slipperiness (m a-1 Pa-3)
m = 3.          # Sliding law stress exponent
sec2yr = 31556926  # Seconds in a year (s)



####################################
# Function used later in the code #
####################################


#This function computes the circular bed topography
def computeBedCircular(x,y,R,Bc,Bl):
    # Input:
    # x = grid x-coordinate in km
    # y = grid y-coordinate in km
    # R = size of the circular domain (m)
    # Bc = maximum bed topography elevation at the center of the domain (m)
    # Bl = minimum bed topography elevation at the edge of the domain (m)
    #
    # Output:
    # B = bed elevation (m)

    r = np.sqrt(x**2 + y**2)
    B = Bc - (Bc - Bl)*r**2/R**2
    return B


#This function computes the Thule bed topography
def computeBedThule(x,y,R,Bc,Bl,Ba):
    # Input:
    # x = grid x-coordinate in km
    # y = grid y-coordinate in km
    # R = size of the circular domain (m)
    # Bc = maximum bed topography elevation at the center of the domain (m)
    # Bl = minimum bed topography elevation at the edge of the domain (m)
    # Ba = minimum bed topography elevation at the edge of the domain (m)
    #
    # Output:
    # B = Bed elevation (m)
    
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y,x)
    l = R-np.cos(2.*theta)*R/2
    a = Bc - (Bc - Bl)*r**2/R**2
    B = Ba*np.cos(3.*np.pi*r/l) + a
    return B


# This function creates a target file identical to a source file.
# Used to create a Thule input file which is the same as the circular input file
# (except for topg, which is to be edited).
# Copied from stackoverflow: https://stackoverflow.com/questions/13936563/copy-netcdf-file-using-python

def createFileFromSource(src_file, trg_file):
    src = Dataset(src_file)
    trg = Dataset(trg_file, mode='w')

    # Create the dimensions of the file
    for name, dim in src.dimensions.items():
        trg.createDimension(name, len(dim) if not dim.isunlimited() else None)

    # Copy the global attributes
    trg.setncatts({a:src.getncattr(a) for a in src.ncattrs()})

    # Create the variables in the file
    for name, var in src.variables.items():
        trg.createVariable(name, var.dtype, var.dimensions)

        # Copy the variable attributes
        trg.variables[name].setncatts({a:var.getncattr(a) for a in var.ncattrs()})

        # Copy the variables values (as 'f4' eventually)
        trg.variables[name][:] = src.variables[name][:]

    # Save the file
    trg.close()


#############
# Main Code #
#############

# Parse options.

parser.add_argument('-c', '--configfile',   type=str,   default='calvingMIP.config.template', help="config file template")
parser.add_argument('-e', '--executable',   type=str,   default='cism_driver', help="path to the CISM executable")
parser.add_argument('-x', '--experiment',   type=str,   default= 'all',   help="CalvingMIP experiment(s) to set up")
parser.add_argument('-t', '--timestep',     type=float, default= 1.0,     help="time step (yr)")
parser.add_argument('-r', '--resolution',   type=float, default= 5000.,   help="grid resolution (m)")
parser.add_argument('-v', '--vertlevels',   type=int,   default= 5,       help="no. of vertical levels")
parser.add_argument('-a', '--approximation',type=str,   default= 'DIVA',  help="Stokes approximation (SSA, DIVA, BP)")
parser.add_argument('-y', '--yearsSpinup',  type=float, default= 10000.,  help="length of spinup run (yr)")
parser.add_argument('-fh', '--outputfreqhi',type=float, default= 10.,     help="high output frequency (yr)")
parser.add_argument('-fm', '--outputfreqmd',type=float, default= 100.,    help="medium output frequency (yr)")
parser.add_argument('-fl', '--outputfreqlo',type=float, default= 1000.,   help="low output frequency (yr)")

args = parser.parse_args()

if args.experiment == 'all':
    experiments = ['SpinupCircular', 'SpinupThule', 'Experiment1', 'Experiment2', 'Experiment3', 'Experiment4', 'Experiment5' ]
    print('Setting up all the CalvingMIP experiments of phase 1')
elif args.experiment == 'allSpinup':
    experiments = ['SpinupCircular', 'SpinupThule']
elif args.experiment == 'allExpt':
    experiments = ['Experiment1', 'Experiment2', 'Experiment3', 'Experiment4', 'Experiment5']
    print('Setting up all the CalvingMIP experiments, excluding spinups')
    #Note: These experiments will not run unless each Spinup directory exists and contains a restart file from an earlier run
elif args.experiment in ['SpinupCircular', 'SpinupThule', 'Experiment1', 'Experiment2', 'Experiment3', 'Experiment4', 'Experiment5']: 
    experiments = [args.experiment]
    print( 'Setting up experiment', args.experiment)
else:
    sys.exit('Please specify experiment(s) from this list: all, allSpinup, allExp, Experiment1, Experiment2, Experiment3, Experiment4, Experiment5') 

# If the executable is not called 'cism_driver', then create a link called 'cism_driver' to the executable.
# Each subdirectory will link to cism_driver in the main directory.
if args.executable != 'cism_driver':
    # Remove the existing link, if present.
    os.unlink('cism_driver')
    # Make the new link.
    os.symlink(args.executable, 'cism_driver')

# Set grid resolution.
if R%args.resolution==0:
    dx = args.resolution
    dy = args.resolution
else:
    sys.exit('Your resolution should be a divider of the domain size of 800 km')

if args.vertlevels >= 2:
    nz = args.vertlevels
else:
    sys.exit('Error: must have at least 2 vertical levels')

print( 'CalvingMIP grid resolution (m) =', args.resolution)
print( 'Number of vertical levels =', nz)

# Set number of grid cells in each direction.
# E.g., if R = 800 km with dx = dy = 10 km, we have 161 cells in each direction
nx = 2*int(R/dx) + 1
ny = 2*int(R/dy) + 1

print('grid dimension in x:',nx)
print('grid dimension in y:',ny)

#sys.exit('coucou')

# Copy the config template to a new master config file.
# Later, the master file will be tailored to each experiment.
masterConfigFile = 'calvingMIP.config'

try:
    shutil.copy(args.configfile, masterConfigFile)
except:
    sys.exit('Could not copy', args.configfile)

print( 'Creating master config file', masterConfigFile)

# Read the master config file.
config = ConfigParser()
config.read(masterConfigFile)

# Set the grid variables
config.set('grid', 'ewn', str(nx))
config.set('grid', 'nsn', str(ny))
config.set('grid', 'upn', str(nz))
config.set('grid', 'dew', str(dx))
config.set('grid', 'dns', str(dy))

# Set the timestep
config.set('time', 'dt',      str(args.timestep))
config.set('time', 'dt_diag', str(args.timestep))

# Set the diagnostic cell
# By default, this cell lies at the center of the x-axis, near the
# northern boundary of the domain.
idiag = int(nx/2) + 1
jdiag = int(idiag*31./16.) - 1
config.set('time', 'idiag', str(idiag))
config.set('time', 'jdiag', str(jdiag))
print('idiag:',idiag)
print('jdiag:',jdiag)

# Set the Stokes approximation
if args.approximation == 'SSA':
    which_ho_approx = 1
    print( 'Using SSA velocity solver')
elif args.approximation == 'DIVA':
    which_ho_approx = 4
    print( 'Using DIVA velocity solver')
elif args.approximation == 'BP':
    which_ho_approx = 2
    print( 'Using Blatter-Pattyn velocity solver')
else:
    which_ho_approx = 4
    print( 'Defaulting to DIVA velocity solver')

config.set('ho_options', 'which_ho_approx', str(which_ho_approx))


# Note: Only one option is supported for the basal sliding law:
#       a Weertman-style power law of the form tau_b = beta * u_b^{1/m}

# Set the number of years for the spin-up
yearsSpinup = float(args.yearsSpinup)
print('years of spinup =', yearsSpinup)
print('spinup restart frequency =', restartfreqSpinup)
print('high output frequency (scalars) =', args.outputfreqhi)
print('medium output frequency (Expt1 to Expt4) =', args.outputfreqmd)
print('low output frequency (Spinup, Expt5) =', args.outputfreqlo)

# Write to the master config file.
with open(masterConfigFile, 'w') as configfile:
    config.write(configfile)


# Create the netCDF input files: calvingMIP.inputCircular.nc and calvingMIP.inputThule.nc

# First create an input file for the circular domain.
# Later, we will copy and edit as needed for the Thule domain.

try:
    parser = ConfigParser()
    parser.read(args.configfile)
    inputfile = parser.get('CF input', 'name')
except:
    sys.exit('Error parsing ' + args.configfile)

#Put 'Circular' in the filename
string1 = 'input'
string2 = 'inputCircular'
inputfileCircular = inputfile.replace(string1,string2)

print('Creating input file', inputfileCircular)
ncfile = Dataset(inputfileCircular, 'w')

# Create dimensions.
# Note: (x0,y0) = staggered (velocity) grid
#       (x1,y1) = unstaggered (scalar) grid
ncfile.createDimension('time', 1)
ncfile.createDimension('x1',  nx)
ncfile.createDimension('y1',  ny)
ncfile.createDimension('x0',  nx-1)
ncfile.createDimension('y0',  ny-1)
ncfile.createDimension('level',         nz)
ncfile.createDimension('staglevel',     nz-1)
ncfile.createDimension('stagwbndlevel', nz+1)

# Create time and grid variables.
# Note: (x1,y1) are loadable and need to be in the input file.
#       (x0,y0) are not loadable, but are derived in CISM from (x1,y1).

ncfile.createVariable('time','f4',('time',))[:] = [0]
x1 = ncfile.createVariable('x1','f4',('x1',))
y1 = ncfile.createVariable('y1','f4',('y1',))
x0 = ncfile.createVariable('x0','f4',('x0',))
y0 = ncfile.createVariable('y0','f4',('y0',))

# Create 2D input fields
thk  = ncfile.createVariable('thk',  'f4', ('time','y1','x1'))
topg = ncfile.createVariable('topg', 'f4', ('time','y1','x1'))
acab = ncfile.createVariable('acab', 'f4', ('time','y1','x1'))
calving_mask = ncfile.createVariable('calving_mask', 'i4', ('time','y1','x1'))

# Compute x and y on each grid.
# The origin (x = y = 0) is placed at the center of the domain.

x = np.arange(-R,R+dx,dx,dtype='float32')   # x = -R, -R+dx,..., -dx, 0, dx,..., R-dx, R
y = np.arange(-R,R+dy,dy,dtype='float32')   # y = -R, -R+dx,..., -dx, 0, dx,..., R-dx, R

#print('x=',x[-4::])
#print('y=',y[-4::])
#sys.exit('coucou')

x1[:] = x[:]
y1[:] = y[:]

x0[:] = x[:-1] + dx/2.   # x = -R+dx/2,..., -dx/2, dx/2,..., R-dx/2
y0[:] = y[:-1] + dy/2.   # x = -R+dy/2,..., -dy/2, dy/2,..., R-dy/2

# Set SMB
acab[:,:,:] = accum

# Set initial thickness and calving mask
thk[:,:,:] = 0.
calving_mask[:,:,:] = 1

for i in range(1,nx-1):
    for j in range(1,ny-1):
        # Find the distance r from the origin to the most distant cell vertex.
        # If r < rcalve, then put ice in the cell and set calving_mask = 0; if not, then set H = 0 and calving_mask = 1.
        rsw = np.sqrt(x0[i-1]**2 + y0[j-1]**2)
        rnw = np.sqrt(x0[i-1]**2 + y0[j]**2)
        rne = np.sqrt(x0[i]**2 + y0[j]**2)
        rse = np.sqrt(x0[i]**2 + y0[j-1]**2)
        rmax = max(rsw,rnw,rne,rse)  # distance from the origin to the most distant vertex
        if rmax <= rcalve:
            thk[:,j,i] = initThickness
            calving_mask[:,j,i] = 0


# Set circular bed topography
topg[:,:,:] = 0.
for i in range(0,nx):
    for j in range(0,ny):
        topg[:,j,i] = computeBedCircular(x1[i], y1[j], R, Bc_circ, Bl_circ)


# Close the file 
ncfile.close()

# Create a Thule input file using the circular input file as the source.
string1 = 'Circular'
string2 = 'Thule'
inputfileThule = inputfileCircular.replace(string1,string2)
print('Creating input file', inputfileThule)

createFileFromSource(inputfileCircular, inputfileThule)
ncfile = Dataset(inputfileThule, 'r+')

# Set Thule bed topography
# Note: The commented lines generate an error, but the uncommented lines give the desired result.
#topg[:,:,:] = 0.
#for i in range(0,nx):
#    for j in range(0,ny):
#        topg[:,j,i] = computeBedThule(x1[i], y1[j], R, Bc_thule, Bl_thule, Ba_thule)
ncfile['topg'][:,:,:] = 0.
for i in range(0,nx):
    for j in range(0,ny):
        ncfile['topg'][:,j,i] = computeBedThule(ncfile['x1'][i], ncfile['y1'][j], R, Bc_thule, Bl_thule, Ba_thule)


# Close the file 
ncfile.close()

print('Experiments:', experiments)

# Loop through experiments.
# For each experiment, mMake a suitable config file and set up a subdirectory.

for expt in experiments:

    # Make a copy of the calvingMIPInit config file.
    # Below, this copy will be tailored for the chosen experiment.
    # Note: Experiment 4 needs two config files: one for the retreat phase (years 0 to 500)
    #       and one for the advance phase (years 500 to 1000).
    #       We make the retreat file first, then the advance file later.

    newConfigFile = 'calvingMIP.' + expt + '.config'
#    print('Config file for this experiment:', newConfigFile)
    shutil.copy(masterConfigFile, newConfigFile)

    # Read the new config file.
    config = ConfigParser()
    config.read(newConfigFile)

    # Experiment-specific settings
    # Note: The standard experiments are SpinupCircular, SpinupThule, and Experiments 1 through 5.
    #       The Spinup runs end at the year given by 'yearsSpinup'.
    #       The other experiments start from the end of the Spinup and finish at year 1000.

    # TODO: Set up these experiments as hybrid restart (restart = 2), starting at year 0
    #       by reading the last time slice from the Spinup restart file.
    #       First need to merge the damage/calving branch to main.
    
    if expt == 'SpinupCircular':
        tstart      = 0.0
        tend        = yearsSpinup
        inputdir    = '../'
        inputfile   = inputfileCircular
        inputslice  = 1
        outputfreq  = args.outputfreqlo
        restartfreq = min(restartfreqSpinup, yearsSpinup)
    elif expt == 'SpinupThule':
        tstart      = 0.0
        tend        = yearsSpinup
        inputdir    = '../'
        inputfile   = inputfileThule
        inputslice  = 1
        outputfreq  = args.outputfreqlo
        restartfreq = min(restartfreqSpinup, yearsSpinup)
    elif expt == 'Experiment1':
        tstart      = 0.0
        tend        = 1000.0
        inputdir    = '../SpinupCircular/'
        inputfile   = 'calvingMIP.SpinupCircular.restart.nc'
        inputslice  = int(yearsSpinup/restartfreqSpinup)
        outputfreq  = args.outputfreqmd
        restartfreq = 1000.
    elif expt == 'Experiment2':
        tstart      = 0.0
        tend        = 1000.0
        inputdir    = '../SpinupCircular/'
        inputfile   = 'calvingMIP.SpinupCircular.restart.nc'
        inputslice  = int(yearsSpinup/restartfreqSpinup)
        outputfreq  = args.outputfreqmd
        restartfreq = 1000.
    elif expt == 'Experiment3':
        tstart      = 0.0
        tend        = 1000.0
        inputdir    = '../SpinupThule/'
        inputfile   = 'calvingMIP.SpinupThule.restart.nc'
        inputslice  = int(yearsSpinup/restartfreqSpinup)
        outputfreq  = args.outputfreqmd
        restartfreq = 1000.
    elif expt == 'Experiment4':
        # parameters for the retreat phase (first 500 years)
        tstart      = 0.0
        tend        = 500.0
        inputdir    = '../SpinupThule/'
        inputfile   = 'calvingMIP.SpinupThule.restart.nc'
        inputslice  = int(yearsSpinup/restartfreqSpinup)
        outputfreq  = args.outputfreqmd
        restartfreq = 500.
    elif expt == 'Experiment5':
        tstart      = 0.0
        tend        = 5000.0
        inputdir    = '../SpinupThule/'
        inputfile   = 'calvingMIP.SpinupThule.restart.nc'
        inputslice  = int(yearsSpinup/restartfreqSpinup)
        outputfreq  = args.outputfreqlo
        restartfreq = 1000.

    # Set the prescribed advance/retreat rate for Experiments 2 and 4
    if expt == 'Experiment2':
        cf_advance_retreat_amplitude = -300.
        cf_advance_retreat_period = 1000.
    elif expt == 'Experiment4':
        # for retreat500
        cf_advance_retreat_amplitude = -750.
        cf_advance_retreat_period = 1000.
        # for advance_500_1000
        cf_advance_retreat_amplitude = 5000.
        cf_advance_retreat_period = 0.
        
    #TODO: Modify for Experiment 5
            
    # Set other parameters specific to certain experiments
    # TODO: Do we need to read in the input temperature? Or do we always want temp_init = 1?

    if expt == 'SpinupCircular' or expt == 'SpinupThule':
        config.set('options', 'temp_init', '1')
        config.set('options', 'marine_margin', '5')
        config.set('ho_options', 'which_ho_calving_front', '0')
    elif expt == 'Experiment2' or expt == 'Experiment4':
        config.set('parameters', 'cf_advance_retreat_amplitude', str(cf_advance_retreat_amplitude))
        config.set('parameters', 'cf_advance_retreat_period', str(cf_advance_retreat_period))
    elif expt == 'Experiment5':
        config.set('options', 'marine_margin', '7')
        calvingMinthck = 325.
        config.set('parameters', 'calving_minthck', str(calvingMinthck))
            
    # Set the start and end times
    config.set('time', 'tstart', str(tstart))
    config.set('time', 'tend',   str(tend))

    # Set the hybrid restart option for Experiments 1 to 5
    if expt != 'SpinupCircular' and expt != 'SpinupThule':
        config.set('options', 'restart', '2')

    # Change the default comment
    comment = 'CalvingMIP ' + expt
    config.set('CF default', 'comment', comment)

    # Set input file and time slice in the section '[CF input]'.
    # Note: This method may not be robust for Spinup runs that start and restart.
    #       For this reason, the script calvingMIPRun.py makes sure the 'time' entry
    #       in [CF input] corresponds to the final time slice.
    config.set('CF input', 'name', inputfile)
    config.set('CF input', 'time', str(inputslice))
#    print('Input file:', inputfile)

    # Set the output files
    # Note: The first output file is already present in calvingMIP.config.template;
    #        here we just rewrite the name according to the experiment.
    #        The output frequency (typically 100 yr) is set in the template file.
    #       The other output files are added here with variable lists tailored to calvingMIP,
    #        with an output frequency that can be passed to this python script.

    # Set the output filename in the section '[CF output]'.
    outputfile = 'calvingMIP.' + expt + '.out.nc'
    config.set('CF output', 'name', outputfile)
    config.set('CF output', 'freq', str(outputfreq))
    #    print('Output file:', outputfile)

    # Specify additional output files for CalvingMIP experiments.
    # Use the high output frequency for scalars, but use the frequency set above for 2D fields.
    # Note: Each section must have a unique name, since configparser does not allow duplicate section names.
    #       CISM handles this by treating sections called 'CF output*' the same as 'CF output',
    #        where * is any string.
    
    if expt != 'SpinupCircular' and expt != 'SpinupThule':

        # fields on the x1 grid
        config.add_section('CF output1')
        outputfile = 'calvingMIP.' + expt + '.out.x1.nc'
        config.set('CF output1', 'name',      outputfile)
        config.set('CF output1', 'frequency', str(outputfreq))
        config.set('CF output1', 'variables',
                   'thk topg usurf lsurf calving_flux_tavg ice_mask grounded_mask floating_mask ice_mask_stag f_ground_cell effective_areafrac')
#        print('Output file:', outputfile)

        # fields on the x0 grid
        config.add_section('CF output2')
        outputfile = 'calvingMIP.' + expt + '.out.x0.nc'
        config.set('CF output2', 'name',      outputfile)
        config.set('CF output2', 'frequency', str(outputfreq))
        config.set('CF output2', 'variables', 'uvel_mean vvel_mean btract f_ground')
#        print('Output file:', outputfile)

        # scalars
        config.add_section('CF output3')
        outputfile = 'calvingMIP.' + expt + '.out.scalars.nc'
        config.set('CF output3', 'name',      outputfile)
        config.set('CF output3', 'frequency', str(args.outputfreqhi))
        config.set('CF output3', 'variables', 'iareaf iareag imass imass_above_flotation total_calving_flux total_gl_flux')
#        print('Output file:', outputfile)


    # Set restart info in the section '[CF restart]'.
    # Note: Each experiment (except the Spinup experiments) writes only one time slice to a restart file.
    restartfile = 'calvingMIP.' + expt + '.restart.nc'
    config.set('CF restart', 'name',       restartfile)
    config.set('CF restart', 'variables', 'restart')
    config.set('CF restart', 'frequency',  str(restartfreq))
    config.set('CF restart', 'xtype',     'double')
    config.set('CF restart', 'write_init', 'False')
#    print('Restart file:', restartfile)


    # Write to the new config file.
    with open(newConfigFile, 'w') as configfile:
        config.write(configfile)


    # Create a subdirectory named for the experiment, and stage the run there.
    try:
        os.mkdir(expt)
        print('Created subdirectory', expt)
    except:
        print('Subdirectory', expt, 'already exists')

    os.chdir(expt)

    # Move the config file from the parent directory to the subdirectory.
    shutil.move('../' + newConfigFile, newConfigFile)
    print('Created config file', newConfigFile)


    # For Experiment 4, create a second config file for the advance phase (years 500 to 1000).
    if (expt == 'Experiment4'):
        
        # Rename the current file to indicate the retreat phase
        ConfigNameRetreat = 'calvingMIP.' + expt + '.retreat_500.config'
        os.rename(newConfigFile,ConfigNameRetreat)

        # Make a copy; this will be the advance file                                                                                                       
        newConfigFile = 'calvingMIP.' + expt + '.advance_1000.config'
        shutil.copy(ConfigNameRetreat,newConfigFile)

        # Read the new config file
        config.read(newConfigFile)

        # Set some config values that differ from the retreat file
        tend = 1000.0
        cf_advance_retreat_amplitude = 5000.
        cf_advance_retreat_period = 0.

        config.set('time', 'tend',   str(tend))
        config.set('options', 'apply_calving_mask', 'True')
        config.set('options', 'restart', '1')
        config.set('parameters', 'cf_advance_retreat_amplitude', str(cf_advance_retreat_amplitude))
        config.set('parameters', 'cf_advance_retreat_period', str(cf_advance_retreat_period))

        # Write to the new config file.
        with open(newConfigFile, 'w') as configfile:
            config.write(configfile)

    # Link to the cism_driver executable in the parent directory.
    try:
        os.symlink('../cism_driver', 'cism_driver')
    except:
        pass   # link to cism_driver already exists

    # Link to the derecho-intel-modules file in the parent directory.
    #TODO: Make this machine-dependent?
    try:
        os.symlink('../derecho-intel-modules', 'derecho-intel-modules')
    except:
        pass   # link to derecho_modules already exists

    # Link to the input file in the appropriate directory.
    try:
        os.symlink(inputdir + inputfile, inputfile)
    except:
        pass   # link to the input file already exists

    # Go back to the parent directory and continue.
    os.chdir('..')
