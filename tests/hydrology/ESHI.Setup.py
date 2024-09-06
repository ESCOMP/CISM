#!/usr/bin/env python

# This script sets up Suite A SHMIP experiments for the Efficient Subglacial Hydrology for Ice Sheets (ESHI).
# More information about ESHI (will) be found on github:
# https://github.com/link.
# For now, the script sets up Experiments 1 to X of ESHI Phase 1.

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

# the base case topography for ESHI experiments is 100 km long and 20 km wide
# later: thk[:,:,i] = 6*(np.sqrt(x+5000)- np.sqrt(5000)+1) 

max_x = 100.e3 # (m) length of ice-sheet
max_y = 20.e3 # (m) width of ice-sheet

# Physical constants - ESHI isn't linked with Ice Dynamics, so some constants won't be used
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
# accum = 0.3     # Surface mass balance (m a-1)
rhoi = 917.     # Ice density (kg m-3)
rhoo = 1028.    # Sea water density (kg m-3)
# A = 2.9377e-18  # Ice rate factor (Pa-3 a-1)
# n = 3           # Flow law stress exponent
# C = 10000.      # Basal slipperiness (m a-1 Pa-3) --- well this is what might be informed by ESHI, we aren't evolving the ice-dynamics yet
# m = 3.          # Sliding law stress exponent
sec2yr = 31556926  # Seconds in a year (s)


####################################
# Function used later in the code #
####################################


# A future bed topography function in which the bed has a linear slope
def computeBedLinear(x,Bc,Bl, alpha):
    # Input:
    # x = grid x-coordinate in km
    # Bl = minimum bed topography elevation at the edge of the domain (m)
    # alpha = the slope coefficient (m/km)
    #
    # Output:
    # B = bed elevation (m)
    B = (alpha*x) + Bl


#############
# Main Code #
#############

# Parse options.
parser.add_argument('-c', '--configfile',   type=str,   default='ESHI.config.template', help="config file template")
parser.add_argument('-e', '--executable',   type=str,   default='cism_driver', help="path to the CISM executable")
parser.add_argument('-x', '--experiment',   type=str,   default= 'all',   help="ESHI experiment(s) to set up")
parser.add_argument('-r', '--resolution',   type=float, default= 5000.,   help="grid resolution (m)")
parser.add_argument('-h', '--hydromodel',   type=str,   default= 'bwat_constant',  help="Hydrology model (bwat_constant, local_till, ss_flux, mp_sheet, cav_sheet)")
# parser.add_argument('-t', '--timestep',     type=float, default= 1.0,     help="time step (yr)")
# parser.add_argument('-v', '--vertlevels',   type=int,   default= 5,       help="no. of vertical levels")
# parser.add_argument('-a', '--approximation',type=str,   default= 'DIVA',  help="Stokes approximation (SSA, DIVA, BP)")
# parser.add_argument('-fl', '--outputfreq',type=float,   default= 1.0,   help="low output frequency (yr)")

args = parser.parse_args()

if args.experiment == 'all':
    experiments = ['SteadyStateA','SteadyStateB','SteadyStateC','SteadyStateD','SteadyStateE','SteadyStateF'] # future experiment, steadystate w/ linear bed slope
    print('Setting up all the ESHI experiments')
elif args.experiment in ['SteadyStateA','SteadyStateB','SteadyStateC','SteadyStateD','SteadyStateE','SteadyStateF']: 
    experiments = [args.experiment]
    print( 'Setting up experiment', args.experiment)
else:
    sys.exit('Please specify experiment(s) from this list: all, SteadyState') 

# If the executable is not called 'cism_driver', then create a link called 'cism_driver' to the executable.
# Each subdirectory will link to cism_driver in the main directory.
if args.executable != 'cism_driver':
    # Remove the existing link, if present.
    os.unlink('cism_driver')
    # Make the new link.
    os.symlink(args.executable, 'cism_driver')

# Set the basal hydrology model - bwat_constant, local_till, ss_flux, mp_sheet, cav_sheet
# if args.hydromodel == 'bwat_constant':
#     # constant basal water
#     config.set('basal_hydro', 'which_ho_bwat', '0')
#     config.set('basal_hydro', 'const_bwat', '10.0')

# elif args.hydromodel == 'local_till':
#     # local till
#     config.set('basal_hydro', 'which_ho_bwat', '1')
#     config.set('basal_hydro', 'c_drainage', '1.0e-3')
#     # config.set('basal_hydro', 'N_0', '1000')
#     # config.set('basal_hydro', 'e_0', '0.69')
#     # config.set('basal_hydro', 'C_c', '0.12')

# elif args.hydromodel == 'ss_flux':
#     # steady state flux
#     config.set('basal_hydro', 'which_ho_bwat', '2')
#     config.set('basal_hydro', 'const_source', '0.0')

if args.hydromodel == 'mp_sheet':
    # macroporous sheet (new!)
    config.set('basal_hydro', 'which_ho_bwat', '3')
    config.set('basal_hydro', 'which_ho_effecpress', '2')
    config.set('basal_hydro', 'bwat_threshold', '0.1')
    config.set('basal_hydro', 'bwat_gamma', '3.5')

elif args.hydromodel == 'cav_sheet':
    # cavity sheet (new!)
    config.set('basal_hydro', 'which_ho_bwat', '3')
    config.set('basal_hydro', 'which_ho_effecpress', '3')
    config.set('basal_hydro', 'bump_height', '0.1')
    config.set('basal_hydro', 'bump_wavelength', '2.0')
    config.set('basal_hydro', 'sliding_speed_fixed', '31.5')
    config.set('basal_hydro', 'flwa_basal', '1.064e-16')
    config.set('basal_hydro', 'c_close', '0.074')

    # what about..
    # config.set('basal_hydro', 'effecpress_delta','0.02')
    # config.set('basal_hydro', 'bpmp_threshold','0.1')

else:
    sys.exit('Please specify a hydrology model from this list: bwat_constant, local_till, ss_flux, mp_sheet, cav_sheet')

# Set grid resolution.
if max_x%args.resolution==0 and max_y%args.resolution==0:
    dx = args.resolution
    dy = args.resolution
else:
    sys.exit('Your resolution should be a divider of the domain size of 100 km x 20 km')

if args.vertlevels >= 2:
    nz = args.vertlevels
else:
    sys.exit('Error: must have at least 2 vertical levels')

print( 'ESHI grid resolution (m) =', args.resolution)
print( 'Number of vertical levels =', nz)

# Set number of grid cells in each direction.
nx = int(max_x/dx) # 100 km length
ny = int(max_y/dy) # 20 km width

print('grid dimension in x:',nx)
print('grid dimension in y:',ny)

# Copy the config template to a new master config file.
# Later, the master file will be tailored to each experiment.
masterConfigFile = 'ESHI.config.template'
# masterConfigFile = 'ESHI.config'

try:
    shutil.copy(args.configfile, masterConfigFile)
except FileNotFoundError:
    sys.exit(f'File not found: {args.configfile}')
except PermissionError:
    sys.exit(f'Permission denied: {args.configfile}')
except Exception as e:
    sys.exit(f'Could not copy {args.configfile} to {masterConfigFile}: {e}')


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
# By default, this cell lies at the center of the x-axis
idiag = int(nx/2) + 1
jdiag = int(ny/2) + 1
config.set('time', 'idiag', str(idiag))
config.set('time', 'jdiag', str(jdiag))
print('idiag:',idiag)
print('jdiag:',jdiag)

# # Set the Stokes approximation - leaving this for now, but we aren't driving ice dynamics
# if args.approximation == 'SSA':
#     which_ho_approx = 1
#     print( 'Using SSA velocity solver')
# elif args.approximation == 'DIVA':
#     which_ho_approx = 4
#     print( 'Using DIVA velocity solver')
# elif args.approximation == 'BP':
#     which_ho_approx = 2
#     print( 'Using Blatter-Pattyn velocity solver')
# else:
#     which_ho_approx = 4
#     print( 'Defaulting to DIVA velocity solver')

# config.set('ho_options', 'which_ho_approx', str(which_ho_approx))

# Write to the master config file.
with open(masterConfigFile, 'w') as configfile:
    config.write(configfile)


# Create the netCDF input files: ESHI.input.nc

try:
    parser = ConfigParser()
    parser.read(args.configfile)
    inputfile = parser.get('CF input', 'name')
except:
    sys.exit('Error parsing ' + args.configfile)

print('Creating input file', inputfile)
ncfile = Dataset(inputfile, 'w')

# Create dimensions.
# Note: (x0,y0) = staggered (velocity) grid
#       (x1,y1) = unstaggered (scalar) grid
ncfile.createDimension('time', 1)
ncfile.createDimension('x1',  nx)
ncfile.createDimension('y1',  ny)
# ncfile.createDimension('x0',  nx-1)
# ncfile.createDimension('y0',  ny-1)
ncfile.createDimension('level',         nz)
# ncfile.createDimension('staglevel',     nz-1)
ncfile.createDimension('stagwbndlevel', nz+1)

# Create time and grid variables.
# Note: (x1,y1) are loadable and need to be in the input file.
#       (x0,y0) are not loadable, but are derived in CISM from (x1,y1).

ncfile.createVariable('time','f4',('time',))[:] = [0]
x1 = ncfile.createVariable('x1','f4',('x1',))
y1 = ncfile.createVariable('y1','f4',('y1',))
# x0 = ncfile.createVariable('x0','f4',('x0',))
# y0 = ncfile.createVariable('y0','f4',('y0',))

# Create 2D input fields
thk  = ncfile.createVariable('thk',  'f4', ('time','y1','x1'))
topg = ncfile.createVariable('topg', 'f4', ('time','y1','x1'))
# acab = ncfile.createVariable('acab', 'f4', ('time','y1','x1'))

# Compute x and y on each grid.
# The origin (x = y = 0) is placed at the terminus of the domain.

x = np.arange(0,max_x,dx,dtype='float32')
y = np.arange(0,max_y,dy,dtype='float32')

x1[:] = x[:]
y1[:] = y[:]

# x0[:] = x[:-1] + dx/2.   # x = -R+dx/2,..., -dx/2, dx/2,..., R-dx/2
# y0[:] = y[:-1] + dy/2.   # x = -R+dy/2,..., -dy/2, dy/2,..., R-dy/2

# Set SMB
# acab[:,:,:] = accum

# Set initial thickness
thk[:,:,:] = 0.
for i in range(0,nx):
    x = i*dx
    thk[:,:,i] = 6*(np.sqrt(x+5000)- np.sqrt(5000)+1) # this is the ice sheet elevation from SHMIP Suite A


# Set bed topography
topg[:,:,:] = 0. # simple case

# another bed topo case can be to vary the slope

# Close the file 
ncfile.close()

print('Experiments:', experiments)

# Loop through experiments.
# For each experiment, Make a suitable config file and set up a subdirectory.

for expt in experiments:
    newConfigFile = 'ESHI.' + expt + '.config'
#    print('Config file for this experiment:', newConfigFile)
    shutil.copy(masterConfigFile, newConfigFile)

    # Read the new config file.
    config = ConfigParser()
    config.read(newConfigFile)

    if expt == 'SteadyStateA':
        tstart      = 0.0
        tend        = args.timestep # run for one timestep
        inputdir    = '../'
        inputfile   = inputfile
        inputslice  = 1
        outputfreq  = args.outputfreq
        m           = 7.93e-11 # water input, m s-1
    elif expt == 'SteadyStateB':
        tstart      = 0.0
        tend        = args.timestep # run for one timestep
        inputdir    = '../'
        inputfile   = inputfile
        inputslice  = 1
        outputfreq  = args.outputfreq
        m           = 1.59e-9 # water input, m s-1
    elif expt == 'SteadyStateC':
        tstart      = 0.0
        tend        = args.timestep # run for one timestep
        inputdir    = '../'
        inputfile   = inputfile
        inputslice  = 1
        outputfreq  = args.outputfreq
        m           = 5.79e-9 # water input, m s-1
    elif expt == 'SteadyStateD':
        tstart      = 0.0
        tend        = args.timestep # run for one timestep
        inputdir    = '../'
        inputfile   = inputfile
        inputslice  = 1
        outputfreq  = args.outputfreq
        m           = 2.5e-8 # water input, m s-1
    elif expt == 'SteadyStateE':
        tstart      = 0.0
        tend        = args.timestep # run for one timestep
        inputdir    = '../'
        inputfile   = inputfile
        inputslice  = 1
        outputfreq  = args.outputfreq
        m           = 4.5e-9 # water input, m s-1
    elif expt == 'SteadyStateF':
        tstart      = 0.0
        tend        = args.timestep # run for one timestep
        inputdir    = '../'
        inputfile   = inputfile
        inputslice  = 1
        outputfreq  = args.outputfreq
        m           = 5.79e-7 # water input, m s-1

        
    # Set the start and end times
    config.set('time', 'tstart', str(tstart))
    config.set('time', 'tend',   str(tend))

    # Change the default comment
    comment = 'ESHI ' + expt
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
    outputfile = 'ESHI.' + expt + '.out.nc'
    config.set('CF output', 'name', outputfile)
    config.set('CF output', 'freq', str(outputfreq))
    #    print('Output file:', outputfile)

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