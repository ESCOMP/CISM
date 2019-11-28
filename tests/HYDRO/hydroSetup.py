#!/usr/bin/env python

# This script sets up initial conditions for a test case using the Shakti
# hydrology model. For more information about the Shakti model:
#   SOMMERS, Aleah, RAJARAM, Harihar, et MORLIGHEM, Mathieu. SHAKTI: 
#   Subglacial Hydrology and Kinetic, Transient Interactions v1. 0.
#   Geoscientific Model Development, 2018, vol. 11, no 7, p. 2955-2974.


import sys, os
import shutil
import numpy as np
from netCDF4 import Dataset
from configparser import ConfigParser
from optparse import OptionParser

#print 'version of numpy =',np.version.version

#############
# Constants #
#############


xDomain = 1000.0     # domain x-dimension (m)
yDomain = 1000.0     # domain y-dimension (m)
MinutePerY = 365.*24.*60. # Number of minutes in one year 
DayPerY    = 365.         # Number of days in one year 


####################################
# Function used later in the code #
####################################


#This function set up a linear bed.
def computeBed(x,y,b0,bx):
    x  = x/1.e3     # km
    y  = y/1.e3     # km
    nx = np.size(x)
    ny = np.size(y)
    bed = np.zeros((ny,nx))
    slope = bx
    eps_b = 1e-10    # small regularization number
    abs_x = np.sqrt(x**2 + eps_b**2)  # regularizing to avoid problems at the divide

    b = b0 + slope*abs_x
    for yy in range(ny):
        bed[yy,:] = b
    return bed
    


########
# Code #
########

# Parse options.
optparser = OptionParser()

optparser.add_option('-c', '--config', dest='configfile',      type='string', default='hydro.config.template', help="config file template", metavar="FILE")
optparser.add_option('-e', '--exec',   dest='executable',      type='string', default='cism_driver', help="path to the CISM executable")
optparser.add_option('-t', '--tstep',  dest='timestep',        type='float',  default= 15.,     help="time step (min)",         metavar="TSTEP")
optparser.add_option('-y', '--tend',   dest='endtime',         type='float',  default= 30.,     help="length of simulation run (days)")
optparser.add_option('-r', '--res',    dest='resolution',      type='float',  default= 50,      help="grid resolution (m)",    metavar="RES")
optparser.add_option(      '--b0',     dest='bedheight',       type='float',  default= 100,     help="bed height at divide (m)")
optparser.add_option(      '--bx',     dest='bedslope',        type='float',  default= 2,       help="bed slope (m/km)")
optparser.add_option(      '--H0',     dest='thickinit',       type='float',  default= 500,     help="Uniform ice thickness (m)")
optparser.add_option(      '--v0',     dest='velinit',         type='float',  default= -30,     help="Constant sliding velocity (m/yr)")
 


optparser.add_option 

for option in optparser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = optparser.parse_args()


# If there is not already a link to cism_driver in the main directory, then make one.
# Each subdirectory will link to cism_driver in the main directory.
if options.executable != 'cism_driver':
    # Remove the existing link, if present.
    os.unlink('cism_driver')        
    # Make the new link.
    os.symlink(options.executable, 'cism_driver')

# initial uniform ice thikcness (m)
initThickness = options.thickinit   

vel0 = options.velinit



# Converting the time entries in years.
# The time step is given in minutes.
# Note: we will keep the rounding to 2 decimal points.
dt_inY    = options.timestep/MinutePerY        # command line input for time step in years
dec_place = -int(np.floor(np.log10(dt_inY)))   # first non zero decimal place
round_dt  = dec_place+2                        # keeping 2 digits after significant figure
dt = round(dt_inY,round_dt)                    # CISM friendly input time step

print("The time step in years is ",dt)

# Now the end time of the simulation.
# Note: We need the time step to divide the end time as an integer.
tend_inY  = options.endtime/DayPerY   # command line input for end time in years
nt        = int(tend_inY/dt)          # Number of time steps
tend      = nt*dt                     # CISM friendly end time.

print("The lenght of the simulation is now tend = ",tend)


# Set grid resolution.
dx = options.resolution
dy = dx

# Checking that the resolution divides evenly w.r.t the lenght of the domain
remx = xDomain%dx
remy = yDomain%dy
if (remx!=0):
    print("Remainder of xDomain and dx =", remx)
    print("Your resolution does not divide evenly with your x-dimension length.")
    sys.error("Modify your choice of resolution or length of domain along x dimension")
elif (remy!=0):
    print("Remainder of xDomain and dy =", remy)
    print("Your resolution does not divide evenly with your y-dimension length.")
    sys.error("Modify your choice of resolution or length of domain along y dimension")
else:
    print('Grid resolution (m) =', dx)


# Set number of grid cells in each direction.
# Include a few extra cells in the x direction to handle boundary conditions.
nx = int(xDomain/dx) + 1
ny = int(yDomain/dy)


# Copy the config template to a new master config file.
masterConfigFile = 'hydromaster.config'

try:
    shutil.copy(options.configfile, masterConfigFile)
except:
    sys.exit('Could not copy', options.configfile)

print('Creating master config file', masterConfigFile)

# Read the master config file.
config = ConfigParser()
config.read(masterConfigFile)

# Set the grid variables in the master config file.
config.set('grid', 'ewn', str(nx))
config.set('grid', 'nsn', str(ny))
config.set('grid', 'dew', str(dx))
config.set('grid', 'dns', str(dy))
nz = int(config.get('grid','upn'))

config.set('time', 'tend', str(tend))

# Write to the master config file.
with open(masterConfigFile, 'w') as configfile:
    config.write(configfile)


# Create the netCDF input file.
try:
    parser = ConfigParser()
    parser.read(options.configfile)
    initfile = parser.get('CF input', 'name')
except:
    sys.exit('Error parsing ' + options.configfile)

print('Creating input file', initfile)

# If the file exists in the path, remove it first
if os.path.exists(initfile):
    os.remove(initfile)

ncfile = Dataset(initfile, 'w')

# Create dimensions.
# Note: (x0,y0) = staggered (velocity) grid.
#       (x1,y1) = unstaggered (scalar) grid.
ncfile.createDimension('time', 1)
ncfile.createDimension('x1',  nx)
ncfile.createDimension('y1',  ny)
ncfile.createDimension('x0',  nx-1)
ncfile.createDimension('y0',  ny-1)
ncfile.createDimension('level',         nz)

# Create time and grid variables.
# Note: (x1,y1) are loadable and need to be in the input file.  
#       (x0,y0) are not loadable, but are derived in CISM from (x1,y1). May not be needed.

ncfile.createVariable('time','f4',('time',))[:] = [0]
x1 = ncfile.createVariable('x1','f4',('x1',)) 
y1 = ncfile.createVariable('y1','f4',('y1',))
x0 = ncfile.createVariable('x0','f4',('x0',))
y0 = ncfile.createVariable('y0','f4',('y0',))

# Create 2D input fields
thk  = ncfile.createVariable('thk',  'f4', ('time','y1','x1'))
topg = ncfile.createVariable('topg', 'f4', ('time','y1','x1'))
acab = ncfile.createVariable('acab', 'f4', ('time','y1','x1'))
uvel = ncfile.createVariable('uvel', 'f4', ('time','level','y0','x0'))
vvel = ncfile.createVariable('vvel', 'f4', ('time','level','y0','x0'))
kinbcmask = ncfile.createVariable('kinbcmask', 'i4', ('time','y0','x0'))
ice_hydro_mask = ncfile.createVariable('ice_hydro_mask' , 'f4', ('time','y1','x1'))
head_gradient_mask_east  = ncfile.createVariable('head_gradient_mask_east' , 'f4', ('time','y1','x1'))
head_gradient_mask_north = ncfile.createVariable('head_gradient_mask_north', 'f4', ('time','y1','x1'))


# Compute x and y on each grid.
# Note: (1) The x origin is placed at the center of the second cell from the left.
#           This assumes that kinbcmask = 1 at the first vertex from the left.
#           Thus the left edge of the grid has x = -3*dx/2.
#       (2) The y origin is placed at the bottom edge of the CISM grid.
#           The line of central symmetry runs along cell edges at y = 40 km.

x = dx * np.arange(nx,dtype='float32')   # x = 0, dx, 2*dx, etc.
y = dy * np.arange(ny,dtype='float32')   # y = 0, dy, 2*dy, etc.

x1[:] = x[:] - dx         # x1 = -dx, 0, dx, ..., (nx-2)*dx - dx/2  
y1[:] = y[:] + dy/2       # y1 = dy/2, 3*dy/2, ..., (ny-1)*dy - dy/2

x0[:] = x[:-1] - dx/2.    # x0 = -dx/2, dx/2, 3*dx/2, ..., (nx-2)*dx
y0[:] = y[:-1] + dy       # y0 = dy, 2*dy, ..., (ny-1)*dy

# Initialize thickness.
thk[:,:,:] = 0.
for i in range(nx):
    thk[:,:,i] = initThickness


# Set bed topography.
b0 = options.bedheight
bx = options.bedslope
for i in range(nx):
    for j in range(ny):
        topg[:,j,i] = computeBed(x1[i], y1[j], b0, bx)


# Set initial velocity to zero.
uvel[:,:,:,:] = vel0
vvel[:,:,:,:] = 0.

acab[:,:,:] = 0.

# Set kinematic velocity mask.
# Where kinbcmask = 1, the velocity is fixed at its initial value.
# Note: Although there is no ice on the RHS of the domain, we need kinbcmask =1 there 
#       to preserve symmetry with the LHS (since east-west BCs are formally periodic).
kinbcmask[:,:,:]  = 1   # initialize to 0 everywhere
#kinbcmask[:,:,0]  = 1   # mask out left-most column
#kinbcmask[:,:,-1] = 1   # mask out right-most column

ice_hydro_mask[:,:,:] = 1
ice_hydro_mask[:,:,0:2] = 0 
head_gradient_mask_east[:,:,:] = 0    # initialize to 0 everywhere
head_gradient_mask_north[:,:,:] = 0    # initialize to 0 everywhere
head_gradient_mask_east[:,:,-1] = 1   # mask out right-most column
head_gradient_mask_north[:,-1,:] = 1   # mask out up-most column

ncfile.close()




# Make a copy of the hydro config file.
# Below, this copy will be tailored for the chosen hydro experiment,
#  without changing the settings used for spin-up.

newConfigFile = 'hydro.config'
print('Config file for this experiment:', newConfigFile)
shutil.copy(masterConfigFile, newConfigFile)

# Read the new config file.
config = ConfigParser()
config.read(newConfigFile)

tstart      = 0.0
inputdir    = '../'
inputfile   = initfile
inputslice  = 1
outputfreq  = min(100.*dt, tend)
restartfreq = min(500.*dt, tend)

# Set the time step in the master config file.
# Set the diagnostic interval to the same value (not necessary, but helpful for debugging).
config.set('time', 'dt',      str(dt))
config.set('time', 'dt_diag', str(dt))

# Set the start and end times.
config.set('time', 'tstart', str(tstart))
config.set('time', 'tend',   str(tend))

# Change the default comment.
comment = 'Hydrology experiment '
config.set('CF default', 'comment', comment)

# Set input file and time slice in the section '[CF input]'.
print('Input file:', inputfile)
config.set('CF input', 'name', inputfile)
config.set('CF input', 'time', str(inputslice))

# Set the output filename in the section '[CF output]'.
outputfile = 'hydro.fields.out.nc'
print('Output file:', outputfile)
config.set('CF output', 'name',      outputfile)
config.set('CF output', 'frequency', str(outputfreq))

# Set restart info in the section '[CF restart]'.
restartfile = 'hydro.restart.nc'
print('Restart file:', restartfile)
config.set('CF restart', 'name',       restartfile)
config.set('CF restart', 'variables', 'restart')
config.set('CF restart', 'frequency',  str(restartfreq))
config.set('CF restart', 'xtype',     'double')
config.set('CF restart', 'write_init', 'False')

# Write to the new config file.
with open(newConfigFile, 'w') as configfile:
    config.write(configfile)

# Create a subdirectory named for the experiment, and stage the run there.
caseFolder = str(dx)+'m'
try:
    os.mkdir(caseFolder)
    print('Created subdirectory', caseFolder)
except:
    print('Subdirectory', caseFolder, 'already exists')

os.chdir(caseFolder)

# Move the config file from the parent directory to the subdirectory.
shutil.move('../' + newConfigFile, newConfigFile)
print('Created config file', newConfigFile)

# Link to the cism_driver executable in the parent directory.
try:
    os.symlink('../cism_driver', 'cism_driver')
except:
    pass   # link to cism_driver already exists

# Link to the input file in the appropriate directory.
try:
    os.symlink(inputdir + inputfile, inputfile)
except:
    pass   # link to the input file already exists

# Go back to the parent directory and continue.
os.chdir('..')
