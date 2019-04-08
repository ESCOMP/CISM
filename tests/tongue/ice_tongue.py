#!/usr/bin/env python
# This script runs a one-dimensional ice tongue test case used for verifying damage.
# Pre-loaded ice shelves: Erebus & Drygalski (unconfined); Amery & Ross (confined)
# Last edited by Morgan Whitcomb at the University of Michigan on April 8, 2019.

# Imports
from optparse import OptionParser
from ConfigParser import ConfigParser
from netCDF import *
import numpy as np
import sys,os,datetime

# User-defined inputs ##
const_damage = 0.1   # Constant initial damage value
#amery_melt = 0.6   # Melt rate across the Amery ice shelf (m/a)
amery_melt = 0.8    #
#amery_melt = 1.    #
#ross_melt = 0.2   # Melt rate across the Ross ice shelf (m/a)
ross_melt = 0.25   #
#ross_melt = 0.5   #
##

# Classes #####
class constants:
   # Structure for holding general constants
   def __init__(s):
      s.n = 3.   # Glen's Flow Law exponent
      s.rhoi = 910.   # Ice density (kg/m^3)
      s.rhoo = 1028.   # Seawater density (kg/m^3)
      s.grav = 9.81   # Gravitational acceleration (m/s^2)

class geometry_vals:
   # Structure for holding parameters specific to the geometry of interest
   def __init__(s,options,c,nx,dx,ny,dy,nz,flwa,dinit,amery_melt,ross_melt):
      # Initialize resolution values
      s.nx = nx
      s.nxst = nx-1
      s.dx = dx
      s.ny = ny
      s.nyst = ny-1
      s.dy = dy
      s.nz = nz
      s.nzst = nz-1

      # Compute vectors and matrices of spatial grid points
      s.x = s.dx*np.arange(s.nx,dtype='float32')
      s.xst = s.x[:-1] + s.dx/2.
      s.y = s.dy*np.arange(s.ny,dtype='float32')
      s.yst = s.y[:-1] + s.dy/2.
      s.xgrid,s.ygrid = np.meshgrid(s.x,s.y,indexing='xy')
      s.xstgrid,s.ystgrid = np.meshgrid(s.xst,s.yst,indexing='xy')

      # Initialize geometry-specific parameters
      if options.shelf == 'erebus':
         s.set_shelf_params(340.,100.,2.2)
      elif options.shelf == 'drygalski':
         s.set_shelf_params(575.,940.,6.8)
      elif options.shelf == 'amery':
         s.set_shelf_params(1090.,390.,amery_melt,conflen=505.)
      elif options.shelf == 'ross':
         s.set_shelf_params(660.,340.,ross_melt,conflen=650.)
      s.flwa = flwa   # Ice viscosity (Pa^[-n] yr^[-1])
      s.C = s.flwa*(c.rhoi*c.grav*(c.rhoo-c.rhoi)/(4.*c.rhoo))**c.n

      # Set the initial damage field
      s.dinit = dinit

      # Specify matrix indices
      s.ptnum = int((s.nx+2)/3)   # Number of points along the initial ice tongue
      s.gl = 2   # Grounding line location
      s.cf = s.ptnum+2   # Calving front location

   def set_shelf_params(s,hgl,ugl,mdot,conflen=0.):
      # Load the ice shelf grounding line flux and melt rate into the struct
      s.hgl = hgl   # Grounding line thickness (m)
      s.ugl = ugl   # Grounding line velocity (m/a)
      s.mdot = mdot   # Melt rate (m/a)
      s.conflen = conflen   # Downstream extent of lateral embayment walls (km)
      s.confend = int(s.conflen*1000./s.dx)   # Matrix index corresponding to conflen

   def construct_ice_tongue(s,options,c):
      # Compute the necessary fluid fields and establish boundary conditions
      # Initialize matrices
      s.thk = np.zeros([1,s.ny,s.nx],dtype='float32')   # Thickness (m)
      s.acab = np.zeros([1,s.ny,s.nx],dtype='float32')   # Mass balance (m/a)
      s.topg = np.zeros([1,s.ny,s.nx],dtype='float32')   # Bed elevation (m)
      s.beta = np.zeros([1,s.nyst,s.nxst],dtype='float32')   # Higher-order bed stress
      s.kbc = np.zeros([1,s.nyst,s.nxst],dtype='int')   # Kinetic boundary condition mask
      s.no_advance_mask = np.zeros([1,s.ny,s.nx],dtype='int')   # No advance mask
      s.const_thk_mask = np.zeros([1,s.ny,s.nx],dtype='int')   # Constant thickness mask
      # We need to use the extended staggered mesh for our velocities to avoid boundary
      # condition problems - the fact that the standard staggered mesh is smaller than the
      # unstaggered mesh by 1 point per dimension causes numerical issues when using
      # periodic global boundaries with nonzero velocities at the walls.
      s.uvel_extend = np.zeros([1,s.nz,s.ny,s.nx],dtype='float32')   # X-velocity (m/a)

      # Set the climate forcing
      s.acab[0,:,s.gl+1:] = -s.mdot   # Mass balance is applied everywhere but the gl

      # Prevent the ice from flowing backwards from the gl and from advancing beyond the cf
      s.no_advance_mask[0,:,:s.gl] = 1   # Backwards from gl
      if options.length == 'const':
         s.no_advance_mask[0,:,s.cf+1:] = 1   # Advancing from cf

      # Set the kinetic boundary conditions to hold the gl velocity constant in time
      s.kbc[0,:,:s.gl+1] = 1

      # Hold the gl thickness constant in time
      s.const_thk_mask[0,:,s.gl] = 1

      # Compute the ice thickness
      s.thk[0,:,s.gl] = s.hgl
      for j in range(0,s.ny):
         s.thk[0,j,s.gl+1:s.cf+1] = s.compute_thk(c,s.x[s.gl+1:s.cf+1],s.x[s.gl])

      # Compute the ice velocity
      s.compute_uvel(c)

      # Lower the bedrock elevation everywhere to make the ice shelf float everywhere but the gl
      # We have to raise topg by a small amount in order to account for floating point inaccuracies,
      # which we've chosen here to use half the difference between the grounding line thickness and the
      # thickness in the grid cell immediately downstream from it
      for i in range(0,s.nx):
         s.topg[0,:,i] = -(c.rhoi/c.rhoo)*s.hgl + 0.5*(s.hgl-s.thk[0,:,s.gl+1])

      # Enforce lateral confinement
      if s.conflen > 0.:   # If laterally confined
         s.enforce_confinement(s.acab[0],0.,'acab',s.confend)
         s.enforce_confinement(s.thk[0],0.,'thk',min(s.cf,s.confend))
         s.enforce_confinement(s.uvel_extend[0],0.,'uvel',min(s.cf,s.confend))
         s.enforce_confinement(s.topg[0],s.hgl,'topg',s.confend)
         s.enforce_confinement(s.kbc[0],1,'kbc',s.confend)

   def compute_thk(s,c,x,xgl):
      # Compute the analytic ice thickness, via van der Veen Ch. 5
      if s.mdot == 0.:   # Accumulation = ablation
         thk = ((c.n+1.)*s.C*(x-xgl)/(s.hgl*s.ugl)+s.hgl**(-(c.n+1.)))**(-1./(c.n+1.))

      else:   # Nonzero accumulation
         thk = (s.ugl**(c.n+1.)*(s.C*s.hgl**(c.n+1.)/s.mdot+1.) \
                /(s.hgl*s.ugl-s.mdot*(x-xgl))**(c.n+1.)-s.C/s.mdot)**(-1./(c.n+1.))

      return thk

   def compute_uvel(s,c):
      # Compute the analytic ice velocity from conservation of flux
      # Calculate the thickness beyond the calving front for smooth interpolation
      thktemp = s.thk[0,0,s.gl:s.cf+2].copy()
      thktemp[-1] = s.compute_thk(c,s.x[s.cf+1],s.x[s.gl])

      # Linearly interpolate the thickness onto the staggered grid
      thk_stag = np.zeros(len(thktemp[:-1]))
      thk_stag = 0.5*(thktemp[:-1]+thktemp[1:])
      
      # Compute the velocity
      for i in range(s.gl,s.cf+1):
         s.uvel_extend[0,:,:,i] = (s.hgl*s.ugl-s.mdot*(s.xst[i]-s.x[s.gl]))/thk_stag[i-s.gl]

   def enforce_confinement(s,field,val,fieldstr,xend):
      # Enforce lateral confinement in the specified fluid field
      twod_fields = ['acab','thk','topg','kbc']
      threed_fields = ['uvel','damage']
      if fieldstr in twod_fields:
         field[0,:xend+1] = val
         field[-1,:xend+1] = val
      elif fieldstr in threed_fields:
         field[:,0,:xend+1] = val
         field[:,-1,:xend+1] = val

         if fieldstr == 'uvel':   # Need to include an extra row for uvel_extend
            field[:,-2,:xend+1] = val

      # Save the altered field to the struct
      if fieldstr == 'acab':
         s.acab[0] = field
      elif fieldstr == 'thk':
         s.thk[0] = field
      elif fieldstr == 'uvel':
         s.uvel_extend[0] = field
      elif fieldstr == 'topg':
         s.topg[0] = field
      elif fieldstr == 'damage':
         s.damage[0] = field

class restart_file:
   # Structure for holding data from the restart file, when specified
   def __init__(s,g,options):
      message('Copying restart data from "'+options.restfile+'"')
      
      # Open the restart file
      infile = open_netcdf_file(options.restfile,'r')

      # Extract the grid and time information from the restart file
      s.x = np.array(infile.variables['x1'][:])
      s.xst = np.array(infile.variables['x0'][:])   # Staggered grid
      s.y = np.array(infile.variables['y1'][:])
      s.yst = np.array(infile.variables['y0'][:])   # Staggered grid
      s.z = np.array(infile.variables['level'][:])
      s.zst = np.array(infile.variables['staglevel'][:])   # Staggered grid
      s.t = np.array(infile.variables['time'][:])

      # Parse the spatial and temporal vectors to find the resolution
      s.nx,s.dx = s.find_res(s.x)
      s.nxst = len(s.xst)
      s.ny,s.dy = s.find_res(s.y)
      s.nyst = len(s.yst)
      s.nz = len(s.z)
      s.nzst = len(s.zst)
      s.nt,s.dt = s.find_res(s.t)

      # Error checking
      if s.nx != g.nx or s.dx != g.dx or s.ny != g.ny or s.dy != g.dy or s.nz != g.nz:
         # If the resolutions don't match, quit
         error('the resolution specified in the configuration file must match that of the restart file')
         
      # Extract the necessary variables from the restart file
      s.thk = np.array(infile.variables['thk'][-1,:,:])   # Ice thickness (m)
      s.uvel = np.array(infile.variables['uvel'][-1,:,:,:])   # X-velocity (m/a)
      s.vvel = np.array(infile.variables['vvel'][-1,:,:,:])   # Y-velocity (m/a)
      s.kbc = np.array(infile.variables['kinbcmask'][-1,:,:])   # Kinetic boundary condition mask
      s.const_thk_mask = np.array(infile.variables['const_thk_mask'][-1,:,:])  # Constant thickness mask
      s.no_advance_mask = np.array(infile.variables['no_advance_mask'][-1,:,:])
      s.efvs = np.array(infile.variables['efvs'][-1,:,:,:])   # Ice viscosity
      s.tau_xx = np.array(infile.variables['tau_xx'][-1,:,:,:])   # Deviatoric X-stress (Pa)
      s.tau_yy = np.array(infile.variables['tau_yy'][-1,:,:,:])   # Deviatoric Y-stress (Pa)
      s.eps_xx = np.array(infile.variables['eps_xx'][-1,:,:,:])   # X-strain rate (1/a)
      s.eps_yy = np.array(infile.variables['eps_yy'][-1,:,:,:])   # Y-strain rate (1/a)
      s.topg = np.array(infile.variables['topg'][-1,:,:])   # Bed topography (m)
      s.acab = np.array(infile.variables['acab'][-1,:,:])   # Mass balance (m/a)

      # Create the higher-order bed stress - zeros everywhere for floating ice
      s.beta = np.zeros([1,s.nyst,s.nxst],dtype='float32')   # Higher-order bed stress

      # Create extended versions of the velocity matrices by directly copying from nearest-neighbor
      # interpolation
      s.uvel_extend = s.create_extended_vel(s.uvel)
      s.vvel_extend = s.create_extended_vel(s.vvel)

   def find_res(s,i):
      # Pull the number of points and resolution from a vector
      ni = len(i)
      di = i[1]-i[0]

      return ni,di

   def create_extended_vel(s,vel):
      # Create an extended version of a velocity matrix
      vel_extend = np.zeros([1,s.nz,s.ny,s.nx],dtype='float32')

      # Copy the non-extended velocity into the extended matrix
      vel_extend[0,:,:-1,:-1] = vel

      # Fill the remaining values by nearest-neighbor interpolation
      vel_extend[0,:,-1,:-1] = vel[:,-1,:]
      vel_extend[0,:,:-1,-1] = vel[:,:,-1]
      vel_extend[0,:,-1,-1] = vel[:,-1,-1]

      return vel_extend
#####

# Functions #####
def cmdline_opts(cfgfile='ice_tongue.config'):
   # Parse command-line options
   optparser = OptionParser()
   optparser.add_option('-c','--config',dest='configfile',type='string',default=cfgfile, \
                        help='Name of .config file to use for the run', metavar='FILE')
   optparser.add_option('-m','--parallel',dest='parallel',type='int', \
                        help='Number of processors to run the model with; if specified, then execute'+\
                             ' the run in parallel [default: perform a serial run]', metavar="NUMPROCS")
   optparser.add_option('-e','--exec',dest='executable',type='string',default='./cism_driver', \
                        help='Set path to the CISM executable',metavar='FILE')
   optparser.add_option('-s','--shelf',dest='shelf',type='string',default='erebus',help='Ice shelf to'+\
                        ' use representative parameters for; the supported options are: erebus,'+\
                        ' drygalski, ross, amery')
   optparser.add_option('-l','--length',dest='length',type='string',default='var',help="Mode for"+\
                        " determining the final length of the ice tongue; const constrains the ice"+\
                        " tongue's maximum length, and var allows the ice tongue to advance freely")
   optparser.add_option('-r','--restart',dest='restfile',type='string', \
                        help='Specify a netCDF to restart from with a fresh damage field')
   for option in optparser.option_list:
       if option.default != ("NO", "DEFAULT"):
           option.help += (" " if option.help else "") + "[default: %default]"
   options,_ = optparser.parse_args()

   # Error checking
   if os.path.isfile(options.configfile) is False:   # Config file doesn't exist
      error('file "'+options.configfile+'" does not exist')
   if os.path.isfile(options.executable) is False:   # CISM executable doesn't exist
      error('file "'+options.executable+'" does not exist')
   # Restart file was specified but doesn't exist
   if options.restfile is not None and os.path.isfile(options.restfile) is False:
      error('file "'+options.restfile+'" does not exist')
   shelf_list = ['erebus','drygalski','ross','amery']
   if options.shelf not in shelf_list:
      error('Invalid input for shelf; accepted entries are: erebus, drygalski, ross, amery')
   length_list = ['const','var']
   if options.length not in length_list:
      error('Invalid input for length; accepted entries are: const, var')

   return options

def error(string):
   # Send an error to the user and exit the program
   print 'ERROR: '+string
   sys.exit()

def message(string):
   print string+' ... ', datetime.datetime.now().time()

def parse_config_file(options,amery_melt,ross_melt,no_dam_floor=0,nye_dam_floor=1):
   # Read resolution data from the config file
   parser = ConfigParser()
   parser.read(options.configfile)
   nx = int(parser.get('grid','ewn'))
   ny = int(parser.get('grid','nsn'))
   nz = int(parser.get('grid','upn'))
   dx = float(parser.get('grid','dew'))
   dy = float(parser.get('grid','dns'))
   filename = str(parser.get('CF input','name'))
   flwa = float(parser.get('parameters','default_flwa'))
   dam_floor = parser.get('options','damage_floor')
   dam_floor = int(dam_floor.split(' ')[0])   # Remove post-option comments
   if dam_floor == no_dam_floor:   # No damage floor (constant input value)
      dinit = const_damage 
   elif dam_floor == nye_dam_floor:   # Nye damage floor
      dinit = 0.

   # Create a geometry class instance and save the resolution data to it
   c =  constants()
   geom = geometry_vals(options,c,nx,dx,ny,dy,nz,flwa,dinit,amery_melt,ross_melt)

   return geom,c,filename

def open_netcdf_file(filename,action):
   # Open a netCDF file for the specified action
   try:
      openfile = NetCDFFile(filename,action,format='NETCDF3_CLASSIC')
   except TypeError:
      openfile = NetCDFFile(filename,action)

   return openfile

def init_netcdf(s,filename):
   # Initialize input netCDF file for CISM, using information from the configuration file
   # Create the netCDF file
   message('Writing '+filename)
   infile = open_netcdf_file(filename,'w')

   # Create netCDF dimensions for the resolution values
   infile.createDimension('time',1)
   infile.createDimension('x1',s.nx)
   infile.createDimension('x0',s.nxst)   # Staggered grid 
   infile.createDimension('y1',s.ny)
   infile.createDimension('y0',s.nyst)   # Staggered grid
   infile.createDimension('level',s.nz)
   infile.createDimension('staglevel',s.nzst)

   # Create netCDF vectors for the corresponding dimensions
   infile.createVariable('time','f',('time',))[:] = [0]
   infile.createVariable('x1','f',('x1',))[:] = s.x
   infile.createVariable('x0','f',('x0',))[:] = s.xst   # Staggered grid
   infile.createVariable('y1','f',('y1',))[:] = s.y
   infile.createVariable('y0','f',('y0',))[:] = s.yst   # Staggered grid

   return infile

def initialize_damage(s,g):
   # Initialize the damage field across the ice tongue
   # Create the damage matrix
   s.damage = np.zeros([1,s.nzst,s.ny,s.nx],dtype='float32')   # Damage

   # Set the initial damage
   s.damage[0,:,:,g.gl:g.cf+1] = g.dinit

   # Enforce lateral confinement
   if g.conflen > 0.:   # If laterally confined
      g.enforce_confinement(s.damage[0],0.,'damage',min(g.cf,g.confend))

def finalize_netcdf(s,infile,restart=False):
   # Finalize the input netCDF file for outputting to CISM
   message('Saving parameter matrices to input netCDF file')

   # Create the required variables in the netCDF file
   infile.createVariable('thk','f',('time','y1','x1'))[:] = s.thk
   infile.createVariable('uvel_extend','f',('time','level','y1','x1'))[:] = s.uvel_extend
   infile.createVariable('damage','f',('time','staglevel','y1','x1'))[:] = s.damage
   infile.createVariable('acab','f',('time','y1','x1'))[:] = s.acab
   infile.createVariable('topg','f',('time','y1','x1'))[:] = s.topg
   infile.createVariable('beta','f',('time','y0','x0'))[:] = s.beta
   infile.createVariable('kinbcmask','i',('time','y0','x0'))[:] = s.kbc
   infile.createVariable('no_advance_mask','i',('time','y1','x1'))[:] = s.no_advance_mask
   infile.createVariable('const_thk_mask','i',('time','y1','x1'))[:] = s.const_thk_mask
   if restart == True:
      infile.createVariable('vvel_extend','f',('time','level','y1','x1'))[:] = s.vvel_extend
      infile.createVariable('efvs','f',('time','staglevel','y1','x1'))[:] = s.efvs
      infile.createVariable('tau_xx','f',('time','staglevel','y1','x1'))[:] = s.tau_xx
      infile.createVariable('tau_yy','f',('time','staglevel','y1','x1'))[:] = s.tau_yy
      infile.createVariable('eps_xx','f',('time','staglevel','y1','x1'))[:] = s.eps_xx
      infile.createVariable('eps_yy','f',('time','staglevel','y1','x1'))[:] = s.eps_yy

   infile.close()

def run_cism(options,test_type='unconfined'):
   # Submit the input netCDF file for CISM to run
   print 'Running CISM for the '+test_type+' ice-tongue experiment'
   print '==============================================\n'
   if options.parallel == None:
      # Perform a serial run
      runstring = options.executable + ' ' + options.configfile
      print 'Executing serial run with:  ' + runstring + '\n\n'
      os.system(runstring)
   else:
      # Perform a parallel run
      if options.parallel <= 0:
         error('Number of processors specified for parallel run is <=0.')
      else:
         # These calls to os.system will return the exit status:
         # 0 for success (the command exists), some other integer for failure
         if os.system('which openmpirun > /dev/null') == 0:
            mpiexec = 'openmpirun -np ' + str(options.parallel)
         elif os.system('which mpirun > /dev/null') == 0:
            mpiexec = 'mpirun -np ' + str(options.parallel)
         elif os.system('which aprun > /dev/null') == 0:
            mpiexec = 'aprun -n ' + str(options.parallel)
         elif os.system('which mpirun.lsf > /dev/null') == 0:
            # mpirun.lsf does NOT need the number of processors (options.parallel)
            mpiexec = 'mpirun.lsf'
         else:
            error('Unable to execute parallel run.  Please edit the script to use your MPI run'+\
                  ' command, or run the model manually with something like: mpirun -np 4 ./cism_driver'+\
                  ' ice-tongue.config')
         runstring = mpiexec + ' ' + options.executable + ' ' + options.configfile
         print 'Executing parallel run with:  ' + runstring + '\n\n'
         os.system(runstring)  # Here is where the parallel run is actually executed!
#####

# MAIN CODE #####
options = cmdline_opts()
g,c,filename = parse_config_file(options,amery_melt,ross_melt)

# When we're not restarting from a previous file, create an input netCDF file and construct an ice
# tongue; when we are restarting, copy the input data to a new input netCDF file and create a new
# damage field.
if options.restfile is None:   # Not a restart
   infile = init_netcdf(g,filename)
   g.construct_ice_tongue(options,c)
   initialize_damage(g,g)
   finalize_netcdf(g,infile)
else:   # Restarting
   r = restart_file(g,options)
   infile = init_netcdf(r,filename)
   initialize_damage(r,g)
   finalize_netcdf(r,infile,restart=True)

run_cism(options)
#####

