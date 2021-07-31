#!/usr/bin/env python

"""
Run an experiment with an ice "slab". 
"""

# Authors
# -------
# Modified from dome.py by Matt Hoffman, Dec. 16, 2013
#    Test case described in sections 5.1- 5.2 of:
#    J.K. Dukowicz, 2012. Reformulating the full-Stokes ice sheet model for a
#    more efficient computational solution. The Cryosphere, 6, 21-34,
#    https://doi.org/10.5194/tc-6-21-2012.
# Reconfigured by Joseph H Kennedy at ORNL on April 27, 2015 to work with the regression testing.
#
# Revised by William Lipscomb in 2021 to support more options.
# CISM results are described in this paper:
#    Robinson, A., D. Goldberg, and W. H. Lipscomb, 2022, A comparison of the
#    stability and performance of depth-integrated ice-dynamics solvers.
#    The Cryosphere, 16, 689-709, doi:10.5194/tc-16-689-2022.

import os
import sys
import errno
import subprocess
import configparser
import numpy as np
import netCDF

from math import sqrt, sin, cos, tan, pi

# Parse the command line options
# ------------------------------
import argparse
parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# small helper function so argparse will understand unsigned integers
def unsigned_int(x):
    x = int(x)
    if x < 1:
        raise argparse.ArgumentTypeError("This argument is an unsigned int type! Should be an integer greater than zero.")
    return x

parser.add_argument('-c','--config', default='./slab.config', 
        help="The configure file used to setup the test case and run CISM")
parser.add_argument('-e','--executable', default='./cism_driver', 
        help="The CISM driver")
parser.add_argument('--hpc', nargs='?', const='aprun',
        help=" ".join(["Flag to Shortcut parallel run command lookup for High Performance Computing Systems.", 
                       "If flag apears without an argument, it will set run command to `aprun`,", 
                       "otherwise it will use the argument given."]))
parser.add_argument('-m', '--modifier', metavar='MOD', default='',
        help="Add a modifier to file names. FILE.EX will become FILE.MOD.EX")
parser.add_argument('-n','--parallel', metavar='N', type=unsigned_int, default=0, 
        help="Run in parallel using N processors.")
parser.add_argument('-o', '--output-dir',default='./output',
        help="Write all created files here.")
parser.add_argument('-q', '--quiet', action='store_true',
        help="Run the CISM process quietly.")
parser.add_argument('-s','--setup-only', action='store_true',
        help="Set up the test, but don't actually run it.")

# Additional options to set grid, solver, physics parameters, etc.:
#Note: The default mu_n = 0.0 is not actually used.
#      Rather, mu_n is computed below, unless mu_n > 0 is specified in the command line.
#      For n = 1, the default is mu_1 = 1.0e6 Pa yr.
parser.add_argument('-a','--approx', default='DIVA',
        help="Stokes approximation (SIALOC, SIA, SSA, BP, L1L2, DIVA, HYBRID)")
parser.add_argument('-beta','--beta', default=2000.0,
        help="Friction parameter beta (Pa (m/yr)^{-1})")
parser.add_argument('-dh','--delta_thck', default=0.0,
        help="Thickness perturbation (m)")
parser.add_argument('-dt','--tstep', default=0.01,
        help="Time step (yr)")
parser.add_argument('-gn','--glen_exponent', default=1,
        help="Exponent in Glen flow law")
parser.add_argument('-l','--levels', default=10,
        help="Number of vertical levels")
parser.add_argument('-mu','--mu_n', default=0.0,
        help="Viscosity parameter mu_n (Pa yr^{1/n})")
parser.add_argument('-nt','--n_tsteps', default=0,
        help="Number of timesteps")
parser.add_argument('-nx','--nx_grid', default=50,
        help="Number of grid cells in x direction")
parser.add_argument('-ny','--ny_grid', default=5,
        help="Number of grid cells in y direction")
parser.add_argument('-r','--resolution', default=100.0,
        help="Grid resolution (m)")
parser.add_argument('-theta','--theta', default=5.0,
        help="Slope angle (deg)")
parser.add_argument('-thk','--thickness', default=1000.0,
        help="Ice thickness (m)")


# Some useful functions
# ---------------------

# function to make a directory, and not worry if it exists.
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

# prep the command line functions
def prep_commands(args, config_name):
    driver = os.path.abspath(args.executable)
   
    quiet_mod = ''
    if args.quiet:
        quiet_mod = ' > '+config_name+'.oe'

    commands = []
    mkdir_p(args.output_dir)
    commands.append("cd "+os.path.abspath(args.output_dir))
    
    if args.hpc and (args.parallel > 0):
        mpiexec = args.hpc+' -n ' + str(args.parallel)+" "
    elif (args.parallel > 0):
        # These calls to os.system will return the exit status: 0 for success (the command exists), some other integer for failure
        if os.system('which openmpirun > /dev/null') == 0:
            mpiexec = 'openmpirun -np ' + str(args.parallel)+" "
        elif os.system('which mpirun > /dev/null') == 0:
            mpiexec = 'mpirun -np ' + str(args.parallel)+" "
        elif os.system('which aprun > /dev/null') == 0:
            mpiexec = 'aprun -n ' + str(args.parallel)+" "
        elif os.system('which mpirun.lsf > /dev/null') == 0:
            # mpirun.lsf does NOT need the number of processors
            mpiexec = 'mpirun.lsf '
        else:
            print("Unable to execute parallel run!")
            print("   Please edit the script to use your MPI run command, or run the model manually with")
            print("   something like: mpirun -np 4 ./cism_driver slab.config")
            sys.exit(1)
    else:
        mpiexec = ''

    commands.append(mpiexec+driver+" "+config_name+quiet_mod)

    return commands

# -----------------------------------
# Hard-cosed test case parameters
rhoi = 917.0     # kg/m^3
grav = 9.81      # m^2/s
baseElevation = 1000.0 # arbitrary height to keep us well away from sea level
  
# the main script function
# ------------------------
def main():
    """
    Run the slab test.
    """

    # check that file name modifier, if it exists, starts with a '-'
    if not (args.modifier == '') and not args.modifier.startswith('-') :
        args.modifier = '-'+args.modifier
         
    # get the configuration
    # ---------------------

    dx = float(args.resolution)
    dy = dx
    nx = int(args.nx_grid)
    ny = int(args.ny_grid)
    nz = int(args.levels)
    dt = float(args.tstep)
    tend = float(args.n_tsteps) * dt

    try:
        config_parser = configparser.ConfigParser()
        config_parser.read( args.config )
        file_name = config_parser.get('CF input', 'name')
        root, ext = os.path.splitext(file_name)

    except configparser.Error as error:
        print("Error parsing " + args.config )
        print("   "), 
        print(error)
        sys.exit(1)
    
    res=str(int(dx)).zfill(5)  # 00100 for 100m, 01000 for 1000m, etc.

    if args.parallel > 0:
        mod = args.modifier+'.'+res+'.p'+str(args.parallel).zfill(3)
    else:
        mod = args.modifier+'.'+res
    
    file_name = root+mod+ext
    config_name = root+mod+'.config'
    out_name = root+mod+'.out'+ext

    
    # Set up the domain
    # ----------------

    # Create the new config file
    # --------------------------
    if not args.quiet: 
        print("\nCreating config file: "+config_name)
    
    # Set grid and time parameters
    config_parser.set('grid', 'upn', str(nz))
    config_parser.set('grid', 'ewn', str(nx))
    config_parser.set('grid', 'nsn', str(ny))
    config_parser.set('grid', 'dew', str(dx))
    config_parser.set('grid', 'dns', str(dy))
    config_parser.set('time', 'dt',  str(dt))
    config_parser.set('time', 'tend',str(tend))

    # Set physics parameters that are needed to create the config file and the input netCDF file.
    # Note: rhoi and grav are hardwired above.

    # ice thickness
    thickness = float(args.thickness)

    # friction parameter beta (Pa (m/yr)^{-1})
    beta = float(args.beta)

    # basal inclination angle (degrees)
    theta = float(args.theta)
    theta_rad = theta * pi/180.0   # convert to radians

    # flow law exponent
    # Any value n >= 1 is supported.
    gn = float(args.glen_exponent)

    # Note: Fig. 3 in Dukowicz (2012) uses eta = 18 and theta = 18 deg.
    #       This gives u(1) = 10.0 * u(0), where u(1) = usfc and u(0) = ubas.

    # viscosity coefficient mu_n, dependent on n (Pa yr^{1/n})
    #      The nominal default is mu_n = 0.0, but this value is never used.
    #      If a nonzero value is specified on the command line, it is used;
    #        else, mu_n is computed here.  The goal is to choose a value mu_n(n) that
    #        will result in vertical shear similar to a default case with n = 1 and mu_1,
    #        provided we have similar values of H and theta.
    #      In the Dukowicz unpublished manuscript, the viscosity mu is given by
    #             mu = mu_n * eps_e^[(1-n)/n], where eps_e is the effective strain rate.
    #      For n = 1, we choose a default value of 1.0e6 Pa yr.
    #      For n > 1, we choose mu_n (units of Pa yr^{1/n}) to match the surface velocity
    #       we would have with n = 1 and the same values of H and theta.
    #      The general velocity solution is
    #      u(z') = u_b + du(z')
    #        where u_b = rhoi * grav * sin(theta) * cos(theta) / beta
    #          and du(z') = 2^{(1-n)/2}/(n+1) * sin^n(theta) * cos(theta)
    #                     * (rhoi*grav*H/mu_n)^n * H * [1 - (1 - z'/H)^{n+1}]
    #      For z' = H and general n, we have
    #             du_n(H) = 2^{(1-n)/2}/(n+1) * sin^n(theta) * cos(theta)
    #                     * (rhoi*grav*H/mu_n)^n * H
    #      For z' = H and n = 1, we have
    #             du_1(H) = (1/2) * sin(theta) * cos(theta) * (rhoi*grav*H/mu_1) * H
    #      If we equate du_n(H) with du_1(H), we can solve for mu_n:
    #              mu_n = [ 2^{(3-n)/(2n)}/(n+1) * sin^{n-1}(theta) * (rhoi*grav*H)^{n-1} * mu_1 ]^{1/n}
    #              with units Pa yr^{1/n}
    #      This value should give nearly the same shearing velocity du(H) for exponent n > 1
    #        as we would get for n = 1, given mu_1 and the same values of H and theta.

    if float(args.mu_n) > 0.0:
        mu_n = float(args.mu_n)
        mu_n_pwrn = mu_n**gn
    else:
        mu_1 = 1.0e6   # default value for mu_1 (Pa yr)
        mu_n_pwrn = 2.0**((3.0-gn)/2.0)/(gn+1.0) * sin(theta_rad)**(gn-1.0) \
                  * (rhoi*grav*thickness)**(gn-1.0) * mu_1   # (mu_n)^n
        mu_n = mu_n_pwrn**(1.0/gn)

    # Given mu_n, compute the temperature-independent flow factor A = 1 / [2^((1+n)/2) * mu_n^n].
    # This is how CISM incorporates a prescribed mu_n (with flow_law = 0, i.e. constant flwa).
    # Note: The complicated exponent of 2 in the denominator results from CISM and the Dukowicz papers
    #       having different conventions for the squared effective strain rate, eps_sq.
    #       In CISM:        mu = 1/2 * A^(-1/n) * eps_sq_cism^((1-n)/(2n))
    #        where eps_sq_cism = 1/2 * eps_ij * eps_ij
    #                   eps_ij = strain rate tensor
    #       In Dukowicz:    mu = mu_n * eps_sq_duk^((1-n)/(2n))
    #        where eps_sq_duk = eps_ij * eps_ij = 2 * eps_sq_cism
    #       Equating the two values of mu, we get mu_n * 2^((1-n)/(2n)) = 1/2 * A^(-1/n),
    #        which we solve to find A = 1 / [2^((1+n)/2) * mu_n^n]
    #       Conversely, we have  mu_n = 1 / [2^((1+n)/(2n)) * A^(1/n)]
    #TODO: Modify the Dukowicz derivations to use the same convention as CISM.
    flwa = 1.0 / (2.0**((1.0+gn)/2.0) * mu_n_pwrn)

    # Find the dimensionless parameter eta
    # This is diagnostic only; not used directly by CISM
    eta = (beta * thickness / mu_n**gn) * (rhoi * grav * thickness)**(gn-1)

    # periodic offset; depends on theta and dx
    offset = 1.0 * float(nx)*dx * tan(theta_rad)

    # Print some values
    print('nx   = ' + str(nx))
    print('ny   = ' + str(ny))
    print('nz   = ' + str(nz))
    print('dt   = ' + str(dt))
    print('tend = ' + str(tend))
    print('rhoi = ' + str(rhoi))
    print('grav = ' + str(grav))
    print('thck = ' + str(thickness))
    print('beta = ' + str(beta))
    print('gn   = ' + str(gn))
    print('mu_n = ' + str(mu_n))
    print('flwa = ' + str(flwa))
    print('eta  = ' + str(eta))
    print('theta   = ' + str(theta))
    print('offset  = ' + str(offset))

    # Write some options and parameters to the config file

    config_parser.set('parameters', 'periodic_offset_ew', str(offset))
    config_parser.set('parameters', 'rhoi', str(rhoi))
    config_parser.set('parameters', 'grav', str(grav))
    config_parser.set('parameters', 'n_glen', str(gn))
    config_parser.set('parameters', 'default_flwa', str(flwa))

    if (args.approx == 'SIALOC'):
        approx = -1             # Glissade local SIA; uses basal_tract_const = 1/beta
        config_parser.set('options', 'slip_coeff', str(1))
        config_parser.set('parameters', 'basal_tract_const', str(1./beta))
    elif (args.approx == 'SIA'):
        approx = 0              # Glissade matrix-based SIA
    elif (args.approx == 'SSA'):
        approx = 1
    elif (args.approx == 'BP'):
        approx = 2
    elif (args.approx == 'L1L2'):
        approx = 3
    elif (args.approx == 'DIVA'):
        approx = 4
    elif (args.approx == 'HYBRID'):
        approx = 5

    config_parser.set('ho_options', 'which_ho_approx', str(approx))

    config_parser.set('CF input', 'name', file_name)
    config_parser.set('CF output', 'name', out_name)
    config_parser.set('CF output', 'xtype', 'double')
    config_parser.set('CF output', 'frequency', str(tend))    # write output at start and end of run

    with open(config_name, 'w') as config_file:
        config_parser.write(config_file)

    # create the input netCDF file
    # ----------------------------
    if not args.quiet: 
        print("\nCreating slab netCDF file: "+file_name)
    try:
        nc_file = netCDF.NetCDFFile(file_name,'w',format='NETCDF3_CLASSIC')
    except TypeError:
        nc_file = netCDF.NetCDFFile(file_name,'w')

    nc_file.createDimension('time',1)
    nc_file.createDimension('x1',nx)
    nc_file.createDimension('y1',ny)
    nc_file.createDimension('level',nz)
    nc_file.createDimension('x0',nx-1) # staggered grid 
    nc_file.createDimension('y0',ny-1)

    x = dx*np.arange(nx,dtype='float32')
    y = dx*np.arange(ny,dtype='float32')

    nc_file.createVariable('time','f',('time',))[:] = [0]
    nc_file.createVariable('x1','f',('x1',))[:] = x
    nc_file.createVariable('y1','f',('y1',))[:] = y
    nc_file.createVariable('x0','f',('x0',))[:] = dx//2 + x[:-1] # staggered grid
    nc_file.createVariable('y0','f',('y0',))[:] = dy//2 + y[:-1]

    # Calculate values for the required variables.
    thk  = np.zeros([1,ny,nx],dtype='float32')
    topg = np.zeros([1,ny,nx],dtype='float32')
    artm = np.zeros([1,ny,nx],dtype='float32')
    unstagbeta = np.zeros([1,ny,nx],dtype='float32')

    # Calculate the geometry of the slab of ice
    # Note: Thickness is divided by cos(theta), since thck in CISM is the vertical distance
    #       from bed to surface.  On a slanted bed, this is greater than the distance measured
    #       in the direction perpendicular to the bed.
    #      Why does topg use a tan function?  Is the bed slanted?
    #      Do we need unstagbeta instead of beta?  Compare to ISMIP-HOM tests.

    thk[:] = thickness / cos(theta_rad)
    xmax = x[:].max()
    for i in range(nx):
        topg[0,:,i] = (xmax - x[i]) * tan(theta_rad) + baseElevation
    unstagbeta[:] = beta

    # Optionally, add a small perturbation to the thickness field

    if args.delta_thck:
        dh = float(args.delta_thck)
        for i in range(nx):

            # Apply a Gaussian perturbation, using the Box-Mueller algorithm:
            # https://en.wikipedia.org/wiki/Normal_distribution#Generating_values_from_normal_distribution

            mu = 0.0      # mean of normal distribution
            sigma = 1.0   # stdev of normal distribution

            rnd1 = np.random.random()   # two random numbers between 0 and 1
            rnd2 = np.random.random()

            # Either of the next two lines will sample a number at random from a normal distribution
            rnd_normal = mu + sigma * sqrt(-2.0 * np.log(rnd1)) * cos(2.0*pi*rnd2)
#            rnd_normal = mu + sigma * sqrt(-2.0 * np.log(rnd2)) * sin(2.0*pi*rnd1)

            dthk = dh * rnd_normal

            # Uncomment to make the perturbation a sine wave (alternating 1, -1)
            # This can be useful if we want reproducible results (no random numbers)
##            dthk = dh * sin((float(i) - 0.5)*pi)

            thk[0,:,i] = thk[0,:,i] + dthk
            print(i, dthk, thk[0,ny//2,i])
            thk_in = thk   # for comparing later to final thk

    # Create the required variables in the netCDF file.
    nc_file.createVariable('thk', 'f',('time','y1','x1'))[:] = thk
    nc_file.createVariable('topg','f',('time','y1','x1'))[:] = topg
    nc_file.createVariable('unstagbeta','f',('time','y1','x1'))[:] = unstagbeta

    nc_file.close()
    mkdir_p(args.output_dir)
    subprocess.check_call("cp *rilinosOptions.xml "+args.output_dir, shell=True)
    subprocess.check_call("mv "+file_name+" "+args.output_dir, shell=True)
    subprocess.check_call("mv "+config_name+" "+args.output_dir, shell=True)


    # Run CISM
    # --------
    command_list = prep_commands(args, config_name)
    commands_all = ["# SLAB"+mod+" test"]
    commands_all.extend( command_list )
    
    result_mv = "mv results "+root+mod+".results 2>/dev/null"
    timing_mv = "for file in cism_timing*; do mv $file "+root+mod+".$file 2>/dev/null; done"
    commands_all.append(result_mv)
    commands_all.append(timing_mv)
    commands_all.append(" ")
    
    if not args.setup_only:
        if not args.quiet: 
            print("\nRunning CISM slab test")
            print(  "======================\n")

            print('command_list =' + str(command_list))

        process = subprocess.check_call(str.join("; ",command_list), shell=True)
   
        try:
            subprocess.check_call("cd "+args.output_dir+"; "+result_mv, shell=True)
        except subprocess.CalledProcessError:
            pass 

        try:
            subprocess.check_call("cd "+args.output_dir+"; "+timing_mv, shell=True)
        except subprocess.CalledProcessError:
            pass

        if not args.quiet: 
            print("\nFinished running the CISM slab test")
            print(  "===================================\n")

    else:
        run_script = args.output_dir+os.sep+root+mod+".run" 
        
        with open(run_script,'w') as run_file:
            run_file.write('#!/usr/bin/env bash \n \n')
            for command in commands_all:
                run_file.write(command+" \n")

        os.chmod(run_script, 0o755)   # uses an octal number!

        if not args.quiet:
            print("\nFinished setting up the CISM slab test")
            print(  "======================================")
            print(  "   To run the test, use: "+run_script)


# Run only if this is being run as a script.
if __name__=='__main__':
    
    # get the command line arguments
    args = parser.parse_args()
    
    # run the script
    sys.exit(main())
