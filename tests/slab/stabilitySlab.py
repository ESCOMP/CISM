#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script runs a series of CISM experiments at different resolutions.
At each resolution, it determines the maximum stable time step.
A run is deemed to be stable if the standard deviation of a small thickness perturbation
decreases during a transient run (100 timesteps by default).
 
Used to obtain the CISM stability results described in:
    Robinson, A., D. Goldberg, and W. H. Lipscomb, 2022, A comparison of the
    stability and performance of depth-integrated ice-dynamics solvers.
    The Cryosphere, 16, 689-709, doi:10.5194/tc-16-689-2022.
"""

# Authors
# -------
# Created by William Lipscomb, July 2021

import os
import sys
import errno
import subprocess

import numpy as np
import netCDF
from math import sqrt, log10

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

# The following command line arguments determine the set of resolutions to run the slab test.
# At each resolution, we aim to find the maximum stable time step.
# Note: If args.n_resolution > 1, then args.resolution (see below) is ignored.

parser.add_argument('-nr','--n_resolution', default=1,
        help="number of resolutions")
parser.add_argument('-rmin','--min_resolution', default=10.0,
        help="minimum resolution (m)")
parser.add_argument('-rmax','--max_resolution', default=40000.0,
        help="minimum resolution (m)")

# The following command line arguments are the same as in runSlab.py.
# Not sure how to avoid code repetition.

parser.add_argument('-c','--config', default='./slab.config', 
        help="The configure file used to setup the test case and run CISM")
parser.add_argument('-e','--executable', default='./cism_driver', 
        help="The CISM driver")
parser.add_argument('-m', '--modifier', metavar='MOD', default='',
        help="Add a modifier to file names. FILE.EX will become FILE.MOD.EX")
parser.add_argument('-n','--parallel', metavar='N', type=unsigned_int, default=0, 
        help="Run in parallel using N processors.")
parser.add_argument('-o', '--output_dir',default='./output',
        help="Write all created files here.")
parser.add_argument('-a','--approx', default='BP',
        help="Stokes approximation (SIA, SSA, BP, L1L2, DIVA)")
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
        help="Ice thickness")

                        ############
                        # Functions #
                        ############

def reading_file(inputFile):

    #Check whether a netCDF file exists, and return a list of times in the file

    ReadVarFile = True 
    try:
        filein = netCDF.NetCDFFile(inputFile,'r')
        time = filein.variables['time'][:]
        filein.close()
        print('Was able to read file ' + inputFile)
        print(time)
    except:
        ReadVarFile = False
        time = [0.]
        print('Was not able to read file' + inputFile)

    return time, ReadVarFile


def check_output_file(outputFile, time_end):

    # Check that the output file exists with the expected final time slice

    # Path to experiment
    path_to_slab_output = './output/'

    # File to check
    filename = path_to_slab_output + outputFile

    # Read the output file
    print('Reading file ' + str(filename))
    time_var, VarRead = reading_file(filename)

#    print(time_var)

    # Checking that the last time entry is the same as we expect from time_end
    # Allow for a small roundoff difference.
    if abs(time_var[-1] - time_end) < 1.0e-7:
        check_time_var = True
    else:
        check_time_var = False

    print('time_end = ' + str(time_end))
    print('last time in file = ' + str(time_var[-1]))

    # Creating the status of both checks 
    check_passed = check_time_var and VarRead

    if check_passed:
        print('Found output file with expected file time slice')
    else:
        if (not VarRead):
            print('Output file cannot be read')
        else:
            if not check_time_var:
                print('Output file is missing time slices')
    
    return check_passed


def main():

    print('In main')

    """
    For each of several values of the horizontal grid resolution, determine the maximum
    stable time step for a given configuration of the slab test.
    """

    resolution = []

    # Based on the input arguments, make a list of resolutions at which to run the test.
    # The formula and the default values of rmin and rmax give resolutions agreeing with
    #  those used by Alex Robinson for Yelmo, for the case nres = 12:
    #  resolution = [10., 21., 45., 96., 204., 434., 922., 1960., 4170., 8850., 18800., 40000.]

    print('Computing resolutions')
    print(args.n_resolution)
    if int(args.n_resolution) > 1:
        nres = int(args.n_resolution)
        resolution = [0. for n in range(nres)]
        rmin = float(args.min_resolution)
        rmax = float(args.max_resolution)
        for n in range(nres):
            res = 10.0**(log10(rmin) + (log10(rmax) - log10(rmin))*float(n)/float(nres-1))
            # Round to 3 significant figures (works for log10(res) < 5)
            if log10(res) > 4.:
                resolution[n] = round(res, -2)
            elif log10(res) > 3.:
                resolution[n] = round(res, -1)
            else:
                resolution[n] = round(res)
    else:
        nres = 1
        resolution.append(float(args.resolution))

    print('nres = ' + str(nres))
    print(resolution)

    # Create an array to store max time step for each resolution
    rows, cols = (nres, 2)
    res_tstep = [[0. for i in range(cols)] for j in range(rows)]
    for n in range(nres):
        res_tstep[n][0] = resolution[n]

    for n in range(nres):

        print('output_dir: ' + args.output_dir)

        # Construct the command for calling the main runSlab script
        run_command = 'python runSlab.py'
        run_command = run_command + ' -c '  + args.config
        run_command = run_command + ' -e '  + args.executable
        if args.parallel > 0:
            run_command = run_command + ' -n '  + str(args.parallel)
        run_command = run_command + ' -o '  + args.output_dir
        run_command = run_command + ' -a '  + args.approx
        run_command = run_command + ' -beta ' + str(args.beta)
        run_command = run_command + ' -dh ' + str(args.delta_thck)
        run_command = run_command + ' -gn ' + str(args.glen_exponent)
        run_command = run_command + ' -l '  + str(args.levels)
        run_command = run_command + ' -mu ' + str(args.mu_n)
        run_command = run_command + ' -nt ' + str(args.n_tsteps)
        run_command = run_command + ' -nx ' + str(args.nx_grid)
        run_command = run_command + ' -ny ' + str(args.ny_grid)
        run_command = run_command + ' -theta '+ str(args.theta)
        run_command = run_command + ' -thk '+ str(args.thickness)

        tend = float(args.n_tsteps) * args.tstep

        res = resolution[n]
        run_command = run_command + ' -r '  + str(res)
         
        # Choose the time step.
        # Start by choosing a very small timestep that can be assumed stable
        #  and a large step that can be assumed unstable.
        # Note: SIA-type solvers at 10m resolution can require dt <~ 1.e-6 yr.

        tstep_lo = 1.0e-7
        tstep_hi = 1.0e+5
        tstep_log_precision = 1.0e-4
        print('Initial tstep_lo = ' + str(tstep_lo))
        print('Initial tstep_hi = ' + str(tstep_hi))
        print('Log precision = ' + str(tstep_log_precision))

        while (log10(tstep_hi) - log10(tstep_lo)) > tstep_log_precision:

            # Compute the time step as the geometric mean of the tstep_lo and tstep_hi.
            # tstep_lo is the largest time step known to be stable.
            # tstep_hi is the smallest time step known to be unstable.

            tstep = sqrt(tstep_lo*tstep_hi)
    
            run_command_full = run_command + ' -dt ' + str(tstep)

            print("\nRunning CISM slab test...")
            print('resolution = ' + str(res))
            print('tstep = ' + str(tstep))
            print('run_command = ' + run_command_full)

            process = subprocess.check_call(run_command_full, shell=True)

            print("\nFinished running the CISM slab test")

            # Determine the name of the output file.
            # Must agree with naming conventions in runSlab.py

            file_name = args.config
            root, ext = os.path.splitext(file_name)

            res=str(int(res)).zfill(5)  # 00100 for 100m, 01000 for 1000m, etc.

            if args.parallel > 0:
                mod = args.modifier + '.' + res + '.p' + str(args.parallel).zfill(3)
            else:
                mod = args.modifier + '.' + res

            outputFile = root + mod + '.out.nc'

            # Check whether the output file exists with the desired final time slice.

            time_end = float(args.n_tsteps) * tstep

            print('outputFile = ' + str(outputFile))
            print('n_tsteps = ' + str(float(args.n_tsteps)))
            print('tstep = ' + str(tstep))
            print('time_end = ' + str(time_end))

            check_passed = check_output_file(outputFile, time_end)

            if check_passed:

                print('Compute stdev of initial and final thickness for j = ny/2')
                nx = int(args.nx_grid)
                ny = int(args.ny_grid)

                # Read initial and final thickness from output file
                outpath = os.path.join(args.output_dir, outputFile)
                print('outpath = ' + outpath)
                filein = netCDF.NetCDFFile(outpath,'r')
                thk = filein.variables['thk'][:]

                j = ny//2
                thk_in = thk[0,j,:]
                thk_out = thk[1,j,:]

                # Compute <H>
                Hav_in = 0.0
                Hav_out = 0.0
                for i in range(nx):
                    Hav_in = Hav_in + thk_in[i]
                    Hav_out = Hav_out + thk_out[i]
                Hav_in = Hav_in / nx
                Hav_out = Hav_out / nx

                # Compute <H^2>
                H2av_in = 0.0
                H2av_out = 0.0
                for i in range(nx):
                    H2av_in = H2av_in + thk_in[i]**2
                    H2av_out = H2av_out + thk_out[i]**2
                H2av_in = H2av_in / nx
                H2av_out = H2av_out / nx

                print('H2av_out =' + str(H2av_out))
                print('Hav_out^2 =' + str(Hav_out**2))

                # Compute stdev = sqrt(<H^2> - <H>^2)
                var_in  = H2av_in  - Hav_in**2
                var_out = H2av_out - Hav_out**2

                if var_in > 0.:
                    stdev_in = sqrt(H2av_in - Hav_in**2)
                else:
                    stdev_in = 0.

                if var_out > 0.:
                    stdev_out = sqrt(H2av_out - Hav_out**2)
                else:
                    stdev_out = 0.

                if stdev_in > 0.:
                    ratio = stdev_out/stdev_in
                else:
                    ratio = 0.

                print('stdev_in  = ' + str(stdev_in))
                print('stdev_out = ' + str(stdev_out))
                print('ratio = ' + str(ratio))

                # Determine whether the run was stable.
                # A run is defined to be stable if the final standard deviation of thickness
                #  is less than the initial standard deviation

                if ratio < 1.:
                    tstep_lo = max(tstep_lo, tstep)
                    print('Stable, new tstep_lo =' + str(tstep_lo))
                else:
                    tstep_hi = min(tstep_hi, tstep)
                    print('Unstable, new tstep_hi =' + str(tstep_hi))

            else:   # check_passed = F; not stable
                tstep_hi = min(tstep_hi, tstep)
                print('Unstable, new tstep_hi =' + str(tstep_hi))

            print('Latest tstep_lo = ' + str(tstep_lo))
            print('Latest tstep_hi = ' + str(tstep_hi))

            # Add to the array containing the max stable timestep at each resolution.
            # Take the max stable timestep to be the average of tstep_lo and tstep_hi.
            res_tstep[n][1] = 0.5 * (tstep_lo + tstep_hi)

            print('New res_tstep, res #' + str(n))
            print(res_tstep)

    # Print a table containing the max timestep for each resolution
    for n in range(nres):
        float_res = res_tstep[n][0]
        float_dt  = res_tstep[n][1]
        formatted_float_res = "{:8.1f}".format(float_res)
        formatted_float_dt =  "{:.3e}".format(float_dt)  # exponential notation with 3 decimal places
        print(formatted_float_res + '    ' + formatted_float_dt)

# Run only if this is being run as a script.                                                                                                       
if __name__=='__main__':

    # get the command line arguments                                                                                                               
    args = parser.parse_args()

    # run the script
    sys.exit(main())
