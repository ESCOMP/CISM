#!/usr/bin/env python

"""
This script plots the results of an experiment with an ice "slab" on an
inclined plane. Test case is described in sections 5.1-2 of:
    J.K. Dukowicz, 2012. Reformulating the full-Stokes ice sheet model for a
    more efficient computational solution. The Cryosphere, 6, 21-34.
    www.the-cryosphere.net/6/21/2012/

Blatter-Pattyn First-order solution is described in J.K. Dukowicz, manuscript
in preparation.
"""
#FIXME: Manuscript likely not in prep anymore -- JHK, 08/07/2015
#       Not published as of July 2021 -- WHL

# Written by Matt Hoffman, Dec. 16, 2013
# Reconfigured by Joseph H Kennedy at ORNL on August 7, 2015 to work with the regression testing
#     NOTE: Did not adjust inner workings except where needed.
# Revised by William Lipscomb in 2021 to support more options, including general values of Glen's n.

import os
import sys
import glob

import numpy as np
import matplotlib.pyplot as plt

from netCDF import *
from math import tan, pi, sin, cos, atan

# Get hard-coded parameters from the run script.
from runSlab import rhoi, grav

import ConfigParser

import argparse
parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# The command line options
# ------------------------
parser.add_argument('-o', '--output-dir', default='./output',
        help="Directory containing the tests output files. Warning: if there is a" \
        +"path passed via the -f/--out-file option, this argument will be" \
        +"ignored.")

parser.add_argument('-f', '--output-file',
        help="The tests output file you would like to plot. If a path is" \
        +"passed via this option, the -o/--output-dir option will be ignored.")

parser.add_argument('-c','--config-file',
        help="The configure file used to set up the test case and run CISM")


# ===========================================================
# Define some variables and functions used in the main script
# ===========================================================

#WHL args.output-file with a hyphen?
def get_in_file():
    if args.output_file:
        out_d, out_f = os.path.split(args.output_file)
        if out_d:
            args.output_dir = out_d
            args.output_file = out_f
    
        print("\nUsing "+os.path.join(args.output_dir, args.output_file)+"\n")
        
    else:
        outpath = os.path.join(args.output_dir, '*.out.nc')
        matching = glob.glob(outpath)
        if len(matching) == 1:
            newest = matching[0]
            print("\nUsing "+newest+"\n")
       
        elif len(matching) > 1:
            newest = max(matching, key=os.path.getmtime)
            print("\nWARNING: MULTIPLE *.out.nc FILES DETECTED!")
            print(  "==========================================")
            print(  "Plotting the most recently modified file in the output directory:")
            print(  "    "+newest)
            print(  "To plot another file, specify it with the -f/--outfile option.\n")
            
        else:
            print("\nERROR: NO *.out.nc FILES DETECTED!")
            print(  "==================================")
            print(  "Either specify a location to look for the test output")
            print(  "files with the -o/--output-dir option, or the test output")
            print(  "file with the -f/--output-file option.\n")
            sys.exit(1)

        args.output_file = os.path.basename(newest)

    filein = NetCDFFile(os.path.join(args.output_dir, args.output_file),'r')
     
    return filein

def split_file_name(file_name):
    """
    Get the root name, size, and number of processors from an out.nc filename.
    #WHL - Adapted from plotISMIP_HOM.py
    """
    root = ''
    size = ''
    proc = ''

    file_details = file_name.replace('.out.nc','') .split('.')
#    print(file_details)
#    print('len = ' + str(len(file_details)))

    if len(file_details) > 2:
        proc = '.'+file_details[2]
    size = '.'+file_details[1]
    root = file_details[0]

    return (root, size, proc)

# =========================
# Actual script starts here
# =========================
def main():
    """
    Plot the slab test results.
    """

    filein = get_in_file()

    # Get needed variables from the output file
    x1 = filein.variables['x1'][:]
    y1 = filein.variables['y1'][:]
    x0 = filein.variables['x0'][:]
    y0 = filein.variables['y0'][:]
    level = filein.variables['level'][:]

    # =====================================================================
    # Choose where you want the velocity profile to be taken
    # use integer floor division operator to get an index close to the center 
    xp = len(x0)//2
    yp = len(y0)//2
    # =====================================================================

    print('Using x index of '+str(xp)+'='+str(x0[xp]))
    print('Using y index of '+str(yp)+'='+str(y0[yp]))

    thk = filein.variables['thk'][:]
    if netCDF_module == 'Scientific.IO.NetCDF':
        thk = thk * filein.variables['thk'].scale_factor
    topg = filein.variables['topg'][:]
    if netCDF_module == 'Scientific.IO.NetCDF':
        topg = topg * filein.variables['topg'].scale_factor
    uvel = filein.variables['uvel'][:]
    if netCDF_module == 'Scientific.IO.NetCDF':
        uvel = uvel * filein.variables['uvel'].scale_factor
    beta_2d = filein.variables['beta'][:]
    if netCDF_module == 'Scientific.IO.NetCDF':
        beta_2d = beta_2d * filein.variables['beta'].scale_factor

    # Get the name of the config file
    # If not entered on the command line, then construct from the output filename

    if not args.config_file:
        root, size, proc = split_file_name(args.output_file)
        args.config_file = root + size + proc + '.config'

    configpath = os.path.join(args.output_dir, args.config_file)
    print('configpath = ' + configpath)

    # Get gn and default_flwa from the config file

    try:
        config_parser = ConfigParser.SafeConfigParser()
        config_parser.read( configpath )

        gn = float(config_parser.get('parameters','n_glen'))
        flwa = float(config_parser.get('parameters', 'default_flwa'))

    except ConfigParser.Error as error:
        print("Error parsing " + args.config )
        print("   "),
        print(error)
        sys.exit(1)

    # Derive the viscosity constant mu_n from flwa
    # This expression is derived in the comments on flwa in runSlab.py.
    mu_n = 1.0 / (2.0**((1.0+gn)/(2.0*gn)) * flwa**(1.0/gn))

    # Get the ice thickness from the output file.
    # If thickness = constant (i.e., the optional perturbation dh = 0), it does not matter where we sample.
    # Note: In general, this thickness will differ from the baseline 'thk' that is used in runSlab.py
    #        to create the input file.
    #       This is because the baseline value is measured perpendicular to the sloped bed,
    #        whereas the CISM value is in the vertical direction, which is not perpendicular to the bed.
    thickness = thk[0,yp,xp]

    # Get beta from the output file.
    # Since beta = constant, it does not matter where we sample.
    beta = beta_2d[0,yp,xp]

    # Derive theta from the output file as atan(slope(topg))
    # Since the slope is constant, it does not matter where we sample.
    slope = (topg[0,yp,xp] - topg[0,yp,xp+1]) / (x0[xp+1] - x0[xp])
    thetar = atan(slope)
    theta = thetar * 180.0/pi

    # Compute the dimensionless parameter eta and the velocity scale,
    # which appear in the scaled velocity solution.
    eta = (beta * thickness / mu_n**gn) * (rhoi * grav * thickness)**(gn-1)
    velscale = (rhoi * grav * thickness / mu_n)**gn * thickness

    print('gn   = ' + str(gn))
    print('rhoi = ' + str(rhoi))
    print('grav = ' + str(grav))
    print('thck = ' + str(thickness))
    print('mu_n = ' + str(mu_n))
    print('flwa = ' + str(flwa))
    print('beta = ' + str(beta))
    print('eta  = '  + str(eta))
    print('theta= ' + str(theta))
    print('velscale = ' + str(velscale))

    # === Plot the results at the given location ===
    # Note we are not plotting like in Fig 3 of paper.
    # That figure plotted a profile against zprime.
    # It seemed more accurate to plot a profile against z to avoid interpolating model results (analytic solution can be calculated anywhere).
    # Also, the analytic solution calculates the bed-parallel u velocity, but CISM calculates u as parallel to the geoid,
    #  so we need to transform the analytic solution to the CISM coordinate system.

    #WHL - I think the analytic solution is actually for u(z'), which is not bed-parallel.
    #      The bed-parallel solution would be u'(z'), with w'(z') = 0.

    fig = plt.figure(1, facecolor='w', figsize=(12, 6))

    # define the scaled x & z variables, with an origin at the bed at this cell
    z = (1.0-level)*thk[0,yp,xp]/thickness  # elevation / thickness  (with the Dukowicz coord. sys.)
    #print 'z', z
    x = (x0-x0[xp]) / thickness
    # calculate rotated zprime coordinates for this column (we assume the solution truly is spatially uniform)
    zprime = x[xp] * sin(thetar) + z * cos(thetar)

    # Calculate analytic solution for x-component of velocity (eq. 39 in paper) for the CISM-column
    uvelStokesAnalyticScaled = sin(thetar) * cos(thetar) / eta   \
        - 2**((1.0-gn)/2.0) * sin(thetar)**gn * cos(thetar) / (gn+1) * ( (1.0 - zprime)**(gn+1) - 1.0 )

    # Calculate the BP FO solution for x-component of velocity (Dukowicz, in prep. paper, Eq.30, n=3)
    uvelFOAnalyticScaled = + tan(thetar) / eta  \
                           -  2**((1.0-gn)/2.0) * tan(thetar)**gn /  \
                           ( (gn + 1) * (1.0 + 3.0 * sin(thetar)**2)**((gn+1.0)/2.0) )  \
                           * ( (1.0 - zprime)**(gn+1) - 1.0 )

    ### 1. Plot as nondimensional variables
    # Plot analytic solution 
    fig.add_subplot(1,2,1)
    plt.plot(uvelStokesAnalyticScaled, z, '-kx', label='Analytic Stokes')
    plt.plot(uvelFOAnalyticScaled, z, '-ko', label='Analytic FO')

    # Plot model results
    plt.plot(uvel[0,:,yp,xp] / velscale, z, '--ro', label='CISM') 
    plt.ylim((-0.05, 1.05))
    #plt.gca().invert_yaxis()  # put 0.0 (surface) on top & 1.0 (bed) on bottom
    plt.legend(loc='best')
    plt.xlabel('Nondimensional velocity')
    plt.ylabel('Nondimenstional, unrotated vertical coordinate')
    plt.title('Velocity profile at x=' + str(x0[xp]) + ' m, y=' + str(y0[yp]) + ' m\n(Scaled coordinates)')

    ### 2. Plot as dimensional variables
    # Plot analytic solution for x-component of velocity (eq. 39 in paper)  
    fig.add_subplot(1,2,2)
    plt.plot(uvelStokesAnalyticScaled * velscale, z * thk[0,yp,xp] + topg[0,yp,xp], '-kx', label='Analytic Stokes')
    plt.plot(uvelFOAnalyticScaled * velscale, z * thk[0,yp,xp] + topg[0,yp,xp], '-ko', label='Analytic FO')
    # Plot model results
    plt.plot(uvel[0,:,yp,xp], z * thk[0,yp,xp] + topg[0,yp,xp], '--ro', label='CISM') 
    plt.legend(loc='best')
    plt.xlabel('velocity (m/yr)')
    plt.ylabel('elevation (m)')
    plt.title('Velocity profile at x=' + str(x0[xp]) + ' m, y=' + str(y0[yp]) + ' m\n(Unscaled coordinates)')

    #################
#    print('y0_min:')
#    print(y0.min())
#    print('y0_max:')
#    print(y0.max())

    # Now plot maps to show if the velocities vary over the domain (they should not)
    # For some reason, the y-axis has a greater extent than the range (y0.min, y0.max).
    #TODO - Fix the y-axis extent.  Currently, the extent is too large for small values of ny.
    #TODO - Plot the thickness relative to the initial thickness.

    fig = plt.figure(2, facecolor='w', figsize=(12, 6))
    fig.add_subplot(1,2,1)
    uvelDiff = uvel[0,0,:,:] - uvel[0,0,yp,xp]
    tol = 1.0e-12
    #plt.pcolor(x0,y0,uvelDiff,vmin=uvelDiff.min()-tol, vmax=uvelDiff.max()+tol)
    plt.imshow(uvelDiff,extent=(x0.min(), x0.max(), y0.min(), y0.max()),vmin=uvelDiff.min()-tol, vmax=uvelDiff.max()+tol, interpolation="none")
    plt.xlabel('x position (m)')
    plt.ylabel('y position (m)')
    plt.colorbar()
    plt.plot(x0[xp],y0[yp],'k*',markersize=9)
    plt.title('Map of difference of surface x-velocity from\nvalue at reference location (m/yr)')
    plt.axis('equal')

    fig.add_subplot(1,2,2)
    uvelDiff = uvel[0,-1,:,:] - uvel[0,-1,yp,xp]
    tol = 1.0e-12
    #plt.pcolor(x0,y0,uvelDiff,vmin=uvelDiff.min()-tol, vmax=uvelDiff.max()+tol)
    plt.imshow(uvelDiff,extent=(x0.min(), x0.max(), y0.min(), y0.max()),interpolation="none",vmin=uvelDiff.min()-tol, vmax=uvelDiff.max()+tol)
    plt.xlabel('x position (m)')
    plt.ylabel('y position (m)')
    plt.colorbar()
    plt.plot(x0[xp],y0[yp],'k*',markersize=9)
    plt.title('Map of difference of basal x-velocity from\nvalue at reference location (m/yr)')
    plt.axis('equal')


    # Optional plot for comparing analytic solution to Fig. 3 in the paper (model output is not plotted)
    #fig = plt.figure(9, facecolor='w')
    #plt.plot(level, sin(thetar)**3 * cos(thetar) / 8.0 * (1.0 - (level-1.0)**4 ) + sin(thetar)*cos(thetar)/eta , 'b-', label='nonlinear stokes')
    #plt.plot(level, tan(thetar)**3 / (8.0 * (1.0 + 3.0 * sin(thetar)**2)**2) * (1.0 - (level-1.0)**4 ) + tan(thetar)/eta, 'b--' , label='nonlinear fo')
    #plt.ylim((0.0, 0.04)); plt.xlabel("z'"); plt.ylabel('u'); plt.legend()

    plt.draw()
    plt.show()

    filein.close()

# Run only if this is being run as a script.
if __name__=='__main__':
    
    # get the command line arguments
    args = parser.parse_args()
    
    # run the script
    sys.exit(main())
