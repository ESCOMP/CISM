#!/usr/bin/env python

# This script runs any or all of the various MISMIP3d experiments.
# See this paper for details about the MISMIP3d experiment:
# Pattyn et al., Grounding-line migration in plan-view marine ice-sheet models: results of the ice2sea MISMIP3d intercomparison, J. Glaciol., 59, doi:10.3189/2013JoG12J129, 2013.

import sys, os
import shutil
import fileinput
import numpy as np
from argparse import ArgumentParser
import configparser 
from netCDF4 import Dataset


###################################
# string constants and parameters #
###################################



####################################
# Function used later in the code #
####################################


# This function returns the basal friction perturbation for experiment P75S.
# Note: x, y and xGL (the grounding line position) need to be in km.
def computeBasalPerturbation(x,y,xGL):
    a   = 0.75  # perturbation amplitude
    xc  = 150   # in km
    yc  = 10    # in km
    yGL = 50    # in km (The function is center at the center line of the domain)

    # Calculating the spatial variation.
    xpart    = -(x-xGL)**2/(2*xc**2)
    ypart    = -(y-yGL)**2/(2*yc**2)
    Cspatial = 1-a*np.exp(xpart + ypart)
    return Cspatial


########
# Code #
########


# Parse options.
parser = ArgumentParser()

parser.add_argument('-e', '--exec',     dest='executable', type =str, default='./cism_driver', help='Path to the CISM executable')
parser.add_argument('-x', '--expt',     dest='experiment', type=str,  default='all', help='MISMIP3d experiment(s) to run', metavar='EXPT')
parser.add_argument('-n', '--parallel', dest='parallel',   type=int, help='Number of processors: if specified then run in parallel', metavar='NUMPROCS')

options = parser.parse_args()

if options.experiment == 'all':
    experiments = ['Stnd', 'P75S', 'P75R']
    print( 'Run all the MISMIP3d experiments')
elif options.experiment == 'allP75':
    experiments = ['P75S', 'P75R']
    print( 'Run all the MISMIP3d P75 (S) and (R) experiments')
elif options.experiment in ['Stnd', 'P75S', 'P75R']:
    experiments = [options.experiment]
    print( 'Run experiment', options.experiment)
else:
    sys.exit('Please specify experiment(s) from this list: Stnd, P75S, P75R')


restart = 0  # constant meant to keep track of restart status.

# Loop through experiments.
for expt in experiments:

    # Change to directory for this experiment.
    os.chdir(expt)
    print( 'Changed directory to ', expt)

    # Set and open config file for reading.
    configfile = 'mismip3d' + expt + '.config'
    config = configparser.ConfigParser(delimiters=('=', ':'),
                            comment_prefixes=('#', ';'),
                            inline_comment_prefixes=';',
                            interpolation=None)
    config.read(configfile)

    inputFile   = config.get('CF input',   'name')
    outputFile  = config.get('CF output',  'name')
    outputFreq  = config.get('CF output',  'frequency')
    endTime     = config.get('time',       'tend')
    endTime     = float(endTime)
    
    # Time buffer to remedy the restart pb when switching to different time step within a run.
    buffer = float(outputFreq) - 1.

    # Read output file content information.
    lastTimeEntry     = 0
    lastEntryInternal = 0
    sizeTimeOutput    = 0
    if os.path.exists(outputFile):
        outputData = Dataset(outputFile,'r')
        if outputData['time'].size != 0:
            sizeTimeOutput = outputData['time'].size
            lastTimeEntry  = outputData['time'][-1]
        
        outputData.close()

    # Take action based on the stage of the experimental run.
    if (lastTimeEntry >= (endTime - buffer)) and (sizeTimeOutput > 1):
        # The run for this A value is done, moving to the next one.
        pass
    elif (lastTimeEntry < endTime) and (sizeTimeOutput > 1):
        # The run for this A value is not done and needs to continue.
        
        # Make sure restart is set to 1 in config file.
        config.set('options', 'restart', 'standard')
        
        # Write to config file.
        with open(configfile, 'w') as newconfigfile:
            config.write(newconfigfile)
    
    else:
        print( 'there is nothing to restart from, executing from the beginning')


    # Make sure we are starting from the final time slice of the input file,
    # (Except for Stnd, the input file is a restart file from a previous run.)



    if expt != 'Stnd':

        # Edit the 'time' entry in the [CF input] section.
        inputfile = config.get('CF input', 'name')
        try:
           inputfileOut = config.get('CF input', 'nameout')
        except:
           pass  # no nameout name in config file
        
        infile = Dataset(inputfile,'r')
        ntime = len(infile.dimensions['time']) # reading the last time slice of the previous experiment
        config.set('CF input', 'time', str(ntime))
        print( 'input file =', inputfile)
        print( 'time slice =', str(ntime))
       
        infile.close()


        # For experiment P75S we need to induce the perturbation and modify C_space_factor.
        if (expt == 'P75S') and (restart==0):
            
            # We need to calculate the grounding line location from the input file of the experiment at the center line location.
            infileOut = Dataset(inputfileOut,'r')
            ntimeout  = len(infileOut.dimensions['time']) # reading last time slice of previous ext output file
            f_ground  = infileOut.variables['f_ground'][:,:,:]
            infileOut.close()

            print( 'opening inputfile ',inputfile)
            infile = Dataset(inputfile,'r+')
            nx = len(infile.dimensions['x0'])
            ny = len(infile.dimensions['y0'])
            x0 = infile.variables['x0'][:]
            y0 = infile.variables['y0'][:]

            ycenter = int(ny/2)
            print( 'ycenter=',ycenter)
                
            for i in range(1,nx-1):
                if (f_ground[ntimeout-1,ycenter,i] > 0) and (f_ground[ntimeout-1,ycenter,i+1] == 0):
                     xGL = x0[i]
                     print( 'grounding line at ycenter =', xGL)
        
            # We now write C_space factor to the file.
            print( 'writing peturbation to C_space_factor')
            for i in range(nx):
                for j in range(ny):
                     infile.variables['c_space_factor'][-1,j,i] = computeBasalPerturbation(x0[i]/1000,y0[j]/1000,xGL/1000)

            # We first check that the variable "internal_time" from the Stnd restart file has 0 for its last entry.
            lastentry = infile['internal_time'][-1]
            if lastentry != 0:
                infile['internal_time'][:] = infile['internal_time'][:] - lastentry
                print( 'the new internal_time array is ', infile['internal_time'][:])
            
            infile.close()

    
        # We need to reset C_space_factor to 1 and internal_time to 0.
        if (expt == 'P75R') and (restart==0):
             infile = Dataset(inputfile,'r+')
             print( 'opening inputfile ',inputfile)
             infile.variables['c_space_factor'][:,:,:] = 1

             lastentry = infile['internal_time'][-1]
             if lastentry != 0:
                 infile['internal_time'][:] = infile['internal_time'][:] - lastentry
                 print( 'the new internal_time array is ', infile['internal_time'][:])

             infile.close()
        

        # Write the modified config file.
        with open(configfile, 'w') as newconfigfile:
            config.write(newconfigfile)



    # Run CISM.

    print( 'parallel =', options.parallel)

    if options.parallel == None:
        # Perform a serial run.
        os.system(options.executable + ' ' + configfile)
    else:
        # Perform a parallel run.
        if options.parallel <= 0:
            sys.exit( 'Error: Number of processors specified for parallel run is <=0.' )
        else:
            # These calls to os.system will return the exit status: 0 for success (the command exists), some other integer for failure.
            if os.system('which openmpirun > /dev/null') == 0:
                mpiexec = 'openmpirun -np ' + str(options.parallel)
            elif os.system('which mpirun > /dev/null') == 0:
                mpiexec = 'mpirun -np ' + str(options.parallel)
            elif os.system('which aprun > /dev/null') == 0:
                mpiexec = 'aprun -n ' + str(options.parallel)
            elif os.system('which mpirun.lsf > /dev/null') == 0:
                # mpirun.lsf does NOT need the number of processors (options.parallel).
                mpiexec = 'mpirun.lsf'
            else:
                sys.exit('Unable to execute parallel run.  Please edit the script to use your MPI run command, or run manually with something like: mpirun -np 4 ./cism_driver mismip3dInit.config')

            runstring = mpiexec + ' ' + options.executable + ' ' + configfile
            print( 'Executing parallel run with:  ' + runstring + '\n\n')

            # Here is where the parallel run is actually executed.
            os.system(runstring)

    print( 'Finished experiment', expt)

    # Change to parent directory and continue.
    os.chdir('..')
