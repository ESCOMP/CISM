#!/usr/bin/env python

# This script runs a specific Shakti experiments. The name of the experiment
# is defined by the name of the folder where the experiment is staged.


import sys, os
from optparse import OptionParser
from configparser import ConfigParser


# Parse options.
optparser = OptionParser()

optparser.add_option('-e', '--exec',     dest='executable', type ='string',  default ='./cism_driver', help="Path to the CISM executable")
optparser.add_option('-x', '--expt',     dest='experiment', type ='string',  default = '50m', help="Shakti exp to run", metavar="EXPT")
optparser.add_option('-n', '--parallel', dest='parallel',   type ='int',     default = 1,     help="Number of processors: if specified then run in parallel", metavar="NUMPROCS")
optparser.add_option(      '--restart' , action = "store_true", dest='restart', default = False, help="Number of processors: if specified then run in parallel", metavar="NUMPROCS")


for option in optparser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = optparser.parse_args()


experiments = [options.experiment]


    
# Loop through experiments.
for expt in experiments:
    print('Running experiment', expt)

    # Change to directory for this experiment.
    os.chdir(expt)
    
    # Set config file to pass at command line when launching the run.
    configfiletowrite = 'hydro.config'
    # Read the master config file.
    config = ConfigParser()
    config.read(configfiletowrite)
    if options.restart:
        config.set('options', 'restart', str(1))
    else:
        config.set('options', 'restart', str(0))
    
    # Write to the master config file.
    with open(configfiletowrite, 'w') as configfile:
        config.write(configfile)
        
    sys.exit('let see')
    

    # Run CISM.

    print('parallel =', options.parallel)

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
                sys.exit('Unable to execute parallel run.  Please edit the script to use your MPI run command, or run manually with something like: mpirun -np 4 ./cism_driver mismip+Init.config')

            runstring = mpiexec + ' ' + options.executable + ' ' + configfile
            print('Executing parallel run with:  ' + runstring + '\n\n')

            # Here is where the parallel run is actually executed!
            os.system(runstring)

    print('Finished experiment', expt)

    # Change to parent directory and continue.
    os.chdir('..')
