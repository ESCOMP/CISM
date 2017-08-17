import os
import errno
import sys

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
        mpiexec = 'aprun -n ' + str(args.parallel)+" "
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
            print("   something like: mpirun -np 4 ./cism_driver dome.config")
            sys.exit(1)
    else:
        mpiexec = ''

    commands.append(mpiexec+driver+" "+config_name+quiet_mod)

    return commands
