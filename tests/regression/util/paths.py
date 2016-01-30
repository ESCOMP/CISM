import os
import sys
import errno
import fnmatch
import argparse
 
# small helper function so argparse will understand unsigned integers
def unsigned_int(x):
    x = int(x)
    if x < 0:
        raise argparse.ArgumentTypeError("This argument is an unsigned int type! Should be an integer greater than zero.")
    return x

#NOTE: this parser is used to find the directory names needed for the reg_test
#      tree structure. The DOF and PROCESSOR directory levels are determined by
#      the individual test options defined in the util.dicts module
run_parser = argparse.ArgumentParser(prog='RUN',add_help=False)
run_parser.add_argument('-n','--parallel', metavar='N', type=unsigned_int, default=0, 
        help="Run in parallel using N processors.")
run_parser.add_argument('--scale', type=unsigned_int, default=0, 
        help="Change the degrees of freedom.")


def recursive_glob(tree, pattern):
    matches = []
    for base, dirs, files in os.walk(tree):
        goodfiles = fnmatch.filter(files, pattern)
        matches.extend(os.path.join(base, f) for f in goodfiles)
    return matches

def file_modifier_list(args):
    mod_list = []
    #TODO: turn on args.dycore option.
    #if args.dycore: 
    #    dy_split = str.split(args.dycore,"-")
    #    mod_list.extend(dy_split[1:])
    if args.tmod:
        mod_list.append(args.tmod)

    return mod_list


def path_modifier_list(args):
    mod_list = []
    #TODO: turn on args.dycore option.
    #if args.dycore: 
    #    dy_split = str.split(args.dycore,"-")
    #    mod_list.append(dy_split[0])
    #TODO: turn on args.library option.
    #if args.library:
    #    mod_list.append(args.library)

    return mod_list


def make_absolute(args):
    args.out_dir = os.path.abspath(args.out_dir) 
    args.cism_dir = os.path.abspath(args.cism_dir) 
    args.build_dir = os.path.abspath(args.build_dir) 

    return args


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise


def mkdir_test(args, test_dict):
    """
    Sets up a regression testing data directory that looks like:
        reg_test
        |-- PLATFORM-COMPILER
            |-- ICE-MODEL
                |-- TEST-CASE
                    |-- DOF (degrees of freedom)
                        |-- PROCESSORS
                            |-- files.ext
    """

    mod_list = path_modifier_list(args)

    mod_dir = ''
    if mod_list:
        mod_dir = "-"+str.join("-", mod_list)
   
    #TODO: ice_model will be a indicator for which particular ice-sheet model was used.
    #      right now, this is a rather moot point as BATS only works for CISM-GLISADE
    ice_model = "CISM-glissade"
    data_dir = args.out_dir+os.sep+args.platform+'-'+args.compiler+mod_dir+os.sep+ice_model
    

    # make the output directories softly
    mkdir_p(args.out_dir)
    mkdir_p(data_dir)
   
    # clean up run files, if the exist
    # NOTE: this is needed because runnit.hpc uses recursive_globs over the created
    #       run files in order to determine what run files are new. With force, all 
    #       run files will appear to be old. 
    if args.force:
        all_run_files = recursive_glob(data_dir,"*.run")
        for rf in all_run_files:
            os.remove(rf)


    for case in test_dict:
        case_dir = os.path.normpath(data_dir+os.sep+str.split(str.split(case," ")[0],"/")[-1])
        run_script, mod_dict = test_dict[case]
        
        run_args = run_parser.parse_args(['--scale', '0', '-n', '1'])
        run_case_dir = case_dir+os.sep+'s'+str(run_args.scale)+os.sep+'p'+str(run_args.parallel)

        mkdir_p(run_case_dir)
        if not args.force and os.listdir(run_case_dir):
            print("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n") 
            print(  "WARNING: Test data already exists in:")
            print("\n"+run_case_dir+"\n")
            print(  "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
            print(  "Some data may be overwritten. Either specify a different test directory, or ")
            print(  "re-run with the -f or --force option to ignore existing data.")
            print("\nExiting...")
            sys.exit(1)
        
        if args.performance and mod_dict:
            for mod in mod_dict:
                run_args = run_parser.parse_args(mod_dict[mod].split())
                run_case_dir = case_dir+os.sep+'s'+str(run_args.scale)+os.sep+'p'+str(run_args.parallel)

                mkdir_p(run_case_dir)
                if not args.force and os.listdir(run_case_dir):
                    print("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n") 
                    print(  "WARNING: Test data already exists in:")
                    print("\n"+run_case_dir+"\n")
                    print(  "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
                    print(  "Some data may be overwritten. Either specify a different test directory, or ")
                    print(  "re-run with the -f or --force option to ignore existing data.")
                    print("\nExiting...")
                    sys.exit(1)


    return data_dir


def cmake(args):
    #TODO: turn on args.dycore option.
    #if args.dycore:
    #    cmake_dir = args.cism_dir+os.sep+'builds'+os.sep+args.platform+'-'+args.compiler+'-'+str.split(args.dycore,"-")[0]
    #    cmake_file = args.platform+'-'+args.compiler+'-'+args.dycore+'-cmake.bash'
    #else:
    #    cmake_dir = args.cism_dir+os.sep+'builds'+os.sep+args.platform+'-'+args.compiler
    #    cmake_file = args.platform+'-'+args.compiler+'-cmake.bash'
    cmake_dir = args.cism_dir+os.sep+'builds'+os.sep+args.platform+'-'+args.compiler
    cmake_file = args.platform+'-'+args.compiler+'-cmake.bash'

    

    if not args.skip_build and not os.path.isdir(cmake_dir):
        print("ERROR: cannot find your platform-compiler build directory: "+cmake_dir)
        print("   See the cmake builds directory for supported platform-compiler combinations.")
        sys.exit(1)

    if not args.skip_build and not os.path.isfile(cmake_dir+os.sep+cmake_file):
        print("ERROR: cannot find your cmake file: "+cmake_file)
        print("   Looked for the cmake file in: "+cmake_dir)
        sys.exit(1)
    
    return (cmake_dir, cmake_file)

