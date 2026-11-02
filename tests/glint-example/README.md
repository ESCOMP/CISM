=============
Glint example
=============

INPUT DATA:
-----------

To run this glint-example test case, you will need to download the following
climate files:

* `ncep-doe_6h_climate.64x32.nc` (precip and 2-m temp on a global grid)

* `orog.igcmgrid.64x32.nc` (surface elevation on a global grid)

You will also need a 20-km Greenland input file:

* `gland20.input.nc`

These three netCDF files are available in the CESM inputdata repository (account
required), here:

<https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/glc/cism/gland20.input.nc>

and:

<https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/glc/cism/glint_example/>

To run the 5-km Greenland test case, you will also need the 5-km input file
`gland5.input.nc`, available in the CESM inputdata repo:

<https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/glc/cism/gland5.input.nc>

The 5-km test case and 20-km test case use the same climate files.

On Derecho, all of these files can be copied from the CESM inputdata directory
(no svn checkout needed):

```sh
cp /glade/campaign/cesm/cesmdata/inputdata/glc/cism/gland20.input.nc .
cp /glade/campaign/cesm/cesmdata/inputdata/glc/cism/gland5.input.nc .
cp /glade/campaign/cesm/cesmdata/inputdata/glc/cism/glint_example/*.nc .
```


BUILDING AND RUNNING THE TEST CASES:
------------------------------------

glint_example is run using the `cism_driver` executable with two input
arguments, e.g.:

```sh
./cism_driver ice_sheet.config climate.config
```

First, copy or link the executable to this test directory:

```sh
ln -s ../../builds/mac-gnu/cism_driver/cism_driver
```

(or `../../builds/derecho-intel/cism_driver/cism_driver` for an in-place build on Derecho).

As of October 2014, the repo contains three Greenland config files:
1. `greenland_20km.config.pdd`
2. `greenland_5km.config.pdd`
3. `greenland_20km.config.smb`

These files have settings that are physically reasonable but not scientifically
validated.

Files 1 and 2 are basically identical except for the grid resolution,
timestep and input file.  Both these files use the native Glint positive-degree-day (PDD)
scheme to force the Glide dycore, which uses the shallow-ice approximation.

Files 1 and 3 are identical except that file 3 has acab_mode = 0, which means
that Glint uses an externally computed surface mass balance (SMB) instead of
computing the SMB with a PDD scheme.  Normally the SMB would be computed and
passed in by a climate model such as CESM.  But for the standalone code, the
SMB is estimated crudely based on the input temperature and precip fields.
This allows us to test the Glint interfaces that would be called by a climate
model, without actually coupling to a climate model.

There are two glint_example config files in the repo:
1. `glint_example.config.pdd`
2. `glint_example.config.smb`

These files are set up to use the global NCEP climatological data for
downscaling to the Greenland grid. These files also set the length of the model
run (by default, we run a 10-year smoke test).  The only difference is that
file 1 should be run with the PDD scheme and file 2 with the SMB scheme.

The three standard diagnostic tests in this directory can be launched as
follows:

```sh
./cism_driver greenland_20km.config.pdd glint_example.config.pdd
./cism_driver greenland_5km.config.pdd glint_example.config.pdd 
./cism_driver greenland_20km.config.smb glint_example.config.smb
```

These are currently the only standard tests that exercise the Glint part of the code.

Note (November 2026): The 5-km test case is currently unstable with
`gland5.input.nc`: the initial velocity solve gives surface speeds of order
10^5 m/yr, and the run fails when writing output at year 10.
This was not investigated further; it may require a shorter timestep or other
config changes. The two 20-km test cases run correctly and are used to test
Glint.

Diagnostic output is written to a file called `greenland_20km.config.smb.log` or
something similar.  NetCDF output is written to various files at a frequency
set in the Greenland config file.
