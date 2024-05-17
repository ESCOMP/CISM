Slab test case
==============

Note: For setting up the experiments in an NCAR computing environment,
follow the steps in the README.NCAR_HPC file in the tests directory.

This directory contains python scripts for running an experiment involving a
uniform, infinite ice sheet ("slab") on an inclined plane.

The test case is described in sections 5.1-5.2 of:
    Dukowicz, J. K., 2012, Reformulating the full-Stokes ice sheet model for a
    more efficient computational solution. The Cryosphere, 6, 21-34,
    doi:10.5194/tc-6-21-2012.

Some CISM results from this test case are described in Sect. 3.4 of:
    Robinson, A., D. Goldberg, and W. H. Lipscomb, 2022, A comparison of the
    stability and performance of depth-integrated ice-dynamics solvers.
    The Cryosphere, 16, 689-709, doi:10.5194/tc-16-689-2022.

The test case consists of an ice slab of uniform thickness moving down an
inclined plane by a combination of sliding and shearing.
Analytic Stokes and first-order velocity solutions exist for all values of Glen's n >= 1.
The solutions for n = 1 are derived in Dukowicz (2012), and solutions for n > 1
are derived in an unpublished manuscript by Dukowicz (2013).

The original scripts, runSlab.py and plotSlab.py, were written by Matt Hoffman
with support for n = 1.  They came with warnings that the test is not supported.
The test is now supported, and the scripts include some new features:

* The user may specify any n >= 1 (not necessarily an integer).
  The tests assume which_ho_efvs = 2 (nonlinear viscosity) with flow_law = 0 (constant A).
* Physics parameters are no longer hard-coded.  The user can enter the ice thickness,
  beta, viscosity coefficient (mu_n), and slope angle (theta) on the command line.
* The user can specify time parameters dt (the dynamic time step) and nt (number of steps).
  The previous version did not support transient runs.
* The user can specify a small thickness perturbation dh, which is added to the initial
  uniform thickness via random sampling from a Gaussian distribution.
  The perturbation will grow or decay, depending on the solver stability for given dx and dt.

The run script is executed by a command like the following:

> python runSlab.py -n 4 -a DIVA -theta 0.0573 -thk 1000. -mu 1.e5 -beta 1000.

In this case, the user runs on 4 processors with the DIVA solver, a slope angle of 0.0573 degrees,
Glen's n = 1 (the default), slab thickness H = 1000 m, sliding coefficient beta = 1000 Pa (m/yr)^{-1},
and viscosity coefficient 1.e5 Pa yr.
These parameters correspond to the thick shearing test case described by Robinson et al. (2021).

To see the full set of command-line options, type 'python runSlab.py -h'.

Notes on effective viscosity:
   * For n = 1, the viscosity coefficient mu_1 has a default value of 1.e6 Pa yr in the relation
     mu = mu_1 * eps((1-n)/n), where eps is the effective strain rate.
   * For n > 1, the user can specify a coefficient mu_n; otherwise the run script computes mu_n
     such that the basal and surface speeds are nearly the same as for an n = 1 case with the
     mu_1 = 1.e6 Pa yr and the same values of thickness, beta, and theta.
   * There is a subtle difference between the Dukowicz and CISM definitions of the
     effective strain rate; the Dukowicz value is twice as large. Later, it might be helpful
     to make the Dukowicz convention consistent with CISM.

Run the plotting script, plotSlab.py, by typing 'python plotSlab.py'.  Two plots should appear.
The first plot shows the vertical velocity profile in nondimensional units and in units of m/yr.
There is excellent agreement between higher-order CISM solutions and the analytic solution
for small values of the slope angle theta.  For steep slopes, the answers diverge as expected.

For the second plot, the extent of the y-axis is wrong. This remains to be fixed.

This directory also includes a new script, stabilitySlab.py, to carry out the stability tests
described in Robinson et al. (2021).
For a given set of physics parameters and stress-balance approximation (DIVA, L1L2, etc.),
this script launches multiple CISM runs at a range of grid resolutions.
At each grid resolution, the script determines the maximum stable time step.
A run is deemed stable when the standard deviation of an initial small thickness perturbation
is reduced over the course of 100 time steps.  A run is unstable if the standard deviation
increases or if the model aborts (usually with a CFL violation).

To run the stability script, type a command like the following:

> python stabilitySlab.py -n 4 -a DIVA -theta 0.0375 -thk 1000. -mu 1.e5 -beta 1000.  \
  -dh 0.1 -nt 100 -nr 12 -rmin 10. -rmax 40000.

Here, the first few commands correspond to the thick shearing test case and are passed repeatedly
to the run script.  The remaining commands specify that each run will be initialized
with a Gaussian perturbation of amplitude 0.1 m and run for 100 timesteps.
The maximum stable timestep will be determined at 12 resolutions ranging from 10m to 40 km.
This test takes several minutes to complete on a Macbook Pro with 4 cores.

To see the full set of commmand line options, type 'python stabilitySlab.py -h'.

For questions, please contact William Lipscomb (lipscomb@ucar.edu) or Gunter Leguy (gunterl@ucar.edu).
