# ALMA-FETMS-beameff
Beam efficiency calculator for the ALMA Front End Test and Measurement System

* beameff2_CPP contains C++ source code for the 2.x versions.
* beameff_C contains the C source code and make files for the 1.x versions.
* LabVIEW contains the wrapper GUI source code and a LabVIEW 2013 executable.
* WinExe contains the latest Windows executable.
* doc contains user and reference documentation.
* iniparser is a Git submodule link to iniparser 4.0.

http://stackoverflow.com/questions/3796927/how-to-git-clone-including-submodules 

Recent version history:
 * 2.0.7:  Added input file option [settings] invertphase with options {lsb, usb, all, none}
           and synonymns {yes, y, t, 1} for all, {no, n, f, 0} for none.  All options case-insensitive.
 * 2.0.6:  Using multi-pass search with exact expression for fitted phase:
           1) CG search from {0, 0, probe-zdistance}
           2) Line search along delta_z
           3) CG search, line search, final CG search
           Fix calculation error in squint_percent
 * 2.0.5:  Using exact expression for fitted phase including corrections from Alvaro Gonzalez
           Invert phase and rotate for LSB now.
           Invert signs of nominal pointing angles.  Matches scanner coord system now.
           Added makePhaseFitPlot for visualizing the fitted phase surface.
           Improved documentation of FitPhase and ScanDataRaster::calcPhaseEfficiency
 * 2.0.4:  Don't rotate nearfield scans but do invert phase when appropriate.
           Back to el = +2.48 for band1test dewar
 * 2.0.3:  Adjust labels on pointing angles plot
 * 2.0.2:  Added makefile and minor fixes for compilation on Linux
 * 2.0.1:  Name pointing angle plots for RF if scanset_id not provided
           Name [results_xx] section for RF if scanset_id not provided
           Use el = -2.48 for band1test dewar
           Calculate squint even when 180-degree scan for correcting probe asymmetry is not provided.
           For backwards compatibility with LabVIEW wrapper, squint results in output section for pol1, copol.
 * 2.0.0:  Newly ported to C++ code in beameff2_CPP directory.
 * 1.3.6:  Using iniparser 4.0 from https://github.com/ndevilla/iniparser.  Prior was 1.8.
 * 1.3.5:  Added PointingOptionType instead of int ACA7Meter.  Added band 1 test dewar pointing.
 * 1.3.4:  Removed dead code and unused variables throughout.
