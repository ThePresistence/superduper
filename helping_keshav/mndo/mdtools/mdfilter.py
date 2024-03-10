#!/usr/bin/python
#
# --------------------------
#
# mdfilter.py
#
# Tom Keal, July 2007
#
# --------------------------
#
"""mdfilter - filter a list of sample geometries

 PURPOSE

 mdfilter takes a list of sample geometries in XYZ format and filters them
 according to a number of criteria:

 1. Successful mapping of active space
     - Each geometry is tested against a reference (template)
       (where the active space is assumed to be known for certain).
       If the active space cannot be mapped with a confidence over
       a certain threshold (default: 90%), it is rejected.

       Active space mapping can be switched off with the option --no-map

 2. Energy window
     - The excitation energy of each state is compared to that
       of the reference (or a user defined value). If it does not lie
       within a given tolerance (default: 0.15eV), it is rejected.

       This filtering can be surpressed with the option --no-window

 3. Transition probability
     - The transition probability is calculated (following Ref. Eq. 4) as
                 P = (f/DE^2) / max(f/DE^2)
       where f is the oscillator strength, and DE is the energy gap
       between states. For each state, a random number 0 < n < 1 is
       generated, and if n > P, the state is rejected.
       
       This filtering can be surpressed with the option --no-prob

       Ref: M. Barbatti et al., J. Photochem. Photobiol. A, v190, p228 (2007)

 In addition, filters (2) and (3) are limited by default to excitations
 to a single target excited state (default: 2). This restriction can be
 removed with the option --all-states

 The default output is a list of XYZ geometries, with details in the title
 line of the selected state and mapped active orbitals.

 Additionally, MNDO input files can be generated for subsequent excited state
 dynamics runs (using the option --dynamics), and/or an absorption
 spectrum can be simulated (using the option --spectrum)

 The filtered geometries can be used as the basis of an excited states
 dynamics run or simply to simulate an absorption spectrum.


 REQUIREMENTS

 Requires the mndotools.py library


 USAGE

 mdfilter.py [OPTIONS]


 GENERAL INPUTS (default file names)

 template.inp - MNDO template input file, detailing CI options,
                reference geometry and active space, etc.
                Usually would be a ground state optimised geometry or
                equilibrated geometry of some kind.
 sample.xyz   - List of sampled geometries in XYZ format

 INPUT FOR DYNAMICS PREPARATION (option --dynamics)

 dynvar.in    - Template file for dynamics input 
                (this will be left unchanged except possibly for 
                 init_stat and restart)

 INPUT FOR REFILTERING (option --refilter)

 filtered-info.log       - Output from a previous mdfilter run. Use the raw
                           data from this log file instead of recalculating
                           the data with an MNDO calculation. The refiltering
                           input must have been produced using the same
                           template molecule and sample XYZ geometries.

 GENERAL OUTPUTS (default file names)

 (standard output)       - Basic progress information
 filtered-sample.xyz     - List of geometries that passed through the filter.
 filtered-info.log       - Summary information and further information for
                           each geometry
                           (active space, excitation energy, probability).

 OUTPUT FOR DYNAMICS PREPARATION (option --dynamics)

 For each trajectory (here the first),

 run0001/mndo.inp        - Input file with geometry and standard MNDO options
                           (based on template.inp)
 run0001/dynvar.in       - Input file with dynamics options
                           (unchanged except for init_stat/restart)

 OUTPUT FOR SPECTRUM SIMULATION (option --spectrum)
                           
 filtered-spectrum.dat   - Histogram data (default bin width=0.1eV) for
                           simulated absorption spectra.


 DEFAULTS 

 General calculation options
   - Path to MNDO executable: the output of 'which mndo-md'
   - Target excited state: 2
   - Oscillator strength formalism: rp

 Active space filtering
   - Threshold for successful mapping: 0.9

 Energy window filtering
   - Target excitation energy: that of the template file (state 2)
   - Tolerance for energy window filtering: 0.15eV 

 Dynamics preparation
   - ncigrd gradients: from ground state up to target-state
                       (or up to iroot with --all-states)

 Absorption spectrum
   - Histogram bin size for absorption spectrum: 0.1eV

 All defaults can be changed using the command line options below.


 GENERAL OPTIONS

 -h, --help
        display this help message
        
 -q, --quiet
        don't print anything to standard output

 --mndo-path=(path string)
        specify a new path for MNDO. Default is the result of 'which mndo-md'

 -r, --refilter
        do not rerun sample MNDO calculations, but read quantities directly
        from the mdfilter log file from a previous run (must correspond to the
        same input template and xyz files)

 OPTIONS FOR STANDARD INPUT FILE NAMES

 --template=(file name)
        specify a different file name for the template molecule

 --sample=(file name)
        specify a different file name for the sample geometry list

 OPTIONS FOR FILTERING BY STATE

 --target-state=(int)
        specify a new target excited state (integer).
        By default all states are rejected except the target state
        This option is also used by energy window filtering (see below)
        
 -a, --all-states
        consider possibility of excitation to any state
        (i.e. do not filter by state)

 OPTIONS FOR FILTERING BY ACTIVE SPACE MAPPING

 --mapthr=(float)
        specify a new threshold for filtering by active orbital mapping (0-1)
        This option is also used when preparing dynamics inputs (see below)

 -M, --no-map
        do not attempt to map the active molecular orbitals
        (i.e. imomap=0). With this option values of movo != 1 can be used.

 OPTIONS FOR FILTERING BY ENERGY WINDOW

 --target-state=(int)
        specify a new target excited state (integer)
        The excitation energy from the template molecule for the given state
        will be used as the midpoint of the energy window.

 --window-target=(real)
        specify a target excitation energy for the energy window
        (instead of taking the value from the template molecule calculation) 

 --window-tolerance=(real)
        specify a new tolerance for energy window filtering

 -W, --no-window
        do not reject structures whose excitation energies fall outside the
        energy window. 

 OPTIONS FOR FILTERING BY TRANSITION PROBABILITY

 --f-formalism=r|p|rp
        MNDO provides oscillator strengths calculated using three different
        formalisms (which are identical only in the limit of an exact
        wavefunction).
        This option specifies which formalism should be used in calculating
        the transition probability.
        (note: cannot be changed when using --refilter)

 --random-seed=(int)
        specify a seed for the random numbers used by the transition
        probability filter

 -P, --no-prob
        do not reject structures which have been rejected according to the
        stochastic algorithm based on transition probability. 

 OPTIONS FOR REFILTERING

 -r, --refilter
        do not rerun sample MNDO calculations, but read quantities directly
        from the mdfilter log file from a previous run (must correspond to the
        same input template and xyz files)

 --refilter-info=(file name)
        specify a new name for the refilter info input file
 
 OPTIONS FOR STANDARD OUTPUT FILE NAMES

 --filtered-info=(file name)
        specify a new name for the log file
 
 --filtered-sample=(file name)
        specify a new name for the filtered list of sample geometries 

 OPTIONS FOR PREPARATION OF DYNAMICS RUN FILES

 -d, --dynamics
        prepare excited state dynamics run files

 --md-states=(string)
        specify the gradients to be calculated during the dynamics runs.
        By default, gradients will be calculated for the ground state up to
        the target state (or 'iroot' in the case of --all-states).
        With this option a specific set of states can be given
        (separated by spaces, e.g. --md-states='1 2 3' for the first 3 states)
        Note this may result in dynamics runs that crash, if init_stat is
        not included in --md-states.

 --mapthr=(float)
        specify a new threshold for active orbital mapping during the dynamics
        run (0-1, converted to a percentage in the MNDO input)
        This option is also used for active space filtering (see above)

 -M, --no-map
        if filtering by active space mapping has been switched off, 
        mapping will also be off during the dynamics run (imomap=0) 

 --dynvar=(file name)
        specify a new file name for the dynamics options input file

 --md-prefix=(string)
        specify a new prefix for the run file directories

 --vel
        provide initial velocities for dynamics runs from a sample velocity 
        file. In this case a minimal restart file will be created for each
        run containing the appropriate initial velocities.

 --vel-inp=(file name)
        specify a new file name for the sample velocity input file

 OPTIONS FOR SIMULATION OF ABSORPTION SPECTRA

 -s, --spectrum
        write an absorption spectrum file

 --bin-size=(float)
        specify a new histogram bin size for the absorption spectrum (in eV)

 --filtered-spectrum=(file name)
        specify a new name for the histogram absorption spectrum file
        
"""
import sys
import getopt
import os
import commands
import random
import math
import mndotools
import subprocess
import time


def main():
    """Main routine for filtering.

    1. Parse command line options and set internal options
    2. Read in and run the template file MNDO job
    3. Read in the sample XYZ geometries

    Either:
     4a. Create MNDO input files for each sample geometry, run the job and
         extract mapping information and transition properties from the outputs
    Or:
     4b. Read in the same information from an old mdfilter log file
         (with the option --refilter)

    5. Filter by state, mapping, energy window and transition probability
       in that order
    6. Write out the filtered XYZ geometries and log file
    7. Prepare dynamics input files
    8. Write histogram data for an absorption spectrum

    """
    # Set default global options and parse command line options
    options = GlobalOptions()
    options.parseCommandLine()

    # Read in the template file and calculate
    if options.writeStdOut:
        sys.stdout.write("Reading in template file...\n")
    templateMol = TemplateMol(options.fileInpTemplate, options.filterMapthr)
    if options.writeStdOut:
        sys.stdout.write("Running template job...\n")
    templateMol.createInput(options.filterMapthr)
    templateMol.calc(options.mndoPath)
    templateMol.extractQuantities(options.fFormalism)
    for state in templateMol.excitedStates:
        sys.stdout.write("%4s   %10.6f  %10.6f\n" % (
            state, templateMol.excitedStateEnergies[state],
            templateMol.oscillatorStrengths[state]))
      
    # Delete output from memory - no longer needed
    templateMol.output = None

    # Set the energy window if appropriate
    if options.filterWindow and options.windowTemplateAsRef:
        try:
            for i in templateMol.excitedStates:
                options.windowTarget.append\
                    (templateMol.excitedStateEnergies[i])
        except KeyError:
            sys.stderr.write("Error: target state is not one of the "
                             "calculated states.\n")
            sys.exit(1)

    # Read in the sample geometries
    if options.writeStdOut:
        sys.stdout.write("Reading in sample geometries file...\n")
    sampleXYZ = mndotools.XYZFile(options.fileInpSample)
    # sampleMols is a list of objects of type SampleMol
    sampleMols = []
    for i in range(sampleXYZ.noOfGeoms):
        sampleMols.append(SampleMol(sampleXYZ, i))

    # Read in velocity file if necessary
    if options.writeVelocities:
        if options.writeStdOut:
            sys.stdout.write("Reading in sample velocities file...\n")
        velocityFile = mndotools.VelocityFile(options.fileInpVelocities)

    # For each sample, create an MNDO job, run it, and extract information
    # (unless refiltering)
    if not options.refilter:
        if options.writeStdOut:
            sys.stdout.write("Running sample calculations... ")
        for i in range(sampleXYZ.noOfGeoms):
            if options.writeStdOut:
                sys.stdout.write("%s \n" % i)
                sys.stdout.flush()
            sampleMols[i].createInput(templateMol, options.filterMapthr)
            sampleMols[i].calc(options.mndoPath)
            sampleMols[i].extractQuantities(options.fFormalism, options.filterMapthr)
            # Delete output from memory - no longer needed
            sys.stdout.write(" State     Energy      Ediff       OscStr       Ratio\n")
            for state in templateMol.excitedStates:
              sys.stdout.write("%4s   %10.6f  %10.6f  %10.6f  %10.6f  \n" % (
                    state, sampleMols[i].excitedStateEnergies[state],
                    sampleMols[i].excitedStateEnergies[state]-
                    templateMol.excitedStateEnergies[state],
                    sampleMols[i].oscillatorStrengths[state],
                    sampleMols[i].ose2Ratio[state]))
            sampleMols[i].output = None
            if options.writeStdOut:
                sys.stdout.write("\n")
    else:
        # Refilter - read in sampleMols quantities from an old log file
        if options.writeStdOut:
            sys.stdout.write("Refilter: reading in log information...\n")
        refilterFile = FilterLogFile(options.fileInpRefilter)
        if refilterFile.fFormalism != options.fFormalism:
            sys.stderr.write("Cannot change f-formalism with refiltering.\n")
            sys.exit(1)
        if refilterFile.noOfGeoms != sampleXYZ.noOfGeoms:
            sys.stderr.write("Inconsistent number of geometries in refilter "
                             "log file.\n")
            sys.exit(1)
        for i in range(sampleXYZ.noOfGeoms):
            refilterFile.extractLogQuantities(i, sampleMols[i])

    # Initialise filter, and filter the sample geometries by state
    if options.writeStdOut:
        sys.stdout.write("Filtering by state...\n")
    for i in range(sampleXYZ.noOfGeoms):
        sampleMols[i].initFilter(options.allStates, options.targetState)

    # Filter by active space mapping
    if options.filterMapthr:
        if options.writeStdOut:
            sys.stdout.write("Filtering by active space mapping...\n")
        for i in range(sampleXYZ.noOfGeoms):
            sampleMols[i].filterMapping(options.mapThreshold)

    # Filter by energy window
    if options.filterWindow:
        if options.writeStdOut:
            sys.stdout.write("Filtering by energy window...\n")
        for i in range(sampleXYZ.noOfGeoms):
            sampleMols[i].filterWindow(options.windowTarget,
                                       options.windowTolerance)

    # Filter by transition probability
    maxRatio = 0.0
    if options.writeStdOut and options.filterProb:
	sys.stdout.write("Filtering by transition probability...\n")
    # First we find the maximum value of the ratio for transition
    # to any state
    # We assume here that each sample molecule has the same number of
    # states as the template molecule (reasonable as that's why it's
    # a template molecule!)
    for i in range(sampleXYZ.noOfGeoms):
	for state in templateMol.excitedStates:
	    if sampleMols[i].selected[state]:
		if sampleMols[i].ose2Ratio[state] > maxRatio:
		    maxRatio = sampleMols[i].ose2Ratio[state]
    # Now we calculate the probability for each state and
    # filter by stochastic algorithm
    # First we set the seed
    if options.probUseSeed:
	random.seed(options.probSeed)
    else:
	random.seed()
    for i in range(sampleXYZ.noOfGeoms):
	sampleMols[i].filterProb(maxRatio,options.filterProb)

    # Summarise the filtering
    sumState = {}
    sumMap = {}
    sumWindow = {}
    sumProb = {}
    sumPass = {}
    for state in templateMol.excitedStates:
        sumState[state] = 0
        sumMap[state] = 0
        sumWindow[state] = 0
        sumProb[state] = 0
        sumPass[state] = 0
    for i in range(sampleXYZ.noOfGeoms):
        for state in templateMol.excitedStates:
            if sampleMols[i].selected[state]:
                sumPass[state] += 1
            elif sampleMols[i].failCause[state] == 'state':
                sumState[state] += 1
            elif sampleMols[i].failCause[state] == 'map':
                sumMap[state] += 1
            elif sampleMols[i].failCause[state] == 'window':
                sumWindow[state] += 1
            elif sampleMols[i].failCause[state][:4] == 'prob':
                sumProb[state] += 1                
            else:
                sys.stderr.write("Warning: failure cause of molecule %s"
                                 " (%s) not recognised.\n" % (i,
                                 sampleMols[i].failCause[state]))
    
    # Write XYZ output
    if options.writeStdOut:
        sys.stdout.write("Writing filtered sample file...\n")
    try:
        spectrumFile = open(options.fileOutSample, 'w')
    except IOError:
        sys.stderr.write("Error: could not open %s" +
                         "for writing.\n" % options.fileOutSample)
    for i in range(sampleXYZ.noOfGeoms):
        for state in sampleMols[i].excitedStates:
            if sampleMols[i].selected[state]:
                # Write geometry
                title = "G" + str(i) + " S" + str(state) 
                if options.filterMapthr:
                        title = title + " O" + str(sampleMols[i].activeOrbitals)
                xyzGeom = sampleMols[i].xyzGeom.writeXYZGeom(title, False)
                for j in range(len(xyzGeom)):
                    spectrumFile.write(xyzGeom[j] + "\n")
    spectrumFile.close()                

    # Write log output
    if options.writeInfo:
        if options.writeStdOut:
            sys.stdout.write("Writing log file...\n")
        try:
            logFile = open(options.fileOutInfo, 'w')
        except IOError:            
            sys.stderr.write("Error: could not open %s" +
                             "for writing.\n" % options.fileOutInfo)
        # Command line options
        logFile.write("COMMAND LINE OPTIONS\n")
        for i in sys.argv[1:]:
            logFile.write(i + " ")
        logFile.write("\n\n")
        # Internal options
        logFile.write("INTERNAL OPTIONS\n")
        intOptions = options.__dict__.keys()
        intOptions.sort()
        for i in intOptions:
            logFile.write(i + ": " + str(options.__dict__[i]) + "\n")
        logFile.write("\n")
        # Summary of filtering
        logFile.write("SUMMARY OF FILTERING\n")
        logFile.write("Geometries: %s\n" % sampleXYZ.noOfGeoms)
        logFile.write("States: %s\n\n" % templateMol.noOfStates)
        logFile.write(" State   F(State)      F(Map)   F(Window)     "
                      "F(Prob)      PASSED\n")
        for i in templateMol.excitedStates:
            causes = {}
            logFile.write("%4s   %10s  %10s  %10s  %10s  %10s\n" % (i,
                                                              sumState[i],
                                                              sumMap[i],
                                                              sumWindow[i],
                                                              sumProb[i],
                                                              sumPass[i])) 
        logFile.write("\n")
        logFile.write("Maximum f(KO)/dE^2(KO) ratio:  %10.6f\n" % maxRatio)
        logFile.write("\n")
        # Individual geometries
        logFile.write("INDIVIDUAL GEOMETRIES\n")
        for i in range(sampleXYZ.noOfGeoms):
            logFile.write("Geom: %s\n" % i)
            if options.filterMapthr:
                logFile.write(" AcOrb    Mapping\n")
                for orb in range(len(sampleMols[i].activeOrbitals)):
                    logFile.write("%4s      %7.3f\n" % (
                            sampleMols[i].activeOrbitals[orb],
                            sampleMols[i].mapOverlap[orb]))
            logFile.write(" State     Energy      OscStr       Ratio   "
                          "TransProb      FCause\n")
            for state in templateMol.excitedStates:
                logFile.write("%4s   %10.6f  %10.6f  %10.6f  %10.6f  " % (
                    state, sampleMols[i].excitedStateEnergies[state],
                    sampleMols[i].oscillatorStrengths[state],
                    sampleMols[i].ose2Ratio[state],
                    sampleMols[i].transProb[state]))
                if sampleMols[i].selected[state]:
                    logFile.write("    PASSED\n")
                else:
                    logFile.write("%10s\n" % sampleMols[i].failCause[state])
            logFile.write("\n")
        logFile.close()

    # Prepare dynamics input files
    if options.writeDynamics:
        if options.writeStdOut:
            sys.stdout.write("Writing dynamics input files... ")
        # Read and parse dynvar.in
        dynvar = mndotools.DynvarFile(options.fileInpDynvar)
        (dynvarKeyOrder, dynvarKeyVals) = dynvar.parseNamelist(True)
        # Add an init_stat entry if there isn't one already
        if not dynvarKeyVals.has_key('INIT_STAT'):
            dynvarKeyOrder.append('INIT_STAT')
        # Add a restart entry if there isn't one already
        if not dynvarKeyVals.has_key('RESTART'):
            dynvarKeyOrder.append('RESTART')
        # Overwrite init_stat even it is already present, so we know
        # that there will be a crash if it isn't set properly
        dynvarKeyVals['INIT_STAT'] = -1
        # Set restart according to whether initial velocities are to be 
        # provided
        if options.writeVelocities:
            dynvarKeyVals['RESTART'] = 'T'
        else:
            dynvarKeyVals['RESTART'] = 'F'
        # For each selected state, create a directory with dynamics
        # input files.
        runCounter = 0
        for i in range(sampleXYZ.noOfGeoms):
            for state in sampleMols[i].excitedStates:
                if sampleMols[i].selected[state]:        
                    dynvarKeyVals['INIT_STAT'] = state
                    # Create a directory
                    runCounter += 1
                    if options.writeStdOut:
                        sys.stdout.write("%s " % runCounter)
                        sys.stdout.flush()
                    dirname = options.fileOutMDPrefix + '%04d' % runCounter
                    if os.path.exists(dirname):
                        sys.stderr.write("Warning: directory %s already "
                                         "exists - skipping.\n" % dirname)
                        continue
                    try:
                        os.mkdir(dirname)
                    except OSError:
                        sys.stderr.write("Warning: could not create "
                                  "directory %s - skipping.\n" % dirname)
                        continue
                    # Write dynvar.in
                    try:
                        dynvarFileName = '%s/dynvar.in' % dirname
                        df = open(dynvarFileName, 'w')
                    except IOError:
                        sys.stderr.write("Error: could not open %s" +
                                  "for writing.\n" % dynvarFileName)
                        continue
                    df.write(' &DYNVAR\n')
                    for key in dynvarKeyOrder:
                        df.write(' %s = %s,\n' % (key, dynvarKeyVals[key]))
                    df.write(' /\n')
                    df.close()
                    # Write mndo input
                    if len(options.mdUserStates) > 0:
                        gradients = options.mdUserStates
                    elif options.allStates:
                        gradients = range(1, templateMol.excitedStates[-1] + 1)
                    else:
                        gradients = range(1, options.targetState + 1)
                    inpFileName = '%s/mndo.inp' % dirname
                    title = "Geom %i, State %i\n" % (i, state) + \
                            "Generated automatically by mdfilter.py\n"
                    sampleMols[i].createDynamicsRunInput(inpFileName,
                                                         templateMol,
                                                         title,
                                                         gradients,
                                                         options.mapThreshold,
                                                         options.filterMapthr)
                    # Write dynam.restart file if providing initial velocities
                    if options.writeVelocities:
                        try:
                            restartFileName = '%s/dynam.restart' % dirname
                            df = open(restartFileName, 'w')
                        except IOError:
                            sys.stderr.write("Error: could not open %s" +
                                      "for writing.\n" % restartFileName)
                            continue
                        tmpGeom = sampleMols[i].xyzGeom
                        df.write('#restart   MD_only\n')
                        df.write('#N_atoms    %s\n' % tmpGeom.noOfAtoms)
                        df.write('#coordinates\n')
                        for atm in range(tmpGeom.noOfAtoms):
                            geomLine = "  %20.10f %20.10f %20.10f\n" % (
                                tmpGeom.atomCoords[atm].x,
                                tmpGeom.atomCoords[atm].y,
                                tmpGeom.atomCoords[atm].z)
                            df.write(geomLine)
                        df.write('#velocity\n')
                        tmpVel = velocityFile.getVelocity(i)
                        for atm in range(tmpGeom.noOfAtoms):
                            velLine = "  %20.10f %20.10f %20.10f\n" % (
                                tmpVel.atomVelocities[atm].x,
                                tmpVel.atomVelocities[atm].y,
                                tmpVel.atomVelocities[atm].z)
                            df.write(velLine)
                        df.close()
        if options.writeStdOut:
            sys.stdout.write("\n")

    # Write spectrum
    if options.writeSpectrum:
        if options.writeStdOut:
            sys.stdout.write("Writing spectrum file...\n")
        for state in sampleMols[1].excitedStates:
	    histogram = {}
	    for i in range(sampleXYZ.noOfGeoms):
		if sampleMols[i].selected[state]:
		    # Bins are stored as the integer multiple of the bin
		    # size to avoid floating point issues
		    energy = sampleMols[i].excitedStateEnergies[state]
		    count = int(round(energy / options.spectrumBinSize))
		    if histogram.has_key(count):
			histogram[count] += sampleMols[i].transProb[state]
		    else:
			histogram[count] = sampleMols[i].transProb[state]
	    # Print out the data from the lowest bin to the highest
	    if len(histogram)>0:
                bins = histogram.keys()
                bins.sort()
                # Put in dummy zero bins at both ends
                lowest = bins[0] - 1
                if lowest >= 0:
                    histogram[lowest] = 0
                else:
                    lowest = bins[0]
                highest = bins[-1] + 1
                histogram[highest] = 0
                fileName = "%s_%1d" % (options.fileOutSpectrum,state)
                try:
                    spectrumFile = open(fileName, 'w')
                except IOError:
                    sys.stderr.write("Error: could not open %s" +
                                    "for writing.\n" % fileName)
                for i in range(lowest, highest + 1):
                    realBin = i * options.spectrumBinSize
                    wavelength = 1239.8419/realBin
                    if histogram.has_key(i):
                        spectrumFile.write("%10.5f %10.5f %10.5f\n"% (realBin, wavelength, histogram[i]))
                    else:
                        spectrumFile.write("%10.5f %10.5f %10.5f\n"% (realBin, wavelength, 0.0))
                spectrumFile.close()

    # Remove tmp files if present
    if os.path.exists('tmp.inp'):
        os.remove('tmp.inp')
    if os.path.exists('tmp.out'):
        os.remove('tmp.out')


class GlobalOptions:
    
    """Set and store global options."""
    
    def __init__(self):
        """Set default options."""
        # MNDO path
        self.mndoPath = None
        p = commands.getstatusoutput('which mndo-md')
        if p[0] == 0:
            self.mndoPath = p[1]
        # Set default options
        # Program features
        self.refilter = False
        self.filterMapthr = True
        self.filterWindow = True
        self.filterProb = True
        self.writeStdOut = True
        self.writeInfo = True
        self.writeDynamics = False
        self.writeVelocities = False
        self.writeSpectrum = False
        # File I/O
        self.fileInpRefilter = None
        self.fileInpTemplate = "template.inp"
        self.fileInpSample = "sample.xyz"
        self.fileInpDynvar = "dynvar.in"
        self.fileInpVelocities = "sample.vel"
        self.fileOutSample = "filtered-sample.xyz"
        self.fileOutInfo = "filtered-info.log"
        self.fileOutMDPrefix = "run"
        self.fileOutSpectrum = "filtered-spectrum.dat"
        # General calculation options
        self.allStates = False
        self.targetState = 2
        self.fFormalism = 2  # 0 = r, 1 = p, 2 = rp
        # Specific filter options
        self.mapThreshold = 0.9
        self.windowTemplateAsRef = True
        self.windowTarget = []
        self.windowTolerance = 0.15 # eV
        self.probUseSeed = False
        self.probSeed = None
        # Dynamics options
        self.mdUserStates = []
        # Spectrum options
        self.spectrumBinSize = 0.1 # eV

    def parseCommandLine(self):
        """Parse command line options."""
        shortOptions = "adhMPqrsW"
        longOptions = ["all-states",
                       "help", "no-map", "no-window", "no-prob",
                       "quiet", "dynamics", "vel", "spectrum",
                       "bin-size=", "dynvar=", "vel-inp=", "f-formalism=",
                       "filtered-info=", "filtered-sample=",
                       "filtered-spectrum=", "mapthr=",
                       "md-prefix=", "md-states=", "mndo-path=",
                       "random-seed=", "refilter", "refilter-info=",
                       "sample=", "target-state=",
                       "template=", "tolerance=", "window-state=",
                       "window-target=", "window-tolerance="]
        try:
            opts, args = getopt.gnu_getopt(sys.argv[1:],
                                           shortOptions, longOptions)
        except getopt.error, msg:
            sys.stderr.write("%s\n" % msg)
            sys.stderr.write("for help use --help\n")
            sys.exit(2)
        # Process options
        for o, a in opts:
            # Help
            if o in ("-h", "--help"):
                print __doc__
                sys.exit(0)
            # Suppress or switch on features
            if o in ("-r", "--refilter"):
                self.refilter = True
                if self.fileInpRefilter is None:
                    self.fileInpRefilter = 'filtered-info.log'
            if o in ("-a", "--all-states"):
                self.allStates = True
            if o in ("-M", "--no-map"):
                self.filterMapthr = False
            if o in ("-W", "--no-window"):
                self.filterWindow = False
            if o in ("-P", "--no-prob"):
                self.filterProb = False
            if o in ("-q", "--quiet"):
                self.writeStdOut = False
            if o in ("-d", "--dynamics"):
                self.writeDynamics = True
            if o in ("--vel"):
                self.writeVelocities = True
            if o in ("-s", "--spectrum"):
                self.writeSpectrum = True
            # Change default files
            if o == "--refilter-info":
                self.fileInpRefilter = a
            if o == "--template":
                self.fileInpTemplate = a
            if o == "--sample":
                self.fileInpSample = a
            if o == "--dynvar":
                self.fileInpDynvar = a
            if o == "--vel-inp":
                self.fileInpVelocities = a
            if o == "--filtered-sample":
                self.fileOutSample = a
            if o == "--filtered-info":
                self.fileOutInfo = a
            if o == "--md-prefix":
                self.fileOutMDPrefix = a
            if o == "--filtered-spectrum":
                self.fileOutSpectrum = a
            # Change default parameters
            if o == "--mndo-path":
                self.mndoPath = a
            if o == "--target-state":
                try:
                    self.targetState = int(a)
                except ValueError:
                    sys.stderr.write("Error: state should be an integer.\n")
                    sys.exit(2)
            if o == "--f-formalism":
                if a == "r":
                    self.fFormalism = 0
                elif a == "p":
                    self.fFormalism = 1
                elif a == "rp":
                    self.fFormalism = 2
                else:
                    sys.stderr.write("Error: f-formalism should be 'r', 'p' "
                                     "or 'rp'.\n")
                    sys.exit(2)
            if o == "--mapthr":
                try:
                    self.mapThreshold = float(a)
                except ValueError:
                    sys.stderr.write("Error: mapthr must be a real number.\n")
                    sys.exit(2)                    
            if o == "--window-target":
                self.windowTemplateAsRef = False
                try:
                    self.windowTarget = {float(a)}
                except ValueError:
                    sys.stderr.write("Error: window-target must be a real "
                                     "number.\n")
                    sys.exit(2)
            if o == "--window-tolerance":
                try:
                    self.windowTolerance = float(a)
                except ValueError:
                    sys.stderr.write("Error: window-tolerance must supply a "
                                     "real number.\n")
                    sys.exit(2)
            if o == "--random-seed":
                self.probUseSeed = True
                try:
                    self.probSeed = int(a)            
                except ValueError:
                    sys.stderr.write("Error: random-seed must be an integer."
                                     "\n")
                    sys.exit(2)
            if o == "--md-states":
                for state in a.split():
                    try:
                        self.mdUserStates.append(int(state))
                    except ValueError:
                        sys.stderr.write("Error: md-states must be a space-"
                                     "separated list of integers\n")
                        sys.exit(2)
            if o == "--bin-size":
                try:
                    self.spectrumBinSize = float(a)
                except ValueError:
                    sys.stderr.write("Error: bin-size must supply a "
                                     "real number.\n")
                    sys.exit(2)
        # There should be no remaining arguments
        for arg in args:
            sys.stderr.write("Warning: argument '%s' not recognised.\n" % arg)
        # Check there is a valid MNDO path
        if self.mndoPath is None:
            sys.stderr.write("Error: no result from 'which mndo-md' and no "
                             "MNDO path specified.\n")
            sys.exit(2)
        # Check that the input files exist
        if not os.path.exists(self.fileInpTemplate):
            sys.stderr.write("Error: template file %s does not exist.\n"
                             % self.fileInpTemplate)
            sys.exit(2)
        if not os.path.exists(self.fileInpSample):
            sys.stderr.write("Error: input geometry file %s does not exist.\n"
                             % self.fileInpSample)
            sys.exit(2)
        if self.refilter and not os.path.exists(self.fileInpRefilter):
            sys.stderr.write("Error: refilter input file %s does not exist.\n"
                             % self.fileInpRefilter)
            sys.exit(2)            
        if self.writeDynamics and not os.path.exists(self.fileInpDynvar):
            sys.stderr.write("Error: dynamics input file %s does not exist.\n"
                             % self.fileInpDynvar)
            sys.exit(2)            


class Mol:

    """Base class for TemplateMol and SampleMol."""

    def __init__(self):
        pass

    def calc(self, mndoPath):
        """Run the molecule's input file as an MNDO job.

        mndoPath - path to the MNDO executable

        The output is redirected to the file tmp.out, and then read in as
        self.output

        Exit if the MNDO job failed
           (probably should do something less drastic!)

        """
        # The shell command to run the MNDO job.
        # We send the output to tmp so it's
        # easier to see what's gone wrong if things go wrong
        cmd = str(mndoPath) + " < " + self.fileName + " > tmp.out"
        c = commands.getstatusoutput(cmd)
        # Unfortunately a failed MNDO run usually (or possibly always) give a
        # return value of 0
        # so I'm not sure we will ever see this error
        if c[0] != 0:
            sys.stderr.write("Error: MNDO job did not run successfully.\n")
            sys.exit(1)
        # Now let's read the output in and check that it really has worked
        self.output = mndotools.OutputFile('tmp.out')
        # The following seems to be a reasonable thing to search for
        # But it does mean we need mprint>=0
        searchString = 'GRADIENTS FOR GEOMETRY OPTIMIZATION'
        lineNo = self.output.searchForString(searchString, False)
        if lineNo == -1:
            sys.stderr.write("Error: MNDO job appears to have failed.\n")
            sys.stderr.write("Could not find '%s'\n" % searchString)
            sys.exit(1)


class TemplateMol(Mol):
    
    """Structure from template file."""
    
    def __init__(self, fileInpTemplate, map=True):
        """Load in the template file and check that it is set up correctly."""
        Mol.__init__(self)
        # Dodgy hack for new style input creation where the runnable input
        # file is different to the template input file ...
        self.fileName = 'tmp.inp'
        self.inputFile = mndotools.InputFile(fileInpTemplate)
        (self.keyOrder, self.keyVals) = self.inputFile.parseKeywords(True)
        self.output = None
        sane = self.checkSanity(map)
        if not sane:
            # It's safer to ask for a rewrite than to silently correct
            # an insane input file
            sys.stderr.write("Error: rewrite your template file with "
                             "the correct mandatory options.\n")
            sys.exit(1)

    def checkSanity(self, map=True):
        """Check that correct keywords have been given in the template file."""
        # It's polite to assume sanity until proven otherwise
        sane = True
        # jop=-2, imomap=1, movo=1 removed - added in later as required
        mandatoryVals = { 'igeom':1,
                          'iform':1,
                          'kci':5,
                          'mciref':0 }
        if map:
            mandatoryVals['movo'] = 1
        # TODO: Handle mciref when != 0
        # ioutci>=1, iuvcd>=2 removed - added in later as required
        mandatoryMinVals = { 'iroot':2 }
        # This is not an exhaustive list of options that shouldn't be
        # present. In particular checking for jop=-2 and no ief makes
        # most optimization options irrelevant anyway.
        # Most of the options below are excluded because they require
        # extra input (don't want to have to deal with extra lines!).
        # Note some options (icross, keepci, etc.) are added in later
        # when they are needed.
        mandatoryAbsent = [ 'nexmol', 'mplib', 'ief', 'inrefd', 'iparok',
                            'mminp', 'icosmo', 'middle', 'ingeom', 'inpfrg',
                            'inp21', 'inp22', 'inp24', 'inp25',
                            'kmass', 'nmrnuc', 'kgeom',
                            'ksym', 'lroot', 'icimap', 'keepci', 'ncigrd',
                            'icross', 'jop', 'mapthr' ]
        for testKey in mandatoryVals.keys():
            if not self.keyVals.has_key(testKey):
                sys.stderr.write("Error: mandatory option %s missing from "
                                 "template file.\n" % testKey)
                sane = False
            elif int(self.keyVals[testKey]) != mandatoryVals[testKey]:
                sys.stderr.write("Error: value for mandatory option %s "
                                 "is wrong (should be %s).\n"
                                 % testKey, mandatoryVals[testKey])
                sane = False
        for testKey in mandatoryMinVals.keys():
            if not self.keyVals.has_key(testKey):
                sys.stderr.write("Error: mandatory option %s missing from "
                                 "template file.\n" % testKey)
                sane = False
            elif int(self.keyVals[testKey]) < mandatoryMinVals[testKey]:
                sys.stderr.write("Error: value for mandatory option %s "
                                 "is too low (should be at least %s).\n"
                                 % testKey, mandatoryMinVals[testKey])
                sane = False
        for testKey in mandatoryAbsent:
            if self.keyVals.has_key(testKey):
                sys.stderr.write("Error: option %s must not be present in "
                                 "template file.\n" % testKey)                
                sane = False
        return sane

    def createInput(self, map=True):
        """Create a runnable MNDO input file for the template molecule."""
        try:
            f=open(self.fileName, 'w')
        except IOError:
            sys.stderr.write("Error: could not open %s" +
                             "for writing.\n" % self.fileName)
            return -1
        # The following keywords should be added
        # jop=-2, ioutci=1, iuvcd=2
        # mprint=0 necessary for calc check
        mandatoryVals = { 'jop':'-2',
                          'ioutci':'1',
                          'iuvcd':'2',
                          'mprint':'0'}
        # Set the correct value of imomap
        if map:
            mandatoryVals['imomap'] = 1
        else:
            mandatoryVals['imomap'] = -1            
        tempKeyOrder = list(self.keyOrder)
        tempKeyVals = dict(self.keyVals)
        # Remove the plus sign later if there are no additions
        tempKeyOrder.append('+')
        for testKey in mandatoryVals:
            if not tempKeyVals.has_key(testKey):
                tempKeyOrder.append(testKey)  
            tempKeyVals[testKey] = mandatoryVals[testKey]
        if tempKeyOrder[-1] == '+':
            tempKeyOrder.pop()
        # Write keywords
        for key in tempKeyOrder:
            if key == '+':
                f.write('+\n')
            else:
                f.write('%s=%s ' % (key, tempKeyVals[key]))
        f.write('\n')
        # Write template file after keywords
        sLine = self.inputFile.getLastOptionLine() + 1
        for i in self.inputFile.fileContents[sLine:]:
            f.write(i)

    def extractQuantities(self, fFormalism):
        """Return calculated quantities from the MNDO output file.

        Quantities returned:
        Excited state energies (relative to ground state, in eV)
        Oscillator strengths (according to the formalism in fFormalism)

        """
        # Must be some output to analyse!
        if self.output is None:
            sys.stderr.write("Warning: could not extract quantities - "
                             "no output!\n")
            return -1
        (self.excitedStates, self.excitedStateEnergies,
         self.oscillatorStrengths) = \
             self.output.getTransitionProperties(fFormalism, True)
        self.noOfStates = len(self.excitedStates)


class SampleMol(Mol):
    
    """Structure of sample geometry"""
    
    def __init__(self, xyzFileObject, geomNumber):
        Mol.__init__(self)
        self.xyzGeom = xyzFileObject.getGeom(geomNumber)
        self.fileName = 'tmp.inp'
        self.selected = {}
        self.failCause = {}

    def createInput(self, templateMol, map=True):
        """Create an MNDO input file for the sample geometry."""
        # Write first part of job (template molecule) if mapping
        if map:
            templateMol.createInput(map)
            writeMode = 'a'
        else:
            writeMode = 'w'
        # Write second (keepci) part of job (sample geometry) if mapping
        # if not mapping this is the only part of the job
        try:
            f=open(self.fileName, writeMode)
        except IOError:
            sys.stderr.write("Error: could not open %s" +
                             "for writing.\n" % self.fileName)
            return -1
        # The following keywords should be added
        # jop=-2, ioutci=1, iuvcd=2
        mandatoryVals = { 'jop':'-2',
                          'ioutci':'1',
                          'iuvcd':'2',
                          'mprint':'0'}
        # Set the correct value of imomap
        if map:
            mandatoryVals['imomap'] = 1
            mandatoryVals['keepci'] = 1
        else:
            mandatoryVals['imomap'] = -1
            mandatoryVals['keepci'] = 0
        sampleKeyOrder = list(templateMol.keyOrder)
        sampleKeyVals = dict(templateMol.keyVals)
        # Remove the plus sign later if there are no additions
        sampleKeyOrder.append('+')
        for testKey in mandatoryVals:
            if not sampleKeyVals.has_key(testKey):
                sampleKeyOrder.append(testKey)  
            sampleKeyVals[testKey] = mandatoryVals[testKey]
        if sampleKeyOrder[-1] == '+':
            sampleKeyOrder.pop()
        # Write keywords
        for key in sampleKeyOrder:
            if key == '+':
                f.write('+\n')
            else:
                f.write('%s=%s ' % (key, sampleKeyVals[key]))
        f.write('\n')
        # Write title
        f.write('Sample geometry job\nAutomatically generated\n')        
        # Write sample geometry
        mndoGeom = self.xyzGeom.writeMNDOGeom(False)
        for i in range(len(mndoGeom)):
            f.write(mndoGeom[i] + "\n")
        # Write template file lines after geometry if present
        sLine = templateMol.inputFile.getLastGeomLine() + 1
        for i in templateMol.inputFile.fileContents[sLine:]:
            f.write(i)
        f.close()

    def createDynamicsRunInput(self, fileName, templateMol, title, gradients,
                               mapthr, map=True):
        """Create an MNDO input file for a dynamics run."""
        try:
            f=open(fileName, 'w')
        except IOError:
            sys.stderr.write("Error: could not open %s" +
                             "for writing.\n" % fileName)
            return -1
        # The following keywords should be added
        # ipubo=1, ktrial=11, ncisym=-1(?), imomap=2, icross=6, ncigrd
        # jop=-2, mapthr
        ncigrd = str(len(gradients))
        mapthrPercent = str(int(mapthr * 100))
        mandatoryVals = { 'jop':'-2',
                          'ipubo':'1',
                          'ktrial':'11',
                          'ncisym':'-1',
                          'icross':'6',
                          'ncigrd':ncigrd }
        # Set the correct value of imomap and mapthr
        if map:
            mandatoryVals['imomap'] = 1
            mandatoryVals['mapthr'] = mapthrPercent
        else:
            mandatoryVals['imomap'] = -1   
        dynKeyOrder = list(templateMol.keyOrder)
        dynKeyVals = dict(templateMol.keyVals)
        # Remove the plus sign later if there are no additions
        dynKeyOrder.append('+')
        for testKey in mandatoryVals:
            if not dynKeyVals.has_key(testKey):
                dynKeyOrder.append(testKey)  
            dynKeyVals[testKey] = mandatoryVals[testKey]
        if dynKeyOrder[-1] == '+':
            dynKeyOrder.pop()
        # Write keywords
        for key in dynKeyOrder:
            if key == '+':
                f.write('+\n')
            else:
                f.write('%s=%s ' % (key, dynKeyVals[key]))
        f.write('\n')
        # Write title
        f.write(title)
        # Write geometry
        mndoGeom = self.xyzGeom.writeMNDOGeom(False)
        for i in range(len(mndoGeom)):
            f.write(mndoGeom[i] + "\n")        
        # Write the mapped active orbitals and gradients
        if map:
            for a in self.activeOrbitals:
                f.write(' ' + str(a))
            f.write('\n')
        for g in gradients:
            f.write(' ' + str(g))
        f.write('\n')
        f.close()

    def extractQuantities(self, fFormalism, map=True):
        """Extract calculated quantities from the MNDO output file.

        Quantities extracted:
         self.noOfStates - total calculated excited states
         self.excitedStates - list of calculated excited states
         self.excitedStateEnergies - dictionary of excited state energies
                                     (relative to ground state, in eV)
         self.oscillatorStrengths - dictionary of oscillator strengths
                                    (according to the formalism in fFormalism)
              NB: The dictionary keys correspond to self.excitedStates
         self.activeOrbitals - list of active space MOs
         self.mapOverlap - list of mapping 'overlaps' (0.0 if mapping failed)
         self.ose2Ratio - dictionary of calculated f(ko)/dE(ko)^2 ratios

        Return -1 if there was no output to analyse

        """
        # Must be some output to analyse!
        if self.output is None:
            sys.stderr.write("Warning: could not extract quantities - "
                             "no output!\n")
            return -1
        (self.excitedStates, self.excitedStateEnergies,
         self.oscillatorStrengths) = \
             self.output.getTransitionProperties(fFormalism, True)
        if map:
            (self.activeOrbitals, self.mapOverlap) = \
                self.output.getMappingInfo(True)
        self.noOfStates = len(self.excitedStates)
        # Ratio f_ko/dE^2
        self.ose2Ratio = {}
        for i in range(self.noOfStates):
            s = self.excitedStates[i]
            ratio = self.oscillatorStrengths[s] / \
               (self.excitedStateEnergies[s] * self.excitedStateEnergies[s])
            self.ose2Ratio[s] = ratio

    def initFilter(self, allStates, targetState):
        """Initialise filter and filter out unwanted states.

        self.selected - a dictionary of states, for each state,
                        True: state is selected (has not been filtered out)
                        False: state is not selected (has been filtered out)

        self.failCause - if the state is filtered out in any of the filters,
                         the reason will be stored here as a string

        """
        for i in range(self.noOfStates):
            if allStates or self.excitedStates[i] == targetState:
                self.selected[self.excitedStates[i]] = True
            else:
                self.selected[self.excitedStates[i]] = False
                self.failCause[self.excitedStates[i]] = 'state'

    def filterMapping(self, mapThreshold):
        """Filter by active space mapping."""
        if len(self.selected) == 0:
            sys.stderr.write("Error: you must initialise the filter\n")
            return -1
        for i in range(len(self.mapOverlap)):
            if abs(self.mapOverlap[i]) < mapThreshold: 
                for j in range(self.noOfStates):
                    self.selected[self.excitedStates[j]] = False
                    self.failCause[self.excitedStates[j]] = 'map'
                return -1

    def filterWindow(self, target, tolerance): 
        """Filter by energy window."""
        if target is None:
            sys.stderr.write("Error: no window target energy specified.\n")
            return -1
        for i in range(self.noOfStates):
            if not self.selected[self.excitedStates[i]]:
                continue
            windowBottom = target[i] - tolerance
            windowTop = target[i] + tolerance
            eExc = self.excitedStateEnergies[self.excitedStates[i]]
            if eExc < windowBottom or eExc > windowTop:
                self.selected[self.excitedStates[i]] = False
                self.failCause[self.excitedStates[i]] = 'window'

    def filterProb(self, maxRatio, probFilter):
        """Filter by transition probability."""
        # Calculate transition probabilities for all states
        # and for selected states, determine stochastically
        # if they should be filtered out
        self.transProb = {}
        for i in self.excitedStates:
	    if maxRatio > 0.0:
	      self.transProb[i] = self.ose2Ratio[i] / maxRatio
	    else:
	      self.transProb[i] = 0.0
            if not self.selected[i]:
	      continue            
	    if probFilter:
	      rnd = random.random()
	      if rnd > self.transProb[i]:
		  self.selected[i] = False
		  self.failCause[i] = 'prob (rnd={:f})'.format(rnd)
                

class FilterLogFile(mndotools.RawFile):
    def __init__(self, fileName):
        """Load in file and write a contents list"""
        mndotools.RawFile.__init__(self, fileName)
        self.fFormalism = None
        self.noOfGeoms = 0
        self.geomContents = []
        readStage = 0
        for i in range(len(self.fileContents)):
            if self.fileContents[i] == 'SUMMARY OF FILTERING\n':
                readStage = 1
            elif self.fileContents[i] == 'INDIVIDUAL GEOMETRIES\n':
                readStage = 2
            else:
                line = self.fileContents[i].split()
                if len(line) < 1:
                    continue
                if readStage == 0 and line[0] == 'fFormalism:':
                    self.fFormalism = int(line[1])
                elif readStage == 1 and line[0] == 'Geometries:':
                    self.noOfGeoms = int(line[1])
                elif readStage == 2 and line[0] == 'Geom:':
                    self.geomContents.append(i)
        if self.noOfGeoms != len(self.geomContents):
            sys.stderr.write("Error: log file is not self consistent.\n")
            sys.exit(1)

    def extractLogQuantities(self, gNum, mol):
        """Extract transition properties/mapping information from log file."""
        # Mapping
        mol.activeOrbitals = []
        mol.mapOverlap = []
        s = self.geomContents[gNum]
        c = 1
        line = self.fileContents[s + c].split()
        if line[0] != 'AcOrb' or line[1] != 'Mapping':
            sys.stderr.write("Error: could not recognise mapping header.\n")
            sys.exit(1)
        while 1:
            c += 1
            line = self.fileContents[s + c].split()
            if line[0] == 'State':
                break
            mol.activeOrbitals.append(int(line[0]))
            mol.mapOverlap.append(float(line[1]))
        # State properties
        mol.excitedStates = []
        mol.excitedStateEnergies = {}
        mol.oscillatorStrengths = {}
        if line[0] != 'State' or line[1] != 'Energy' or line[2] != 'OscStr':
            sys.stderr.write("Error: could not recognise state properties "
                             "header.\n")
            sys.exit(1)
        while 1:
            c += 1
            line = self.fileContents[s + c].split()
            if len(line) < 3:
                break
            state = int(line[0])
            mol.excitedStates.append(state)
            mol.excitedStateEnergies[state] = float(line[1])
            mol.oscillatorStrengths[state] = float(line[2])
        mol.noOfStates = len(mol.excitedStates)
        # Ratio f_ko/dE^2 - recalculate to avoid loss of precision
        mol.ose2Ratio = {}
        for i in range(mol.noOfStates):
            s = mol.excitedStates[i]
            ratio = mol.oscillatorStrengths[s] / \
               (mol.excitedStateEnergies[s] * mol.excitedStateEnergies[s])
            mol.ose2Ratio[s] = ratio
                

if __name__ == '__main__':
    main()
