#!/usr/bin/python
#
# --------------------------
#
# mdsample.py
#
# Tom Keal, August 2007
#
# --------------------------
#
"""mdsample - sample a molecular dynamics trajectory

 PURPOSE

 mdsample takes a trajectory file (in Molden or XYZ format) and outputs
 a sample of trajectory snapshots.

 The sample can be either random or taken at fixed intervals. By default,
 the output is a randomised list of all the geometries, but normally a
 smaller sample would be specified (using --sample-size).

 The sample can then be used as input to mdfilter.py.

 If the --vel option is selected, the velocities corresponding to the
 sampled structures are taken from an MNDO velocity output file and 
 output to a velocity snapshot file (in the same order as the trajectory 
 snapshot file). This can also be used as an input to mdfilter.py if you 
 want to provide initial velocities as well as initial structures.

 REQUIREMENTS

 Requires the mndotools.py library

 USAGE

 mdsample.py [OPTIONS]

 The default input file is 'traj.out'
 The default output file is 'sample.xyz'
 The default log file is 'mdsample.log'

 The default velocity input file is 'vel.out'
 The default velocity output file is 'sample.vel'

 OPTIONS

 -h, --help
        display this help message

 -q, --quiet
         don't print anything to standard output

 --trajectory=(file name)
        set a new input file name

 --sample=(file name)
        set a new output file name

 --log=(file name)
        set a new log file name

 --vel
        output a sample velocity file corresponding to the sample trajectory
        file.

 --vel-inp=(file name)
        set a new velocity input file name

 --vel-out=(file name)
        set a new velocity output file name

 --sample-size=(int)
        set the number of sample geometries required
        (default: as many as possible)

 --random-seed=(int)
        set the seed for the random number generator

 --step=(int)
        set a fixed interval for sampling instead of taking a random sample
        A step of 1 means that all geometries will be sampled.

 --step-start=(int)
        set the first geometry to be sampled in the case of fixed step
        sampling, numbered from 0. (default: step size - 1)

"""
import sys
import os
import getopt
import random
import mndotools

def main():
    """Routine for sampling geometries."""
    # Set default global options and parse command line options
    options = GlobalOptions()
    options.parseCommandLine()

    # Load in files
    if options.writeStdOut:
        sys.stdout.write('Reading trajectory file...\n')
    trajFile = mndotools.XYZFile(options.fileInpTrajectory)
    if options.sampleVelocities:
        if options.writeStdOut:
            sys.stdout.write('Reading velocity file...\n')
        velFile = mndotools.VelocityFile(options.fileInpVelocities)

    # Determine the sample
    if options.writeStdOut:
        sys.stdout.write('Sampling geometries...\n')
    sample = []
    if options.random:
        if options.userRandomSeed is None:
            random.seed()
        else:
            random.seed(options.userRandomSeed)
        if options.userSampleSize is None:
            sample = range(0, trajFile.noOfGeoms)
            random.shuffle(sample)
        else:
            if options.userSampleSize > trajFile.noOfGeoms:
                sys.stderr.write('Error: number of sample geometries '
                                 'requested is too large.\n')
                sys.exit(1)
            sample = random.sample(range(0, trajFile.noOfGeoms),
                                   options.userSampleSize)
    else:
        if options.step is None or options.stepStart is None:
            sys.stderr.write('Error: step size or start not known.\n')
            sys.exit(1)
        sample = range(options.stepStart, trajFile.noOfGeoms,
                       options.step)
        if options.userSampleSize is not None:
            if options.userSampleSize > len(sample):
                sys.stderr.write('Error: number of sample geometries '
                                 'requested is too large.\n')
                sys.exit(1)
            else:
                sample = sample[:options.userSampleSize]
                     
    # Write out the sample
    if options.writeStdOut:
        sys.stdout.write('Writing sample geometry file...\n')    
    trajFile.writeSnapshots(sample, options.fileOutSample)
    if options.sampleVelocities:
        if options.writeStdOut:
            sys.stdout.write('Writing sample velocity file...\n')    
        velFile.writeSnapshots(sample, options.fileOutVelocities)
    

    # Write out the log file
    if options.writeStdOut:
        sys.stdout.write('Writing log file...\n')
    try:
        logFile = open(options.fileOutLog, 'w')
    except IOError:            
        sys.stderr.write("Error: could not open %s" +
                         "for writing.\n" % options.fileOutLog)    
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
    # Summary of sampling
    logFile.write("SUMMARY OF SAMPLING\n")
    logFile.write("Sample size: %i\n" % len(sample))
    logFile.close()
    

class GlobalOptions:
    
    """Set and store global options."""
    
    def __init__(self):
        """Set default options."""
        # Program features
        self.writeStdOut = True
        self.sampleVelocities = False
        # File I/O
        self.fileInpTrajectory = 'traj.out'
        self.fileOutSample = 'sample.xyz'
        self.fileOutLog = 'mdsample.log'
        self.fileInpVelocities = 'vel.out'
        self.fileOutVelocities = 'sample.vel'
        # Sampling options
        self.userSampleSize = None
        self.random = True
        self.userRandomSeed = None
        self.step = None
        self.stepStart = None

    def parseCommandLine(self):
        """Parse command line options."""
        shortOptions = "hq"
        longOptions = ["help", "quiet",
                       "trajectory=", "sample=", "log=",
                       "vel", "vel-inp=", "vel-out=",
                       "sample-size=", "random-seed=",
                       "step=", "step-start=" ]
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
            if o in ("-q", "--quiet"):
                self.writeStdOut = False
            if o in ("--vel"):
                self.sampleVelocities = True
            # Change default files
            if o == "--trajectory":
                self.fileInpTrajectory = a
            if o == "--sample":
                self.fileOutSample = a
            if o == "--log":
                self.fileOutLog = a
            if o == "--vel-inp":
                self.fileInpVelocities = a
            if o == "--vel-out":
                self.fileOutVelocities = a
            # Change default parameters
            if o == "--sample-size":
                try:
                    self.userSampleSize = int(a)
                except ValueError:
                    sys.stderr.write("Error: sample-size should be an "
                                     "integer.\n")
                    sys.exit(2)
                if self.userSampleSize < 1:
                    sys.stderr.write("Error: sample-size must be greater "
                                     "than zero.\n")
                    sys.exit(2)                    
            if o == "--random-seed":
                try:
                    self.userRandomSeed = int(a)            
                except ValueError:
                    sys.stderr.write("Error: random-seed must be an integer."
                                     "\n")
                    sys.exit(2)
            if o == "--step":
                self.random = False
                try:
                    self.step = int(a)            
                except ValueError:
                    sys.stderr.write("Error: step must be an integer.\n")
                    sys.exit(2)
                if self.step < 1:
                    sys.stderr.write("Error: step must be greater "
                                     "than zero.\n")
                    sys.exit(2)                   
            if o == "--step-start":
                try:
                    self.stepStart = int(a)            
                except ValueError:
                    sys.stderr.write("Error: step-start must be an integer.\n")
                    sys.exit(2)
                if self.stepStart < 0:
                    sys.stderr.write("Error: step-start must be positive.\n")
                    sys.exit(2)  
        # There should be no remaining arguments
        for arg in args:
            sys.stderr.write("Warning: argument '%s' not recognised.\n" % arg)
        # Set stepStart if not already
        if self.step is not None and self.stepStart is None:
            self.stepStart = self.step - 1
        # Check that the input file exists
        if not os.path.exists(self.fileInpTrajectory):
            sys.stderr.write("Error: trajectory file %s does not exist.\n"
                             % self.fileInpTrajectory)
            sys.exit(2)


if __name__ == '__main__':
    main()
