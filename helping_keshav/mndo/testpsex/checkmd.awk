#
# checkmd.awk
#
# Tom Keal 21.5.2007
#
# Checking Born-Oppenheimer MD at SCF level added in July 2010 by
# Axel Koslowski
#
# ---------------------------------------------------------------------
# Checks the highly structured outputs of the MNDO dynamics driver 
# against a pre-calculated reference
#
# Usage: awk -f checkmd.awk name=<testjob name> <reference file>
#
# Format of the reference file:  
# CHECK: <type> 
# <reference output data>
# ENDCHECK
#
# Notes:
# - <testjob name> should have no file ending
# - <type> can currently only be hopping
# - <reference output data> should be formatted as a normal output
#   file but need not contain all the output data, e.g. for checking
#   hopping you only need to include the MD steps you want to check,
#   but each MD step output must be complete. The MD steps in the
#   reference should also be in sequential order.
# - Multiple checks on a single testjob can be carried out by adding
#   extra CHECK...ENDCHECK blocks.
# - Comments can be written between CHECK blocks but not within them. 
#
# Output:
# The output is similar to the output of check.awk, i.e. if all checks
# are passed, a 0-sized .pass file will be written, otherwise a .fail
# file will be written with details of the failure
#
# ---------------------------------------------------------------------
#
BEGIN{
  # bad contains any errors that occur during the check
  bad = "";
}

# The .mdchk file consists of reference files separated by the line
# CHECK: <filetype> <calculation file name>
#
/CHECK:/ {
  if (NF != 2) {
    bad = bad "Wrong number of arguments in CHECK statement\n";
  } else if ($2 == "hopping") {
    calcFile = name ".hopping.out";
    CheckHopping(calcFile);
  } else if ($2 == "scfbomd") {
    calcFile = name ".out";
    CheckSCFBOMD(calcFile);
  } else {
    bad = bad "Unrecognised file type: " $2 "\n";
  }
}

END{
  # Following the standard procedure of check.awk,
  # Print out a 0 size .pass file if the check was passed
  # Print out a .fail file with error message if the check failed
  if (bad != "") {
    printf "%s", bad;
    print bad > (name ".fail");
  }
  else {
    print "Passed" ;
    printf "" > (name ".pass");
  }
}


function CheckHopping(calcFile)
{
  # Tolerance settings for comparison
  tolheat = 0.00001;
  tolgrad = 0.00001;
  tolpop = 0.00001;
  tolprob = 0.00001;
  tolrnd = 0.00001;

  refMDStep = -1; 
  calcMDStep = -1;
  refState = 0;
  calcState = 0;

  # loop over each mdstep.
  # note we return immediately from the function only for unrecoverable errors
  # (bad formatting of files etc.)
  while (1) {

    getline;
    if ($1 == "ENDCHECK") {
      return;
    } else if ($1 != "mdstep:") {
      bad = bad "CheckHopping: mdstep label expected in ref file but not found\n";
      return;
    }

    refMDStep = $2;
    refState = $4;

    getline < calcFile;

    if ($1 != "mdstep:") {
      bad = bad "CheckHopping: mdstep label expected in calc file but not found\n";
      return;
    }

    calcMDStep = $2;
    calcState = $4;

    # It is only necessary to compare against those MD steps which are present
    # in the check file. Skip all other MD steps in the calculation output.
    if (calcMDStep < refMDStep) {
      calcMDStep = SkipMDSteps(refMDStep);
      if (calcMDStep == -1) {
	bad = bad sprintf("CheckHopping: Could not find calc MD step %d\n", refMDStep);
	return;
      }
    }

    if (calcMDStep != refMDStep) {
      bad = bad sprintf("CheckHopping: calc MD step %d and ref MD step %d do not match\n",  calcMDStep, refMDStep);
      return;
    }

    # Check electronic state. 
    if (calcState != refState) {
      bad = bad sprintf("CheckHopping: Wrong electronic state at MD step %d: %d instead of %d\n", 
			refMDStep, calcState, refState);
    }

    # Check energies
    # NB: (1) Assumes energies are on one line
    # NB: (2) Number of states taken from the ref energy line
    getline;
    if ($1 != "energy") {
      bad = bad "CheckHopping: energy label expected in ref file but not found\n";
      return;
    }
    numStates = NF - 1;
    # We put the calculated energies into a variable so that reference info isn't lost
    getline calcLine < calcFile;
    split(calcLine, calcVals);
    if (calcVals[1] != "energy") {
      bad = bad "CheckHopping: energy label expected in calc file but not found\n";
      return;
    }
    for (i = 2; i <= NF; i++) {
      if (ToleranceCheck(calcVals[i], $i, tolheat)) {
	bad = bad sprintf("CheckHopping: Wrong energy for state %d in MD step %d: %.5f instead of %.5f\n",
			  i - 1, refMDStep, calcVals[i], $i);
      }
    }

    # Check energy differences
    # NB: Assumes energy differences are arranged in a triangle 
    #     with (numState-1) rows of 1, 2, 3, ... numbers
    getline;
    if ($0 != "energy differences") {
      bad = bad sprintf("CheckHopping: expected energy differences label in ref file but not found\n");
      return;
    }
    getline < calcFile;
    if ($0 != "energy differences") {
      bad = bad sprintf("CheckHopping: expected energy differences label in calc file but not found\n");
      return;
    }
    for (i = 1; i < numStates; i++) {
      getline;
      getline calcLine < calcFile;
      split(calcLine, calcVals);
      for (j = 1; j <= i; j++) {
	if (ToleranceCheck(calcVals[j], $j, tolheat)) {
	  bad = bad sprintf("CheckHopping: Wrong energy difference between states %d and %d in MD step %d: %.5f instead of %.5f\n",
			  j, i+1, refMDStep, calcVals[j], $j);
	}
      }
    }

    # Check couplings
    # NB: Assumes couplings are formatted in the same way as energy differences 
    getline;
    if ($0 != "coupling") {
      bad = bad sprintf("CheckHopping: expected coupling label in ref file but not found\n");
      return;
    }
    getline < calcFile;
    if ($0 != "coupling") {
      bad = bad sprintf("CheckHopping: expected coupling label in calc file but not found\n");
      return;
    }
    for (i = 1; i < numStates; i++) {
      getline;
      getline calcLine < calcFile;
      split(calcLine, calcVals);
      for (j = 1; j <= i; j++) {
	if (ToleranceCheck(calcVals[j], $j, tolgrad)) {
	  bad = bad sprintf("CheckHopping: Wrong coupling between states %d and %d in MD step %d: %.5f instead of %.5f\n",
			  j, i+1, refMDStep, calcVals[j], $j);
	}
      }
    }

    # Check population and probability
    # NB: Assumes one population and probability per line, with numStates in total
    getline;
    if ($0 != "population   probability") {
      bad = bad sprintf("CheckHopping: expected population/probability label in ref file but not found\n");
      return;
    }
    getline < calcFile;
    if ($0 != "population   probability") {
      bad = bad sprintf("CheckHopping: expected population/probability label in calc file but not found\n");
      return;
    }
    for (i = 1; i <= numStates; i++) {
      getline;
      getline calcLine < calcFile;
      split(calcLine, calcVals);
      if (ToleranceCheck(calcVals[1], $1, tolpop)) {
	bad = bad sprintf("CheckHopping: Wrong population for state %d in MD step %d: %.5f instead of %.5f\n",
			  i, refMDStep, calcVals[1], $1);
      }
      if (ToleranceCheck(calcVals[2], $2, tolprob)) {
	bad = bad sprintf("CheckHopping: Wrong probability for state %d in MD step %d: %.5f instead of %.5f\n",
			  i, refMDStep, calcVals[2], $2);
      }
    }

    # Check "random" number - in reality it will have been seeded to ensure reproducibility
    getline;
    if ($1 != "rnd" || $2 != "number") {
      bad = bad sprintf("CheckHopping: expected rnd number label in ref file but not found\n");
      return;
    }      
    refRndNumber = $3;
    getline < calcFile;
    if ($1 != "rnd" || $2 != "number") {
      bad = bad sprintf("CheckHopping: expected rnd number label in calc file but not found\n");
      return;
    }  
    calcRndNumber = $3;
    if (ToleranceCheck(calcRndNumber, refRndNumber, tolrnd)) {
      bad = bad sprintf("CheckHopping: Wrong rnd number (sic!) for MD step %d: %.5f instead of %.5f\n",
			refMDStep, calcRndNumber, refRndNumber);
    }    

    # Finally skip the dashed line
    getline;
    getline < calcFile;

  } # end of mdstep loop

}

function SkipMDSteps(refMDStep) {
  while (1) {
    eof = getline < calcFile;
    if (eof <= 0) {
      return -1;
    }
    if ($1 == "mdstep:") {
      calcMDStep = $2;
      if (calcMDStep == refMDStep) {
	return calcMDStep;
      }
      if (calcMDStep > refMDStep) {
	return -1;
      }
    }
  }
}


function CheckSCFBOMD(calcFile)
{
  # Tolerance settings for comparison
  tolHeat  = 0.00001;  # Tolerance for the SCF heat of formation.
  tolNorm  = 0.00001;  # Tolerance for the Cartesian gradient norm.

  refStep  = -1;
  refHeat  = 0.;
  refNorm  = 0.;

  calcStep = 0;

  # Loop over lines in the standard input.
  # One line corresponds to one MD step.

  while ((readstat = getline) > 0)
  {
    if ($1 == "ENDCHECK")
    {
      return;
    }
    else
    {
      refStep = $1;
      refHeat = $2;
      refNorm = $3;
    }

    if (refStep > 0)
    {
      # Find specified MD step.
      calcStep = FindMDStep(refStep, calcFile);

      if (calcStep < 0)
      {
	bad = bad sprintf("CheckSCFBOMD: Could not find MD step %d in output file.\n", refStep);
	return;
      }
    }

    if (calcStep != refStep)
    {
      bad = bad sprintf("CheckSCFBOMD: MD step numbers %d and %d do not match.\n", refStep, calcStep);
      return;
    }

    # Loop over lines in the output of the calculation.
    # Check SCF heat of formation.

    while ((calcstat = getline < calcFile) > 0)
    {
      if ($1 == "SCF"  &&  $2 == "HEAT"  &&  $3 == "OF"  &&  $4 == "FORMATION")
      {
        if (ToleranceCheck($5, refHeat, tolHeat))
          bad = bad sprintf("CheckSCFBOMD: Wrong SCF heat of formation in MD step %d: %.5f instead of %.5f\n",
                            refStep, $5, refHeat);
        break;
      }
    }

    if (calcstat < 0)
    {
      bad = bad sprintf("Error reading file '%s': %s\n", calcFile, ERRNO);
      return;
    }

    if (calcstat == 0)
    {
      bad = bad "Unexpected end of output file while reading SCF heat of formation.\n";
      return;
    }

    # Loop over lines in the output of the calculation.
    # Check Cartesian gradient norm.

    while ((calcstat = getline < calcFile) > 0)
    {
      if ($1 == "CARTESIAN"  &&  $2 == "GRADIENT"  &&  $3 == "NORM")
      {
        if (ToleranceCheck($4, refNorm, tolNorm))
          bad = bad sprintf("CheckSCFBOMD: Wrong Cartesian gradient norm in MD step %d: %.5f instead of %.5f\n",
                            refStep, $4, refNorm);
        break;
      }
    }

    if (calcstat < 0)
    {
      bad = bad sprintf("Error reading file '%s': %s\n", calcFile, ERRNO);
      return;
    }

    if (calcstat == 0)
    {
      bad = bad "Unexpected end of output file while reading Cartesian gradient norm.\n";
      return;
    }
  }

  if (readstat < 0)
    bad = bad sprintf("Error reading standard input: %s\n", ERRNO);
}


function FindMDStep(refStep, calcFile)
{
  while ((status = getline < calcFile) > 0)
  {
    if ($1 == "MD"  &&  $2 == "step:")
    {
      calcStep = $3;

      if (calcStep == refStep)
	return calcStep;

      if (calcStep > refStep)
	return -1;
    }
  }

  if (status < 0)
    bad = bad sprintf("Error reading file '%s': %s\n", calcFile, ERRNO);

  return -1;
}


function ToleranceCheck(calcNumber, refNumber, tolerance)
{
  # Return values:
  # 0 - OK
  # 1 - calculated value outside tolerance
  # 2 - numbers missing
  if (calcNumber == "" || refNumber == "") {
    return 2;
  }
  if (calcNumber < refNumber - tolerance || calcNumber > refNumber + tolerance) {
    return 1;
  }
  return 0;
}
