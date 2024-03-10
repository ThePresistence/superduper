#!/usr/bin/php
<?php
// Check number of arguments.
if ($argc < 2  ||  $argc > 3  ||  ($argc == 3  &&  $argv[1] !== "--gugaci"))
{
  fprintf(STDERR, "USAGE: %s [--gugaci] <target file name>\n", $argv[0]);
  return 1;
}

// Check use of the GUGA-CI module.
$gugaci = $argc == 3;

// Analyze target file name.
$target = basename($argv[$argc - 1], ".inp");

if (sscanf($target, "job%d%s", $jobNumber, $jobLetter) != 2)
{
  fprintf(STDERR, "Invalid target file name »%s«.\n", $argv[$argc - 1]);
  return 1;
}

// Write target input to STDOUT.
PrintInput($jobNumber, $jobLetter, "MNDO", $gugaci);
PrintInput($jobNumber, $jobLetter, "AM1",  $gugaci);
PrintInput($jobNumber, $jobLetter, "PM3",  $gugaci);
PrintInput($jobNumber, $jobLetter, "OM1",  $gugaci);
PrintInput($jobNumber, $jobLetter, "OM2",  $gugaci);
PrintInput($jobNumber, $jobLetter, "OM3",  $gugaci);
PrintInput($jobNumber, $jobLetter, "ODM2", $gugaci);
PrintInput($jobNumber, $jobLetter, "ODM3", $gugaci);


class Options
{
  public  $jobNumber;
  public  $jobLetter;
  private $methodKeyword;
  private $pointChargesFlag;
  private $symmetryFlag;
  private $gugaciFlag;
  public  $listOfStates;
  public  $jop;
  public  $iform;
  public  $kharge;
  public  $imult;
  public  $mminp;
  public  $numatm;
  public  $mmcoup;
  public  $kci;
  public  $ici1;
  public  $ici2;
  public  $nciref;
  public  $mciref;
  public  $levexc;
  public  $iroot;
  public  $multci;
  public  $iuvcd;
  public  $maxdav;
  public  $kitdav;


  function __construct($jobNumber, $jobLetter, $method)
  {
    $this->jobNumber        = $jobNumber;
    $this->jobLetter        = $jobLetter;
    $this->methodKeyword    = $method;
    $this->pointChargesFlag = false;
    $this->symmetryFlag     = true;
    $this->gugaciFlag       = false;
    $this->jop              = -2;
    $this->iform            =  1;
    $this->kharge           =  0;
    $this->imult            =  0;
    $this->mminp            =  0;
    $this->numatm           =  0;
    $this->mmcoup           =  0;
    $this->kci              =  0;
    $this->ici1             =  0;
    $this->ici2             =  0;
    $this->ioutci           =  2;
    $this->nciref           =  0;
    $this->mciref           =  0;
    $this->levexc           =  0;
    $this->iroot            =  0;
    $this->multci           =  0;
    $this->iuvcd            =  0;
    $this->maxdav           =  0;
    $this->kitdav           =  0;
  }


  function SetPointCharges()
  {
    $this->pointChargesFlag = true;
    $this->mminp            = 2;
    $this->mmcoup           = 2;
  }


  function IsSetPointCharges()
  {
    return $this->pointChargesFlag;
  }


  function UnsetSymmetry()
  {
    $this->symmetryFlag = false;
  }


  function IsSetSymmetry()
  {
    return $this->symmetryFlag;
  }


  function SetGUGACI($gugaciFlag)
  {
    $this->gugaciFlag = $gugaciFlag;
  }


  function IsSetGUGACI()
  {
    return $this->gugaciFlag;
  }


  function SetCIS()
  {
    if ($this->IsSetGUGACI())
    {
      $this->kci    = 5;
      $this->nciref = 1;
      $this->levexc = 1;
    }
    else
    {
      $this->kci    = 6;
    }
  }


  function SetSASFCIS()
  {
    // Not possible with GUGA-CI.
    $this->SetGUGACI(false);
    $this->kci    = 7;
    $this->nciref = 1;
  }


  function SetSFXCIS()
  {
    if ($this->IsSetGUGACI())
    {
      $this->kci    = 5;
      $this->nciref = 3;
      $this->levexc = 1;
    }
    else
    {
      $this->kci    = 7;
      $this->nciref = 3;
    }
  }


  function SetRPA()
  {
    // Not possible with GUGA-CI.
    $this->SetGUGACI(false);
    $this->kci = 8;
  }


  function SetSingletCI()
  {
    if ($this->imult > 1)
      $this->multci = 1;
  }


  function SetTripletCI()
  {
    if ($this->imult !== 3)
      $this->multci = 3;
  }


  function SetMultipleStates($jobNumber, $jobLetter, $method, ...$istate)
  {
    if ($this->jobNumber     === $jobNumber  &&
        $this->jobLetter     === $jobLetter  &&
	$this->methodKeyword === $method)
    {
      $this->listOfStates = $istate;
      var_dump($istate);
    }
  }


  function IsTypeMNDO()
  {
    return
      $this->methodKeyword === "MNDO"  ||
      $this->methodKeyword === "AM1"   ||
      $this->methodKeyword === "PM3";
  }


  function PrintReferenceConfigurations()
  {
    $norbs = $this->ici1 + $this->ici2;

    switch ($this->imult)
    {
      case 0:
        $nelec = 2 * $this->ici1;
        $ihomo = $this->ici1 - 1;
        $ilumo = $this->ici1;
        break;
      case 1:
      case 3:
        $nelec = 2 * $this->ici1 - 2;
        $ihomo = $this->ici1 - 2;
        $ilumo = $this->ici1 - 1;
        break;
      case 2:
        $nelec = 2 * $this->ici1 - 1;
        $ihomo = $this->ici1 - 1;
        $ilumo = $this->ici1;
        break;
      default:
        printf("Unexpected value (%d) of option imult.\n", $this->imult);
        exit(1);
    }

    if ($this->nciref === 1  ||  $this->nciref === 3)
    {
      for ($i=0; $i<$norbs; ++$i)
      {
        $occ[$i]  = min($nelec, 2);
        $nelec   -= $occ[$i];
      }

      $this->PrintReferenceConfiguration($occ);
    }

    if ($this->nciref === 3)
    {
      $occ[$ihomo] = 1;
      $occ[$ilumo] = 1;
      $this->PrintReferenceConfiguration($occ);

      $occ[$ihomo] = 0;
      $occ[$ilumo] = 2;
      $this->PrintReferenceConfiguration($occ);
    }
  }


  function PrintReferenceConfiguration($occ)
  {
    $norbs = count($occ);

    for ($k=0; $k<$norbs; $k+=20)
    {
      $limit = min($k + 20, $norbs);

      for ($i=$k; $i<$limit; ++$i)
        printf("%4d", $occ[$i]);

      printf("\n");
    }
  }


  function PrintOptions()
  {
    $line = "";
    $this->AppendOption($line, $this->methodKeyword);

    if ($this->IsSetSymmetry())   $this->AppendOption($line, "SYM");
    if ($this->jop)               $this->AppendOption($line,    "jop={$this->jop}");
    if ($this->iform)             $this->AppendOption($line,  "iform={$this->iform}");
    if ($this->kharge)            $this->AppendOption($line, "kharge={$this->kharge}");
    if ($this->imult)             $this->AppendOption($line,  "imult={$this->imult}");
    if ($this->mminp)             $this->AppendOption($line,  "mminp={$this->mminp}");
    if ($this->numatm)            $this->AppendOption($line, "numatm={$this->numatm}");
    if ($this->mmcoup)            $this->AppendOption($line, "mmcoup={$this->mmcoup}");
    
    if (!$this->IsTypeMNDO())
    {
      $this->AppendOption($line, "icuts=-1");
      $this->AppendOption($line, "icutg=-1");
    }

    if ($this->kci)               $this->AppendOption($line,    "kci={$this->kci}");
    if ($this->ici1)              $this->AppendOption($line,   "ici1={$this->ici1}");
    if ($this->ici2)              $this->AppendOption($line,   "ici2={$this->ici2}");
    if ($this->ioutci)            $this->AppendOption($line, "ioutci={$this->ioutci}");
    if ($this->nciref)            $this->AppendOption($line, "nciref={$this->nciref}");
    if ($this->mciref)            $this->AppendOption($line, "mciref={$this->mciref}");
    if ($this->levexc)            $this->AppendOption($line, "levexc={$this->levexc}");
    if ($this->iroot)             $this->AppendOption($line,  "iroot={$this->iroot}");
    if ($this->multci)            $this->AppendOption($line, "multci={$this->multci}");
    if ($this->iuvcd)             $this->AppendOption($line,  "iuvcd={$this->iuvcd}");
    if ($this->maxdav)            $this->AppendOption($line, "maxdav={$this->maxdav}");
    if ($this->kitdav)            $this->AppendOption($line, "kitdav={$this->kitdav}");

    printf("%s\n", $line);
  }


  function AppendOption(&$line, $optstr)
  {
    $line2 = $line . " " . $optstr;

    if (strlen($line2) > 77)
    {
      printf("%s +\n", $line);
      $line = " " . $optstr;
    }
    else
      $line = $line2;
  }
};


function PrintInput($jobNumber, $jobLetter, $method, $gugaciFlag)
{
  $options = new Options($jobNumber, $jobLetter, $method);
  $options->SetGUGACI($gugaciFlag);
  $options->iroot  = 15;
  $options->iuvcd  =  2;

  switch ($jobNumber)
  {
    case 1:
      ClosedShellMoleculeSettings($options);
      MultipleStateSettings($options);
      Benzene($options);
      break;
    case 2:
      ClosedShellMoleculeSettings($options);
      MultipleStateSettings($options);
      $options->SetPointCharges();
      Benzene($options);
      break;
    case 3:
      ClosedShellMoleculeSettings($options);
      MultipleStateSettings($options);
      CyclopentadienylAnion($options);
      break;
    case 4:
      ClosedShellMoleculeSettings($options);
      MultipleStateSettings($options);
      $options->SetPointCharges();
      CyclopentadienylAnion($options);
      break;
    case 5:
      ClosedShellMoleculeSettings($options);
      MultipleStateSettings($options);
      CycloheptadienylCation($options);
      break;
    case 6:
      ClosedShellMoleculeSettings($options);
      MultipleStateSettings($options);
      $options->SetPointCharges();
      CycloheptadienylCation($options);
      break;
    case 7:
      BiradialSettings($options);
      MultipleStateSettings($options);
      Oxygen($options);
      break;
    case 8:
      BiradialSettings($options);
      MultipleStateSettings($options);
      $options->SetPointCharges();
      Oxygen($options);
      break;
    case 9:
      BiradialSettings($options);
      MultipleStateSettings($options);
      EthyleneD2d($options);
      break;
    case 10:
      BiradialSettings($options);
      MultipleStateSettings($options);
      $options->SetPointCharges();
      EthyleneD2d($options);
      break;
    case 11:
      BiradialSettings($options);
      MultipleStateSettings($options);
      YangsBiradicalTrimmed($options);
      break;
    case 12:
      BiradialSettings($options);
      MultipleStateSettings($options);
      $options->SetPointCharges();
      YangsBiradicalTrimmed($options);
      break;
    case 13:
      RadicalSettings($options);
      MultipleStateSettings($options);
      MethylRadical($options);
      break;
    case 14:
      RadicalSettings($options);
      MultipleStateSettings($options);
      $options->SetPointCharges();
      MethylRadical($options);
      break;
    case 15:
      RadicalSettings($options);
      MultipleStateSettings($options);
      NitricOxide($options);
      break;
    case 16:
      RadicalSettings($options);
      MultipleStateSettings($options);
      $options->SetPointCharges();
      NitricOxide($options);
      break;
    default:
      fprintf(STDERR, "Invalid job number (%d).\n", $jobNumber);
      exit(1);
  }
}


function ClosedShellMoleculeSettings($options)
{
  switch ($options->jobLetter)
  {
    case "a":
      $options->SetCIS();
      $options->SetSingletCI();
      break;
    case "b":
      $options->SetCIS();
      $options->SetTripletCI();
      if ($options->IsSetGUGACI())  $options->mciref = 1;
      break;
    case "c":
      $options->SetSASFCIS();
      $options->SetSingletCI();
      break;
    case "d":
      $options->SetSASFCIS();
      $options->SetTripletCI();
      break;
    case "e":
      $options->SetSFXCIS();
      $options->SetSingletCI();
      break;
    case "f":
      $options->SetSFXCIS();
      $options->SetTripletCI();
      if ($options->IsSetGUGACI())  $options->mciref = 1;
      break;
    case "g":
      $options->SetRPA();
      $options->SetSingletCI();
      break;
    case "h":
      $options->SetRPA();
      $options->SetTripletCI();
      break;
    default:
      fprintf(STDERR, "Invalid job letter (%s) for job number %d.\n",
                      $options->jobLetter, $options->jobNumber);
      exit(1);
  }
}


function BiradialSettings($options)
{
  switch ($options->jobLetter)
  {
    case "a":
      $options->SetCIS();
      $options->imult = 1;
      $options->SetSingletCI();
      break;
    case "b":
      $options->SetCIS();
      $options->imult = 1;
      $options->SetTripletCI();
      break;
    case "c":
      $options->SetCIS();
      $options->imult = 3;
      $options->SetSingletCI();
      break;
    case "d":
      $options->SetCIS();
      $options->imult = 3;
      $options->SetTripletCI();
      break;
    case "e":
      $options->SetSASFCIS();
      $options->imult = 1;
      $options->SetSingletCI();
      break;
    case "f":
      $options->SetSASFCIS();
      $options->imult = 1;
      $options->SetTripletCI();
      break;
    case "g":
      $options->SetSASFCIS();
      $options->imult = 3;
      $options->SetSingletCI();
      break;
    case "h":
      $options->SetSASFCIS();
      $options->imult = 3;
      $options->SetTripletCI();
      break;
    case "i":
      $options->SetSFXCIS();
      $options->imult = 1;
      $options->SetSingletCI();
      break;
    case "j":
      $options->SetSFXCIS();
      $options->imult = 1;
      $options->SetTripletCI();
      break;
    case "k":
      $options->SetSFXCIS();
      $options->imult = 3;
      $options->SetSingletCI();
      break;
    case "l":
      $options->SetSFXCIS();
      $options->imult = 3;
      $options->SetTripletCI();
      break;
    case "m":
      $options->SetRPA();
      $options->imult = 1;
      $options->SetSingletCI();
      break;
    case "n":
      $options->SetRPA();
      $options->imult = 1;
      $options->SetTripletCI();
      break;
    case "o":
      $options->SetRPA();
      $options->imult = 3;
      $options->SetSingletCI();
      break;
    case "p":
      $options->SetRPA();
      $options->imult = 3;
      $options->SetTripletCI();
      break;
    default:
      fprintf(STDERR, "Invalid job letter (%s) for job number %d.\n",
                      $options->jobLetter, $options->jobNumber);
      exit(1);
  }

  if ($options->IsSetGUGACI())
    $options->mciref = 1;
}


function RadicalSettings($options)
{
  switch ($options->jobLetter)
  {
    case "a":
      $options->SetCIS();
      break;
    case "b":
      $options->SetSASFCIS();
      break;
    case "c":
      $options->SetSFXCIS();
      break;
    case "d":
      $options->SetRPA();
      break;
    default:
      fprintf(STDERR, "Invalid job letter (%s) for job number %d.\n",
                      $options->jobLetter, $options->jobNumber);
      exit(1);
  }
}


function MultipleStateSettings($options)
{
      /*
      options->SetMultipleStates("MNDO", 1, 2, 3, 8, 11);
      options->SetMultipleStates("AM1",  1, 2, 3, 6,  9);
      options->SetMultipleStates("PM3",  1, 2, 3, 6,  9);
      options->SetMultipleStates("OM1",  1, 2, 3, 4,  9);
      options->SetMultipleStates("OM2",  1, 2, 3, 4,  7);
      options->SetMultipleStates("OM3",  1, 2, 3, 4,  7);
      options->SetMultipleStates("ODM2", 1, 2, 3, 4,  7);
      options->SetMultipleStates("ODM3", 1, 2, 3, 4,  7);
      */
      /*
      options->SetMultipleStates("MNDO", 1);
      options->SetMultipleStates("AM1",  1);
      options->SetMultipleStates("PM3",  1);
      options->SetMultipleStates("OM1",  1);
      options->SetMultipleStates("OM2",  1);
      options->SetMultipleStates("OM3",  1);
      options->SetMultipleStates("ODM2", 1);
      options->SetMultipleStates("ODM3", 1);
      */
}


function Benzene($options)
{
  if ($options->IsSetPointCharges())
    $options->numatm = 8;

  if ($options->IsSetGUGACI())
  {
    $options->ici1 = 15;
    $options->ici2 = 15;
  }

  $options->PrintOptions();
?>

 Benzene, standard bond lengths (D6h)
 XX   0.00  0      0.00  0      0.00  0      0  0  0
 XX   1.00  0      0.00  0      0.00  0      1  0  0
 C    1.40  1     90.00  0      0.00  0      1  2  0
 C    1.40  0     90.00  0     60.00  0      1  2  3
 C    1.40  0     90.00  0    -60.00  0      1  2  3
 C    1.40  0     90.00  0    120.00  0      1  2  3
 C    1.40  0     90.00  0   -120.00  0      1  2  3
 C    1.40  0     90.00  0    180.00  0      1  2  3
 H    2.48  1     90.00  0      0.00  0      1  2  3
 H    2.48  0     90.00  0     60.00  0      1  2  3
 H    2.48  0     90.00  0    -60.00  0      1  2  3
 H    2.48  0     90.00  0    120.00  0      1  2  3
 H    2.48  0     90.00  0   -120.00  0      1  2  3
 H    2.48  0     90.00  0    180.00  0      1  2  3

 3  1  4  5  6  7  8
 9  1 10 11 12 13 14
 0  0  0
<?php
  if ($options->mciref === 1)
    $options->PrintReferenceConfigurations();

  if ($options->IsSetPointCharges())
  {
    PrintInPlanePointChargesYZ(-0.1, 5.5, 6);
    PrintAxialPointCharges(0.3, 3.0);
  }
}


function CyclopentadienylAnion($options)
{
  if ($options->IsSetPointCharges())
    $options->numatm = 5;

  if ($options->IsSetGUGACI())
  {
    $options->ici1 = 13;
    $options->ici2 = 12;
  }

  $options->kharge = -1;
  $options->PrintOptions();
?>

 Cyclopentadienyl anion (D5h)
 XX   0.00  0      0.00  0      0.00  0      0  0  0
 XX   1.00  0      0.00  0      0.00  0      1  0  0
 C    1.19  1     90.00  0      0.00  0      1  2  0
 C    1.19  0     90.00  0     72.00  0      1  2  3
 C    1.19  0     90.00  0    -72.00  0      1  2  3
 C    1.19  0     90.00  0    144.00  0      1  2  3
 C    1.19  0     90.00  0   -144.00  0      1  2  3
 H    2.27  1     90.00  0      0.00  0      1  2  3
 H    2.27  0     90.00  0     72.00  0      1  2  3
 H    2.27  0     90.00  0    -72.00  0      1  2  3
 H    2.27  0     90.00  0    144.00  0      1  2  3
 H    2.27  0     90.00  0   -144.00  0      1  2  3

 3  1  4  5  6  7
 8  1  9 10 11 12
 0  0  0
<?php
  if ($options->mciref === 1)
    $options->PrintReferenceConfigurations();

  if ($options->IsSetPointCharges())
    PrintInPlanePointChargesYZ(0.1, 5.5, 5);
}


function CycloheptadienylCation($options)
{
  if ($options->IsSetPointCharges())
    $options->numatm = 7;

  if ($options->IsSetGUGACI())
  {
    $options->ici1 = 17;
    $options->ici2 = 18;
  }

  $options->kharge =  1;
  $options->PrintOptions();
?>

 Cycloheptadienyl cation (D7h)
 XX   0.00  0      0.00  0      0.00000000000000000  0      0  0  0
 XX   1.00  0      0.00  0      0.00000000000000000  0      1  0  0
 C    1.61  1     90.00  0      0.00000000000000000  0      1  2  0
 C    1.61  0     90.00  0     51.42857142857142857  0      1  2  3
 C    1.61  0     90.00  0    -51.42857142857142857  0      1  2  3
 C    1.61  0     90.00  0    102.85714285714285714  0      1  2  3
 C    1.61  0     90.00  0   -102.85714285714285714  0      1  2  3
 C    1.61  0     90.00  0    154.28571428571428571  0      1  2  3
 C    1.61  0     90.00  0   -154.28571428571428571  0      1  2  3
 H    2.69  1     90.00  0      0.00000000000000000  0      1  2  3
 H    2.69  0     90.00  0     51.42857142857142857  0      1  2  3
 H    2.69  0     90.00  0    -51.42857142857142857  0      1  2  3
 H    2.69  0     90.00  0    102.85714285714285714  0      1  2  3
 H    2.69  0     90.00  0   -102.85714285714285714  0      1  2  3
 H    2.69  0     90.00  0    154.28571428571428571  0      1  2  3
 H    2.69  0     90.00  0   -154.28571428571428571  0      1  2  3

  3  1  4  5  6  7  8  9
 10  1 11 12 13 14 15 16
  0  0  0
<?php
  if ($options->mciref === 1)
    $options->PrintReferenceConfigurations();

  if ($options->IsSetPointCharges())
    PrintInPlanePointChargesYZ(-0.1, 6.0, 7);
}


function Oxygen($options)
{
  if ($options->IsSetPointCharges())
    $options->numatm = 2;

  if ($options->IsSetGUGACI())
  {
    $options->ici1 = 7;
    $options->ici2 = 1;
  }

  $options->PrintOptions();
?>

 Oxygen molecule, standard bond length (D0h)
 XX   0.000  0      0.00  0      0.00  0      0  0  0
 XX   1.000  0      0.00  0      0.00  0      1  0  0
 O    0.605  1     90.00  0      0.00  0      1  2  0
 O    0.605  0     90.00  0    180.00  0      1  2  3

 3  1  4
 0  0  0
<?php
  if ($options->mciref === 1)
    $options->PrintReferenceConfigurations();

  if ($options->IsSetPointCharges())
    PrintInPlanePointChargesYZ(0.1, 4.0, 2);
}


function EthyleneD2d($options)
{
  if ($options->IsSetPointCharges())
    $options->numatm = 6;

  if ($options->IsSetGUGACI())
  {
    $options->ici1 = 7;
    $options->ici2 = 5;
  }

  $options->PrintOptions();
?>

 Twisted ethylene, standard bond length (D2d)
 XX   0.00  0      0.00  0      0.00  0      0  0  0
 XX   1.00  0      0.00  0      0.00  0      1  0  0
 C    0.73  1     90.00  0      0.00  0      1  2  0
 C    0.73  0     90.00  0    180.00  0      1  2  3
 H    1.08  1    120.00  1      0.00  0      3  1  2
 H    1.08  0    120.00  0    180.00  0      3  1  2
 H    1.08  0    120.00  0     90.00  0      4  1  2
 H    1.08  0    120.00  0    -90.00  0      4  1  2

 3  1  4
 5  1  6  7  8
 5  2  6  7  8
 0  0  0
<?php
  if ($options->mciref === 1)
    $options->PrintReferenceConfigurations();

  if ($options->IsSetPointCharges())
  {
?>
  0.00   5.00   0.00  -0.2
  0.00  -5.00   0.00  -0.2
  0.00   1.00   3.00   0.1
  0.00   1.00  -3.00   0.1
  3.00  -1.00   0.00   0.1
 -3.00  -1.00   0.00   0.1
<?php
  }
}


function YangsBiradicalTrimmed($options)
{
  if ($options->IsSetPointCharges())
    $options->numatm = 3;

  if ($options->IsSetGUGACI())
  {
    $options->ici1 = 54;
    $options->ici2 = 52;
  }

  $options->PrintOptions();
?>
 https://doi.org/10.1021/ja01508a067 https://doi.org/10.1002/ange.19931050313
 Yang's biradical with tert-butyl substituents replaced by hydrogen atoms (D3)
 C    0.0  0      0.0  0      0.0  0      0    0    0
 C    1.4  1      0.0  0      0.0  0      1    0    0
 C    1.4  0    120.0  0      0.0  0      1    2    0
 C    1.4  0    120.0  0    180.0  0      1    2    3
 C    1.4  1    120.0  1     30.0  1      2    1    3
 C    1.4  0    120.0  0     30.0  0      2    1    4
 C    1.4  0    120.0  0     30.0  0      3    1    2
 C    1.4  0    120.0  0     30.0  0      4    1    2
 C    1.4  0    120.0  0     30.0  0      3    1    4
 C    1.4  0    120.0  0     30.0  0      4    1    3
 C    1.4  1    120.0  1    170.0  1      5    2    1
 C    1.4  0    120.0  0    170.0  0      6    2    1
 C    1.4  0    120.0  0    170.0  0      7    3    1
 C    1.4  0    120.0  0    170.0  0      8    4    1
 C    1.4  0    120.0  0    170.0  0      9    3    1
 C    1.4  0    120.0  0    170.0  0     10    4    1
 C    4.2  1    120.0  0    180.0  0      1    3    4
 C    4.2  0    120.0  0    180.0  0      1    2    3
 C    4.2  0    120.0  0    180.0  0      1    2    4
 O    5.6  1    120.0  0    180.0  0      1    3    4
 O    5.6  0    120.0  0    180.0  0      1    2    3
 O    5.6  0    120.0  0    180.0  0      1    2    4
 H    1.1  1    120.0  1     10.0  1      5    2    1
 H    1.1  0    120.0  0     10.0  0      6    2    1
 H    1.1  0    120.0  0     10.0  0      7    3    1
 H    1.1  0    120.0  0     10.0  0      8    4    1
 H    1.1  0    120.0  0     10.0  0      9    3    1
 H    1.1  0    120.0  0     10.0  0     10    4    1
 H    1.1  1    120.0  1    170.0  1     11    5    2
 H    1.1  0    120.0  0    170.0  0     12    6    2
 H    1.1  0    120.0  0    170.0  0     13    7    3
 H    1.1  0    120.0  0    170.0  0     14    8    4
 H    1.1  0    120.0  0    170.0  0     15    9    3
 H    1.1  0    120.0  0    170.0  0     16   10    4

  2   1   3   4
 17   1  18  19
 20   1  21  22
  5   1   6   7   8   9  10
  5   2   6   7   8   9  10
  5   3   6   7   8   9  10
 11   1  12  13  14  15  16
 11   2  12  13  14  15  16
 11   3  12  13  14  15  16
 23   1  24  25  26  27  28
 23   2  24  25  26  27  28
 23   3  24  25  26  27  28
 29   1  30  31  32  33  34
 29   2  30  31  32  33  34
 29   3  30  31  32  33  34
  0   0   0
<?php
  if ($options->mciref === 1)
    $options->PrintReferenceConfigurations();

  if ($options->IsSetPointCharges())
    PrintInPlanePointChargesXY(0.1, 9.0, 3);
}


function MethylRadical($options)
{
  if ($options->IsSetPointCharges())
    $options->numatm = 1;

  if ($options->IsSetGUGACI())
  {
    $options->ici1 = 4;
    $options->ici2 = 3;
  }

  $options->PrintOptions();
?>

 Methyl radical, standard bond lengths (C3v)
 C    0.00  0      0.00  0      0.00  0      0  0  0
 XX   1.00  0      0.00  0      0.00  0      1  0  0
 H    1.08  1    109.47  1      0.00  0      1  2  0
 H    1.08  0    109.47  0    120.00  0      1  2  3
 H    1.08  0    109.47  0   -120.00  0      1  2  3

 3  1  4  5
 3  2  4  5
 0  0  0
<?php
  if ($options->mciref === 1)
    $options->PrintReferenceConfigurations();

  if ($options->IsSetPointCharges())
  {
?>
 3.0   0.0   0.0   0.1
<?php
  }
}


function NitricOxide($options)
{
  if ($options->IsSetPointCharges())
    $options->numatm = 2;

  if ($options->IsSetGUGACI())
  {
    $options->ici1 = 6;
    $options->ici2 = 2;
  }

  $options->UnsetSymmetry();
  $options->PrintOptions();
?>

 Nitric oxide, experimental bond length (C0v)
 N    0.0000  0      0.00  0      0.00  0      0  0  0
 O    1.1509  1      0.00  0      0.00  0      1  0  0

<?php
  if ($options->mciref === 1)
    $options->PrintReferenceConfigurations();

  if ($options->IsSetPointCharges())
  {
?>
 -3.0   0.0   0.0   -0.1
  4.2   0.0   0.0    0.1
<?php
  }
}


function PrintInPlanePointChargesXY($charge, $radius, $n)
{
  for ($i=0; $i<$n; ++$i)
  {
    $x = $radius * cos(2.0 * M_PI * $i / $n);
    $y = $radius * sin(2.0 * M_PI * $i / $n);
    printf("%20.14f%20.14f%20.14f%10.4f\n", $x, $y, 0.0, $charge);
  }
}


function PrintInPlanePointChargesYZ($charge, $radius, $n)
{
  for ($i=0; $i<$n; ++$i)
  {
    $y = $radius * cos(2.0 * M_PI * $i / $n);
    $z = $radius * sin(2.0 * M_PI * $i / $n);
    printf("%20.14f%20.14f%20.14f%10.4f\n", 0.0, $y, $z, $charge);
  }
}


function PrintAxialPointCharges($charge, $distance)
{
  printf("%20.14f%20.14f%20.14f%10.4f\n",  $distance, 0.0, 0.0, $charge);
  printf("%20.14f%20.14f%20.14f%10.4f\n", -$distance, 0.0, 0.0, $charge);
}
?>
