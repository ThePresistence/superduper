#!/usr/bin/tcsh
#
# REMOVE TRAJECTORIES THAT ARE NOT CONTINUOUS WITH RESPECT TO POTENTIAL ENERGIES BETWEEN FORWARD AND BACKWARD STEPS
#
# BE CAREFUL WHAT YOU ARE DOING NOW USING THIS SCRIPT. 
#      
# CUIGL
# 
set rootdir = `dirname $0`
set abs_rootdir = `cd $rootdir && pwd`
set path_prog = "$abs_rootdir/TOOLS"

@ NUM1  = 1
@ TOT   = 0
@ TOT_R = 0

if ( "$1" == "" ) then
   set CUTOFF = "10"
else 
   set CUTOFF = "$1"
endif  


sed -i -e 's/CUTOFF=.*/CUTOFF='$CUTOFF'/' $path_prog/energy_difference_forward_backward.f90
grep "CUTOFF=" $path_prog/energy_difference_forward_backward.f90
ifort -o $path_prog/energy_difference_forward_backward.exe $path_prog/energy_difference_forward_backward.f90

rm -f tmp.txt
touch tmp.txt

while ( $NUM1 <= 1000 ) 
   set NUMF=`printf '%04i' $NUM1`
   
   if ( -d run$NUMF/MD_data/ ) then 
     @ TOT = $TOT + 1
     cd run$NUMF
     rm -rf bad_run_pot
     rm -rf bad_run_tot
     $path_prog/energy_difference_forward_backward.exe
     cd ..
     if ( -e ./run$NUMF/bad_run_tot || -e ./run$NUMF/bad_run_pot ) then
       echo "mv -v ./run$NUMF ./run${NUMF}_bk2">>tmp.txt
       @ TOT_R = $TOT_R + 1
       echo "run${NUMF}: PES crash"
     endif
   endif

   @ NUM1 = $NUM1 + 1

end


echo "$TOT_R of $TOT bad energy difference."
if ($TOT_R > 0) then
  echo "Move them to backup folders? [y/n]"
  set readin = $<
  if ($readin == "y") then
    foreach i ("`cat tmp.txt`")
      $i
    end
endif

rm -f tmp.txt
