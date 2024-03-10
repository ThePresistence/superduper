#!/usr/bin/tcsh
#
# BE CAREFUL WHAT YOU ARE DOING NOW USING THIS SCRIPT. 
#
# CUIGL 
#
@ NUM1   =     1
@ LENGTH =  5000 
@ TOT    = 0
@ TOT_R  = 0

rm -f tmp.txt
touch tmp.txt
if ( -e dynvar.in ) then
  set LENGTH = `grep -I "NSTEP" dynvar.in|sed -e 's/[^0-9]//g'`
else
  echo "Could not find dynvar.in in the folder `pwd`"
  echo "Assume length of 5000"
endif
echo "Number of steps should be: $LENGTH"


@ LENGTH = $LENGTH - 2
while ( $NUM1 <= 1000 )  
   set NUMF=`printf '%04i' $NUM1`

  if ( -d run${NUMF}/MD_data/ ) then
    @ TOT = $TOT + 1

    @ NUMber = `tail -1 ./run$NUMF/MD_data/state.dat | gawk '{ print $1 }'`
    if ( $NUMber == $LENGTH ) then
    else
      @ TOT_R = $TOT_R + 1
      echo "run${NUMF}: lenghth $NUMber is less than $LENGTH"
      echo "mv -v ./run$NUMF ./run${NUMF}_bk">>tmp.txt
    endif

  endif

  @ NUM1 = $NUM1 + 1
end

echo "$TOT_R of $TOT too short."
if ($TOT_R > 0) then
  echo "Move them to backup folders? [y/n]"
  set readin = $<
  if ($readin == "y") then
    foreach i ("`cat tmp.txt`")
      $i
    end
endif
rm -f tmp.txt