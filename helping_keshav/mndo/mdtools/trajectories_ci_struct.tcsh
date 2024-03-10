#!/usr/bin/tcsh
#
# Create xyz-files for all starting, hopping and ending points
#

set NUMATOMS   = `head -1 filtered-sample.xyz`
echo "CI for $NUMATOMS atoms "
set rootdir = `dirname $0`
set abs_rootdir = `cd $rootdir && pwd`

set path_prog = "$abs_rootdir/TOOLS"

ifort -o $path_prog/ci_struct.exe $path_prog/ci_struct.f90

@ NUM1   = 1
@ TOT    = 0
@ TOT_CI = 0
set nonomatch
rm -rf combined_ci_struct*
rm -rf combined_final_structures.xyz
rm -rf combined_start_structures.xyz
rm -rf combined_3structures.xyz

@ NUMTAIL = $NUMATOMS + 2 
@ NUMLOW = $NUMATOMS + 4 
@ NUMHIGH = $NUMATOMS * 2 + 5 
@ NUMXYZ = $NUMATOMS + 2 
while ( $NUM1 <= 1000 ) 
  set NUMF=`printf '%04i' $NUM1`

  if ( -d run$NUMF ) then 
    @ TOT = $TOT + 1
    cd run$NUMF
    rm -f ci_struct*
    $path_prog/ci_struct.exe
    cd ..
    cat ./run$NUMF/traj.out|sed -n "$NUMLOW,$NUMHIGH p" >> combined_start_structures.xyz
    cat ./run$NUMF/traj.out|sed -n "$NUMLOW,$NUMHIGH p" > ./run$NUMF/3structures.xyz
    @ CI_EXISTS = 0
    foreach x (./run$NUMF/ci_struct* )
      if ( `cat $x|wc -l` == $NUMXYZ ) then 
        @ CI_EXISTS = 1
	echo "$x created"
	@ TOT_CI = $TOT_CI + 1
        cat $x >> combined_`basename $x`
        cat $x >> ./run$NUMF/3structures.xyz
      endif
    end
    cat ./run$NUMF/traj.out|tail -n $NUMTAIL -q >> combined_final_structures.xyz
    cat ./run$NUMF/traj.out|tail -n $NUMTAIL -q >> ./run$NUMF/3structures.xyz
    if ( $CI_EXISTS == 1 ) then
      cat ./run$NUMF/3structures.xyz >> combined_3structures.xyz
    endif
  endif

  @ NUM1 = $NUM1 + 1

end

foreach x (./combined_ci_struct*.xyz )
  set state1=`basename $x .xyz|awk -F _ -e '{print $4}'`
  set state2=`basename $x .xyz|awk -F _ -e '{print $5}'`
  sed -n '2~'$NUMXYZ'p' combined_ci_struct_${state1}_${state2}.xyz > IC_TIME_${state1}_${state2}.data
end

echo "@@@ TOTALLY $TOT_CI STRUCTURES FROM $TOT @ @@"

