#!/usr/bin/tcsh

set rootdir = `dirname $0`
set abs_rootdir = `cd $rootdir && pwd`

set mmean_path="$abs_rootdir/mmean.py"

@ NUM1 =   1

rm -rf data.ls

while ( $NUM1 <= 1000 )  
   set NUMF=`printf '%04i' $NUM1`
   if ( -f run$NUMF/MD_data/state.dat ) then
     echo "./run$NUMF/MD_data/state.dat" >> data.ls
   endif

   @ NUM1 = $NUM1 + 1
end
@ NUM1 = $NUM1 - 1
echo "$NUM1"

python $mmean_path -s --s-error

sed -i -e 's/^ *$/\n/' data.mean