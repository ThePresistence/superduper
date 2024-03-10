#!/bin/bash

NPROC=1
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

path="$DIR/MDextr.py"
echo $path

for ((NUM1=1; NUM1 <= 1000; NUM1++)); do
   NUMF="`printf '%04i' $NUM1`"
   if [ -d run$NUMF ]; then
      if [ \( -d run$NUMF/MD_data \) -o \( -d run${NUMF}_bk/MD_data \) ]; then
         echo
      else
         command="${command}run$NUMF,"
      fi
   fi
done
echo $command
echo $command | xargs --max-procs=$NPROC --delimiter="," --max-args=1 --replace='REPL' sh -c 'cd ./'REPL';\
rm -rf MD_data;\
echo "'REPL':";\
'"$path"' > /dev/null;\
cd ..'

