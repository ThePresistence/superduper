#!/bin/sh

# Create MNDO input files for dynamics run from a Chemshell sampling
# Give the .trj file as first argument
# dynamics.trj -> traj.out, vel.out


if (( "$#" == 1 )) 
then
   filename="$1"
else
   filename="dynamics.trj"
fi
if [ ! -e $filename ] 
then
   echo "File $filename does not exist"
   exit 
fi


NAtoms=`grep "block=coordinates records=" $filename | sed -e 's/[^0-9]//g'`
echo "Found dynamics output with $NAtoms atoms."
echo "Write coordinates to traj.out"

echo "[Molden Format]" > traj.out
echo "[Atoms] Angs" >> traj.out
sed -n '3,'$(($NAtoms+2))' p' $filename| cat -n | sed -e 's/\(.*\)H\(.*\)/ H \1 1\2/' -e 's/\(.*\)C\(.*\)/ C \1 6\2/' -e 's/\(.*\)N\(.*\)/ N \1 7\2/' -e 's/\(.*\)O\(.*\)/ O \1 8\2/' -e 's/\(.*\)F\(.*\)/ F \1 9\2/' | awk -e '{fact=0.52917721;print sprintf("%5s%5d%5d%12.6f%12.6f%12.6f",$1,$2,$3,$4*fact,$5*fact,$6*fact)}' >> traj.out
echo "[GEOMETRIES] XYZ" >> traj.out
sed -e '1 d' -e '/^[0-9-]/ d' -e '/vel/d' -e 's/block.*/'$NAtoms'\n/' $filename | awk -e '$2!=""{fact=0.52917721;print sprintf("%5s%12.6f%12.6f%12.6f",$1,$2*fact,$3*fact,$4*fact)}$2==""{print}' >> traj.out

sed -e '1 d' -e '/^[0-9-]/ d' -e '/vel/d' -e 's/block.*/'$NAtoms'\n/' $filename | awk -e '$2!=""{fact=0.52917721;print sprintf("%5s%12.6f%12.6f%12.6f",$1,$2*fact,$3*fact,$4*fact)}$2==""{print}' > traj.xyz

echo "Write velocities to vel.out"
sed -e '1,'$(($NAtoms+3))' d' -e '/^[HCNOF]/ d' -e '/coord/d' -e 's/block.*//' $filename > vel.out
echo "" >>vel.out

echo "Found $(((`cat vel.out|wc -l`+1)/($NAtoms+1))) geometries."