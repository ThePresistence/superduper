import mlatom as ml

# Get the initial guess for the molecules to optimize
initmol = ml.data.molecule.from_xyz_file('init.xyz')

# Choose one of the predifined (automatically recognized) methods
mymodel = ml.models.methods(method='AIQM1')

# Optimize the geometry with the choosen optimizer:
geomopt = ml.optimize_geometry(model=mymodel, initial_molecule=initmol, program='ASE', maximum_number_of_steps=10000)

# Get the optimized molecule
final_mol = geomopt.optimized_molecule

# Do vibration analysis and thermochemistry calculation
freq = ml.thermochemistry(model=mymodel, molecule=final_mol, program='ASE')
# or vibration analysis only
# freq = ml.freq(model=mymodel, molecule=final_mol, program='')

# Save the molecule with vibration analysis and thermochemistry results
final_mol.dump(filename='final_mol.json',format='json')

# Check vibration analysis
print("Mode     Frequencies     Reduced masses     Force Constants")
print("           (cm^-1)            (AMU)           (mDyne/A)")
for ii in range(len(final_mol.frequencies)):
    print("%d   %13.4f   %13.4f   %13.4f"%(ii,final_mol.frequencies[ii],final_mol.reduced_masses[ii],final_mol.force_constants[ii]))

# Check thermochemistry results
print(f"Zero-point vibrational energy: {final_mol.ZPE} Hartree")
print(f"Enthalpy at 298 K: {final_mol.H} Hartree")
print(f"Gibbs Free energy at 298 K: {final_mol.G} Hartree")
print(f"Heat of formation at 298 K: {final_mol.DeltaHf298} Hartree")