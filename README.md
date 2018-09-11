# KitaevHoneycombHigherGenus
MATLAB functions for creating an effective Hamiltonian for the Kitaev Honeycomb model, 
restricted to particular vortex sector of its' Hilbert space, on lattices of genus greater than 1


The function Hamiltonian.m returns the matrix representing the Hamiltonian 
for the Kitaev Honeycomb model on a lattice with a specified genus.

The script Genus_Degeneracy.m dianonalises the Hamiltonian of the model 
in the vortex free sector for a range of couplings beginning in the models' 
Ableian phase and ending in the non-Abelian phase. It then plots the energy 
of the ground state and first excited state from each homology sector as a
function of the coupling.

This code was used to produce the results presented in:
http://iopscience.iop.org/article/10.1088/1367-2630/aabb95
