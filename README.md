# Thermo-visco-elastic BEM

This folder contains the code to generate BEM data for the scattering of two rigid spheres in a linear thermo-visco-elastic medium.

The BEM output can be obtained from the MATLAB script TVEBEM.m which calls the three functions Uij_vector_BIR.m, temp_vector_BIR.m and scalar_rigid_BIR.m. These functions contain the fundamental solutions for displacement and temperature in a thermo-visco-elastic medium which have been obtained using the Mathematica script GreensFunctions.wl. The script TVEBEM.m outputs a .txt file of the traction and thermal gradient at the centre of each element which can be imported into the MATLAB script PlotTVEBEM.m to plot the fields.

This code is adapted from 'Elastodynamic Dual Void Scattering' by Charlotte Charlton that is available at https://doi.org/10.5281/zenodo.19162993 under a CC BY 4.0 license. The scripts share the same underlying structure as the original, with modifications to extend the work from isothermal elasticity to linear thermo-visco-elasticity.
