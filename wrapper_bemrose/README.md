This folder contains the 'wrapper_bemrose.m' script used to run the analysis. As the script controls the overall workflow and calls the relevant computations functions contained within 'update_compute.m', it should be run as the main entry point.

## Modifications
The script has been modified and updated from the original 'wrapper.m' script produced by Rushby et al., (2018) to fit the requirements of this current study.

### Stellar Parameters
Stellar parameter values from line 52 to 58 in the script were modified based upon the planet that was being modelled for the long-term climate evolution simulations. Specifically, 'semi_major_axis', 'e', and 'star_size.' 'Planet_mass' on line 15 was also altered for each run.

### Weathering Options
Lines 71, 74 and 75 were altered when analysing the consequences of varying heat flow output on climatic outcomes. 0.7 and 0.2 were used interchangeably for high and low heat flow outcomes, respectively, in which each planet presented in the results was run under both simulations and the subsequent data compared.
As discussed in the paper, preliminary runs varied 'B' and 'B_kinetic' values, however, for the final analysis the values were kept at 0.5 and 11.1 respectively which were originally provides by Berner (1994). The change was largely due to a lack of justification with using the default values provided by Rushby et al., (2018).

## Output
The output of running this script was a 'trimmed_tracker' file which was exported for visual representations. As not all data produced was required, columns 1, 2, 7, and 16 were specifically extracted and processed. The names of all corresponding columns, including those mentioned here, can be found on the 'updated_compute.m' files from line 370.
