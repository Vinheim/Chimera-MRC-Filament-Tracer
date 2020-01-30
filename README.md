# Chimera-MRC-Filament-Tracer
This repository serves as a host for all of the work done on my funded research project, "Precise Termination of Actin Filaments in Hair Cell Stereocilia". 
This program reads data from binary .mrc files interpreted within the Chimera protein structure visualization environment. 
Within the files is stored density information regarding hair cell stereocilia.
The program reads this information and derives membrane and filament structures to be appended onto the base .mrc density map. 
After initial charting of the filaments and membranes, the filaments are pruned at the points at which their computed local average density values fall beneath a threshold, or at which they meet the terminal end of their corresponding membranes.
