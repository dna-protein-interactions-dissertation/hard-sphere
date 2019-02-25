# hard-sphere
We simulate a protein-sized sphere moving under the influence of a hard collisions potential in a heat bath of point particles. 

The file hardpot.f90 is the central program and contains 3 functions which are called to aid the sampling of distributions and calculation of square distance.

The parameters have been chosen to be a realistic representation of RNA polymerase but can easily be changed to other proteins, providing p_new remains small.
