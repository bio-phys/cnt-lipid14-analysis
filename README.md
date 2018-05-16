# This project contains

Analysis of Lipid14 simulations
 - order14.py: Calculates the deuterium order parameter of the lipids (works without CNT).

Analysis of Lipid14 simulations with a CNT
 - order14distance.py: Calculates the order parameter as a function of the annular shell around the CNT axis.
 - order14distance-manyCNTs: same as above, but adds a category for lipids shared between two CNTs.
 - cntmotion.py: Calculates the COM motion and the tilting angle of a CNT (works with any lipids).
 - rdf14.py: Calculates the radial distribution of a given atom type aroud the CNT axis.

Required geometry functions
 - geometry.py

#  Requirements: 
 - Python 2.7
 - Python packages: sys, argparse, numpy, scipy, MDAnalysis

# Literature
 - M. Vögele, J. Köfinger, G. Hummer: 
   Simulations of Carbon Nanotube Porins in Lipid Bilayers.
   Faraday Discuss., 2018, Accepted Manuscript, DOI: 10.1039/C8FD00011E  
   http://pubs.rsc.org/en/content/articlelanding/2018/fd/c8fd00011e
