# RDF14
# Martin Voegele, 2016-07-29

# Import necessary modulies
import sys
import argparse
import numpy as np
import scipy as sp
import MDAnalysis as mda
import MDAnalysis.analysis.helanal
import geometry



def hist_from_file(prmfile,trjfile,selection,reference,skip,a0,maxrange,nbins):

    # Create universe from files
    u = mda.Universe(prmfile,trjfile)

    # Initialize box dimensions
    boxdim = np.zeros(6)

    # Initialize empty histogram
    dn = np.copy(a0[0])

    # counting variables
    i        = 0
    area     = 0
    num_dist = 0

    for ts in u.trajectory[::int(skip)]:
        
        i += 1
       
        # Bring the reference into a nice region close to the center of the box
        recenter = u.select_atoms(reference).positions[0]
        u.select_atoms('all').positions = u.select_atoms('all').positions - [recenter[0],recenter[1],0]
        # Bring the reference to the center of the box
        recenter = u.select_atoms(reference).center_of_mass(pbc=False)
        u.select_atoms('all').positions = u.select_atoms('all').positions - recenter
        # Get the new coordinates of the selection
        coord    = u.select_atoms(selection).positions

        # Calculate the reference coordinate system
        refcs = u.select_atoms(reference).principal_axes(pbc=False)
        nx    = refcs[2]
        ny    = refcs[1]
        nz    = refcs[0]
        
        # Get the box dimensions
        dim   = ts.dimensions[:3] 
        
        # Correct for PBC
        newcoord = ( coord[:,] + dim/2 )%dim - dim/2

        # Initialize list of distances
        rlist = []

        # rotate to the principal axis system
        for point in newcoord:
            px = np.dot(point,nx)
            py = np.dot(point,ny)
            rlist.append(np.sqrt(px**2+py**2))

        rlist = np.array(rlist)
        
        # Make histogram for this timestep    
        curr_hist  = np.histogram(rlist, range=[0,maxrange], bins=nbins)
        
        # Append histogram to total histogram
        dn += curr_hist[0]
	
        # number of particles
        num_dist += len(rlist)

        # Sum of the areas of the box
        area += dim[0]*dim[1]

    return dn, num_dist, area, i



## MAIN ##

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument( '-p',    dest='prmfile',   default='system.prmtop', help="topology file (prmtop)"           )
parser.add_argument( '-t',    dest='trjfile',   default='trajectory.nc', help="trajectory file",       nargs='+' )
parser.add_argument( '-o',    dest='outfile',   default='rdf.dat',       help="output file"                      )
parser.add_argument( '-sel',  dest='selection', default='resname OL PA and name C*', help="selection command"    )
parser.add_argument( '-ref',  dest='reference', default='resname CNT',   help="reference selection command"      )
parser.add_argument( '-skip', dest='skip',      default='1',             help="use every nth frame"              )
args = parser.parse_args()


# Maximal range of the RDF
maxrange = 35
nbins    = 350 

# Make empty histogram for initialization
a0 = np.histogram([maxrange+1], range=[0,maxrange], bins=nbins)

# initialize the total histogram
tot_hist = np.copy(a0[0])

# r axis for the RDF
r  = a0[1][:-1]
dr = r[1]-r[0]
r  = r + dr/2

# Initialize box properties
tot_area   = 0
tot_ndist  = 0
tot_frames = 0

# Go through all parts of the trajectory
for trjf in args.trjfile:
    print trjf
    dn, num_dist, area, frames = hist_from_file(args.prmfile,trjf,args.selection,args.reference,args.skip,a0,maxrange,nbins)
    tot_hist   += dn
    tot_area   += area
    tot_ndist  += num_dist
    tot_frames += frames

# Calculate the RDF from the histogram
rho = tot_ndist/tot_area
g   = tot_hist/(2*rho*np.pi*dr*r)/tot_frames


output = np.transpose([r, g])
np.savetxt(args.outfile,output,fmt='%.18e %.18e')

