# CNTmotion14
# Martin Voegele, 2016-08-03

# Import necessary modulies
import sys
import argparse
import numpy as np
import scipy as sp
import MDAnalysis as mda
import MDAnalysis.analysis.helanal
import geometry



def motion_from_file(prmfile,trjfile,selection,skip):

    # Create universe from files
    u = mda.Universe(prmfile,trjfile)

    sim_time = []
    tube_com = []
    tilt_cos = []

    for ts in u.trajectory[::skip]:
        sim_time.append(ts.time)
        cnt = u.select_atoms(selection)
        com = cnt.center_of_mass()
        tube_com.append( com )
        princ_ax = cnt.principal_axes()[0]
        tilt_cos.append( np.abs(princ_ax[2]) )

    return sim_time, tube_com, tilt_cos



## MAIN ##

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument( '-p',    dest='prmfile',   default='system.prmtop', help="topology file (prmtop)"           )
parser.add_argument( '-t',    dest='trjfile',   default='trajectory.nc', help="trajectory file",       nargs='+' )
parser.add_argument( '-ocom', dest='outfcom',   default='com.dat',       help="output file for COM"              )
parser.add_argument( '-otil', dest='outftil',   default='tilt.dat',      help="output file for tilting cos"      )
parser.add_argument( '-sel',  dest='selection', default='resname CNT',   help="selection command"                )
parser.add_argument( '-skip', dest='skip',      default=1, type=int,     help="use every nth frame"              )
args = parser.parse_args()


# Initialize the time, the COM, and the tilting cosine
tot_time = []
tot_tilt = []
tot_com  = []

# Go through all parts of the trajectory
for trjf in args.trjfile:
    print trjf
    time, com, tilt = motion_from_file(args.prmfile,trjf,args.selection,args.skip)
    tot_com  += com
    tot_tilt += tilt
    tot_time += time

tot_com  = np.array(tot_com)
tot_tilt = np.array(tot_tilt)
tot_time = np.array(tot_time)


output = np.transpose([tot_time, tot_com[:,0], tot_com[:,1], tot_com[:,2]])
np.savetxt(args.outfcom,output,fmt='%.18e %.18e %.18e %.18e')

output = np.transpose([tot_time, tot_tilt])
np.savetxt(args.outftil,output,fmt='%.18e %.18e')


