# PERMEATION
# Martin Voegele, 2018-03-20

# Import necessary modulies
import sys
import argparse
import numpy as np
import scipy as sp
import MDAnalysis as mda


## FUNCTIONS ##


def compute_correlation_via_fft(x, y=None):
    """
    Correlation of two arrays calculated via FFT. 
    If only one is given, the autocorrelation of the first one is calculated.
    """
    x   = np.array(x)
    l   = len(x)
    xft = np.fft.fft(x, 2*l)
    if y is None:
        yft = xft
    else:
        y   = np.array(y)
        yft = np.fft.fft(y, 2*l)
    corr    = np.real(np.fft.ifft(np.conjugate(xft)*yft))
    norm    = l - np.arange(l)
    corr    = corr[:l]/norm
    return corr


def compute_msd_via_correlation(x,y=None,z=None):
    """
    MSD calculated via FFT-based auto-correlation
    """

    y_is_given = y is not None
    z_is_given = z is not None

    # compute correlation
    corrx  = compute_correlation_via_fft(x) 
    # squared distances
    dsq    = x**2
    # repeat correlation for y and z and add  y/z components to squared distance
    if y_is_given:
        corry = compute_correlation_via_fft(y)
        dsq  += y**2
    if z_is_given:
        corrz = compute_correlation_via_fft(z)
        dsq  += z**2
    # sum up all squared distances
    sumsq  = 2*np.sum(dsq)

    # initialize MSD
    nt     = len(x)
    msd    = np.empty(nt)
    msd[0] = 0
    
    # iterate through the time to subtract squared distances
    for m in xrange(1,nt):
        sumsq  = sumsq - dsq[m-1]-dsq[nt-m]
        msd[m] = sumsq/(nt-m)

    # subtract correlations
    msd[1:] -= 2*corrx[1:]    
    if y_is_given:
        msd[1:] -= 2*corry[1:]
    if z_is_given:
        msd[1:] -= 2*corrz[1:]

    return msd



## MAIN ##

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument( '-p',    dest='prmfile',   default='system.prmtop', help="topology file (prmtop)"           )
parser.add_argument( '-t',    dest='trjfile',   default='trajectory.nc', help="trajectory file",       nargs='+' )
parser.add_argument( '-o',    dest='outfile',   default='permeation.dat',help="output file"                      )
parser.add_argument( '-sel',  dest='selection', default='resname WAT and name O', help="selection command"      )
parser.add_argument( '-ref',  dest='reference', default='resname CNT',   help="reference selection command"      )
args = parser.parse_args()

vdw_r_c = 1.6  # vdW radius of Carbon [Angstrom]
vol_h2o = 29.7 # volume of a water molecule [Angstrom^3]  

# Create the universe
u = mda.Universe(args.prmfile,args.trjfile)

# select water and CNT
water = u.select_atoms(args.selection)
cnt   = u.select_atoms(args.reference)

# Initialize all loop variables 
num_water   = len(water)
num_steps   = len(u.trajectory)
t           = np.zeros(num_steps)
dn          = np.zeros(num_steps)
n           = np.zeros(num_steps)
num_inside  = np.zeros(num_steps)
volume_cnt  = np.zeros(num_steps)
all_z_water = np.zeros([num_water,num_steps],dtype=float)
state_water = np.zeros([num_water,num_steps],dtype=int)

# Loop over the whole trajectory
for ts in u.trajectory:

    # simulation time in ps
    t[ts.frame] = ts.time
    
    # Wrap all coordinates into a box with the CNT at the center and origin
    com_cnt_old = cnt.center_of_mass(pbc=False)
    u.atoms.translate(-com_cnt_old+ts.dimensions[:3]*.5)
    u.atoms.wrap()
    com_cnt_new = cnt.center_of_mass()
    u.atoms.translate(-com_cnt_new)    
    
    # Check for each water particle whether it is above or below the membrane
    upper_water = water.positions[:,2] > 0
    lower_water = water.positions[:,2] < 0
    state_water[ upper_water, ts.frame ] =  1
    state_water[ lower_water, ts.frame ] = -1
    
    # Calculate the principal axis of the CNT
    pa = cnt.principal_axes()
    # ... and make sure n_z always points upwards
    if pa[0][2] < 0: 
        pa = -pa 
        
    # Transform everything to the reference frame of the CNT
    tm = np.zeros([4,4])
    tm[:3,:3] = pa
    u.atoms.transform(tm)

    # Calculate z for each water particle
    all_z_water[ :, ts.frame ] = water.positions[:,0]
    
    # Calculate the geometry of the CNT
    min_cnt = np.min(cnt.positions[:,0])
    max_cnt = np.max(cnt.positions[:,0])
    len_cnt = max_cnt - min_cnt
    r2_cnt  = np.average(cnt.positions[:,1]**2+cnt.positions[:,2]**2)
    volume_cnt[ts.frame] = np.pi*(np.sqrt(r2_cnt)-vdw_r_c)**2*len_cnt

    # Check for each water particle whether it is inside the CNT
    abovemin = water.positions[:,0] > min_cnt
    belowmax = water.positions[:,0] < max_cnt
    incylind = water.positions[:,1]**2 + water.positions[:,2]**2 < r2_cnt
    inside_water =  abovemin*belowmax*incylind
    num_inside[ts.frame] = np.sum(inside_water)
    state_water[ inside_water, ts.frame ] = 0
    
    # Calculate n 
    if ts.frame > 0:
        inside_water_before = state_water[:,ts.frame-1] == 0
        dz = all_z_water[ inside_water*inside_water_before, ts.frame ] - all_z_water[ inside_water*inside_water_before, ts.frame-1 ]
        dn[ts.frame] = np.sum(dz)/len_cnt
        n[ts.frame]  = n[ts.frame-1] + dn[ts.frame] # np.sum(dn[:ts.frame+1])
        
# Write the results from the loop
header = "# time [ns], dn, n, # particles inside CNT, volume of CNT [A^3]"
np.savetxt(args.outfile,np.array([t*0.001,dn,n,num_inside,volume_cnt]).T,header=header,fmt=5*'%.8e ')



# Calculate the diffusion coefficient of n (in ps^-1)
msdn  = compute_msd_via_correlation(n)
par   = np.polyfit(t[0:50]-t[0],msdn[0:50],1,cov=True)
diffc = 0.5*par[0][0]
sigma_diffc = 0.5*par[1][0,0]
print "D_n:", diffc, "+/-", sigma_diffc, "ps^-1"

# Permeabilty (A^3/ps)
perm = diffc*vol_h2o
sigma_perm = sigma_diffc*vol_h2o
print "p_f:", perm, "+/-", sigma_perm, "A^3/ps"
