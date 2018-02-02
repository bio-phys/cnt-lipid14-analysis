# Order14
# Martin Voegele, 2016-05-03


## Import necessary modules
import argparse
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import MDAnalysis as mda



## Residues in Lipid14
# 
## Acyl Chains -
# Lauroyl (12:0)     LA
# Myristoyl (14:0)   MY
# Palmitoyl (16:0)   PA
# Stearoyl (18:0)    ST
# Oleoyl (18:1 n-9)  OL
# 
## Head Groups 
# Phosphatidylcholine      PC
# Phosphatidylethanolamine PE

n_carbon = {'LA': 12,'MY':14,'PA':16,'ST':18,'OL':18}



## Class and function definitions

class Tail(object):
    """Contains information about all residues of a specific lipid tail type.""" 
    
    def __init__(self,name,tail_length):
        
        # Name and length of the tail
        self.name     = name
        self.length   = tail_length
        # All atoms belonging to the tail
        self.atoms    = u.select_atoms('resname '+name)
        # Lists of C and H atom names for selection
        self.list_c   = ["name C1{}".format(i) for i in range(2,self.length)]
        self.list_hr  = ["name H{}R".format(i) for i in range(2,self.length)]
        self.list_hs  = ["name H{}S".format(i) for i in range(2,self.length)]
        # C and H atoms
        self.atoms_c  = [self.atoms.select_atoms(selection) for selection in self.list_c]
        self.atoms_hr = [self.atoms.select_atoms(selection) for selection in self.list_hr]
        self.atoms_hs = [self.atoms.select_atoms(selection) for selection in self.list_hs]
        # Number of carbon atoms in one tail
        self.tail_length = len(self.atoms_c)
        # Number of lipid tails of the chosen species
        self.n_residues  = len(self.atoms_c[0])
        if self.n_residues == 0:
            print "I found no residues of type "+self.name+". No order parameter is calculated."
            self.order_was_calculated = False
            self.order, self.order_t, self.order_mol = None, None, None
        else:
            print "I found "+str(self.n_residues)+" residues of type "+self.name+"."
            # Order parameters: 
            # 1. average, 
            # 2. time series, 
            # 3. separate for each molecule 
            self.order, self.order_t, self.order_mol = self.calculate_order_cd()
            self.order_was_calculated = True
        
        
    def calculate_s(self,a,b,normal):
        """Calculates the order parameter for atoms in two arrays of equal length according to the membrane normal (0,1,2) = (x,y,z)"""
    
        diff = a.positions - b.positions
        dx   = diff[:,normal]
        l    = np.sqrt(np.sum(np.power(diff,2),axis=-1))

        return 1.5*np.square(dx/l)-0.5
        
    
    def calculate_order_cd(self):
        """Calculates the order parameter from distances between tail carbon atoms and the corresponding hydrogen atoms."""

        # Initialize time average for every single tail
        order_z = np.zeros([self.n_residues,self.tail_length]) 
        # Initialize ensemble average for every time step
        order_t = np.zeros([u.trajectory.n_frames,self.tail_length])              

        # Loop through the trajectory
        for t,ts in enumerate(u.trajectory):
            # Loop through the carbon number (for all molecules at once)
            for i in xrange(self.tail_length):
                # Calculate order parameters for a certain carbon number
                sz = self.calculate_s(self.atoms_c[i],self.atoms_hr[i],2)
                # Update lists
                order_z[:,i] += sz
                order_t[t,i] = np.mean(sz)

        # Normalize the time averaged order parameter
        order_z /= u.trajectory.n_frames

        # total average over time and ensemble
        ordercd = np.mean(order_t,axis=0)

        return ordercd, order_t, order_z
    
    

## MAIN ##

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-p', dest='prmfile', default='system.prmtop', help="topology file (prmtop)")
parser.add_argument('-t', dest='trjfile', default='trajectory.nc', help="trajectory file")
parser.add_argument('-o', dest='outfile', default='./system',      help="directory and prefix for output files")
parser.add_argument('-r', dest='resname', default='PA', nargs='+', help="residue name of the lipid tails", choices=['LA', 'MY', 'PA', 'ST', 'OL'])
args = parser.parse_args()
#args = parser.parse_args('-p ../simulation/system.prmtop -t ../simulation/05_Prod.nc -o ./system -r PA OL LA'.split())


# Create universe from files
u = mda.Universe(args.prmfile,args.trjfile)

# Loop over all lipid tail species 
for rn in args.resname:
    
    # Construct the tails and calculate all values
    num  = n_carbon[rn]-2
    tail = Tail(rn,n_carbon[rn])
    
    if tail.order_was_calculated:
        # Save the results to files
        header = "Order parameters for each tail carbon (time and ensemble average)"
        np.savetxt( args.outfile+'_order_'+rn+'.dat',     tail.order,     fmt='%.8e '    , header=header )
        header = "Order parameters for each tail carbon over time (ensemble average)"
        np.savetxt( args.outfile+'_order_'+rn+'_t.dat',   tail.order_t,   fmt=num*'%.8e ', header=header )
        header = "Order parameters for each tail carbon over different residues (time average)"
        np.savetxt( args.outfile+'_order_'+rn+'_mol.dat', tail.order_mol, fmt=num*'%.8e ', header=header )
        

exit
