# Order14
# Martin Voegele, 2016-05-03


## Import necessary modules
import argparse
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import MDAnalysis as mda
import geometry


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
    
    def __init__(self,name,tail_length,reference):
        
        # Name and length of the tail
        self.name     = name
        self.length   = tail_length
        self.ref      = reference
        # Definitions of the four closest shells 
        self.d_range  = [12,16,20,24,28]
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
            self.order = None
        else:
            print "I found "+str(self.n_residues)+" residues of type "+self.name+"."
            self.order = self.calculate_order_cd()
            self.order_was_calculated = True
        
        
    def calculate_s(self,a,b,normal):
        """Calculates the order parameter for atoms in two arrays of equal length according to the membrane normal (0,1,2) = (x,y,z)"""
    
        diff = a.positions - b.positions
        dx   = diff[:,normal]
        l    = np.sqrt(np.sum(np.power(diff,2),axis=-1))

        return 1.5*np.square(dx/l)-0.5
        
    
    def calculate_order_cd(self):
        """Calculates the order parameter from distances between tail carbon atoms and the corresponding hydrogen atoms."""

        # Initialize the list for distance ranges
        order_d = [[],[],[],[],[],[],[]]
        
        # Select the carbon nanotube
        cnt     = u.select_atoms(self.ref)

        # Loop through the trajectory
        for ts in u.trajectory[::10]:
            
            # Box dimensions
            box_dim = ts.dimensions            
            
            # Calculate COM and principal axis of the CNT
            com_cnt = []
            pax_cnt = []
            for i,cnti in enumerate(cnt.residues):
                cnta = cnti.atoms
                com_cnt.append( cnta.center_of_mass() )
                pax_cnt.append( cnta.principal_axes(pbc=True)[0] )
            
            # Calculate the order parameter of each single tail and assign it to the respective distance range
            for i,res in enumerate(self.atoms.residues):
                
                # Select C and H atoms
                atoms_c  = [res.atoms.select_atoms(selection) for selection in self.list_c]
                atoms_hr = [res.atoms.select_atoms(selection) for selection in self.list_hr]
                atoms_hs = [res.atoms.select_atoms(selection) for selection in self.list_hs]
                
                # Loop through the carbon number and calculate the order parameter
                sz = []
                for j in xrange(self.tail_length):
                    # Calculate order parameters for a certain carbon number
                    sz.append(self.calculate_s(atoms_c[j],atoms_hr[j],2))
                order = np.mean(np.array(sz))
                
                # Calculate the distance to the CNTs
                com_lip = res.atoms.center_of_mass()
                
                range_index = len(order_d)-1
                
                for i,cnti in enumerate(cnt.residues):
                    cnta = cnti.atoms
                    dist = geometry.dist_point_line_pbc(com_lip,com_cnt[i],pax_cnt[i],box_dim)
                
                    # Append the order parameter to the respective distance range
                    if dist < self.d_range[0]:
                        if range_index == 1:
                            range_index = 0
                        elif range_index == 0:
                            range_index = 0
                        else:
                            range_index = 1
                    elif dist < self.d_range[1]:
                        range_index = np.minimum(2,range_index)
                    elif dist < self.d_range[2]:
                        range_index = np.minimum(3,range_index)
                    elif dist < self.d_range[3]:
                        range_index = np.minimum(4,range_index)
                    elif dist < self.d_range[4]:
                        range_index = np.minimum(5,range_index)
                    else:
                        range_index = np.minimum(6,range_index)
                        
                order_d[range_index].append(order)
            
        mean_order = []
        std_order  = []
        sem_order  = []
        for l in order_d:
            if len(l) > 0:
                mean_order.append(np.mean(l))
                std_order.append(np.std(l))
                sem_order.append(np.std(l)/len(l))
            else:
                mean_order.append(-1)
                std_order.append(0)
                sem_order.append(0)    
                
        print np.array(mean_order)
        print np.array(std_order)
        print np.array(sem_order)
            
        return np.array(mean_order), np.array(std_order), np.array(sem_order)


    

## MAIN ##

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-p', dest='prmfile', default='system.prmtop', help="topology file (prmtop)")
parser.add_argument('-t', dest='trjfile', default='trajectory.nc', help="trajectory file")
parser.add_argument('-o', dest='outfile', default='./system',      help="directory and prefix for output files")
parser.add_argument('-r', dest='resname', default='PA', nargs='+', help="residue name of the lipid tails", choices=['LA', 'MY', 'PA', 'ST', 'OL'])
parser.add_argument('-ref', dest='reference', default='resname CNT',   help="reference to which the distance is calculated")
args = parser.parse_args()
#args = parser.parse_args('-p ../simulation/system.prmtop -t ../simulation/05_Prod.nc -o ./system -r PA OL LA'.split())


# Create universe from files
u = mda.Universe(args.prmfile,args.trjfile)

# Loop over all lipid tail species 
for rn in args.resname:
    
    # Construct the tails and calculate all values
    num  = n_carbon[rn]-2
    tail = Tail(rn,n_carbon[rn],args.reference)
    
    if tail.order_was_calculated:
        # Save the results to files
        header = "Order parameters in the ranges between "+str(tail.d_range[0])+", "+str(tail.d_range[1])+", "+str(tail.d_range[2])+", and "+str(tail.d_range[3])+" Angstrom with Stdev. and SEM"
        np.savetxt( args.outfile+'_order_'+rn+'.dat',     tail.order,     fmt='%.8e '    , header=header )
        
exit
