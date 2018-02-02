import numpy as np
import numpy.linalg as la


def dist_pbc(r1,r2,boxdim):
    """Calculates the distance vector between two points under periodic boundary conditions
    
    :Arguments:
        *r1: (r1x,r1y,r1z) numpy array, one point in space
        *r2: (r2x,r2y,r2z) numpy array, another point in space
        *box: (bx,by,bz)   numpy array, dimensions of (rectangular) periodic box
        
    :Returns:
        *dist: The distance between the points
    """
    
    dist = np.array([0,0,0])
    
    for i, di in enumerate(dist):
        #print i
        dist[i] = r1[i]-r2[i]
        if dist[i] > boxdim[i]/2:
            dist[i] = dist[i]-np.sign(dist[i])*boxdim[i]
    
    return dist


def dist_point_line(p,a,n):
    """Calculates the distance of a point to a line
    
    :Arguments:
        *p: (x0,y0,z0) numpy array, One point in space
        *a,n: numpy array, point and vector defining the line
        
    :Returns:
        *distance: The distance between the point and the line
    """
    
    distance = (a-p)-np.dot((a-p),n)*n
    
    return np.linalg.norm(distance)


def dist_point_line_pbc(p,a,n,boxdim):
    """Calculates the distance of a point to a line
    
    :Arguments:
        *p: (x0,y0,z0) numpy array, One point in space
        *a,n: numpy array, point and vector defining the line
        
    :Returns:
        *distance: The distance between the point and the line
    """
    #print n, a, p, boxdim
    #print dist_pbc(a,p,boxdim)
    #print np.dot(dist_pbc(a,p,boxdim),n)
    #print np.dot(dist_pbc(a,p,boxdim),n)*n
    distance =  dist_pbc(a,p,boxdim) - np.dot(dist_pbc(a,p,boxdim),n)*n
    
    return np.linalg.norm(distance)



def tilt(calpha):
    """
    Purpose:
    calculate the angle between the principal axis of one alpha
    helix and the z-axis
    The principal axis of the helix is calculated using the MDAnalysis
    function
    
    Arguments:
    calpha: alpha atoms object within MDAnalysis
    """
    P1,P2,P3=calpha.principalAxes()
    zaxis=np.array([0,0,1])
    
    """angle"""
    cosang = np.dot(P1, zaxis)
    sinang = la.norm(np.cross(P1, zaxis))
    return np.arctan2(sinang, cosang)