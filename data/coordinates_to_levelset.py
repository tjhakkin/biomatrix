#
# Given input data folder, forms a level set function containing all the shapes
# in the given data folder (.txt files), in a domain such that the longest axis 
# has given number of nodes within [-1,1]. 
# 
# The domain includes padding to keep the interface at given distance (in %)
# from the boundaries.
#
# Computes only the initial level set function, i.e., direct Euclidean distances
# to mesh nodes from mesh edge crossing interface nodes. No projection distances
# to interface edges.
#

import glob, os, sys
import numpy as np
import scipy.io as sio
from scipy.spatial.distance import cdist
from scipy.spatial import Delaunay
from shapely.geometry.polygon import LineString, Polygon


nargs = len(sys.argv)
if (nargs != 4 and nargs != 5):
    print("Usage: python [script] [data folder] [shortest axis nodes] [pad] "
          "[crop y (optional)]")
    exit()

min_nodes = int(sys.argv[2])    # number of nodes in [-1,1] (the shortest axis)
pad = float(sys.argv[3])        # padding on all sides as (pad*100)% 
cropy = []
if (nargs == 5):
    cropy = float(sys.argv[4])

os.chdir(sys.argv[1])
files = glob.glob("*.txt")
data = []
for file in files:
    d = np.loadtxt(fname=file, delimiter='\t', skiprows=0, usecols=[0,1])
    data.append(d)

# Scale data so that the longest axis fits within [-1,1] with padding.
dcon = np.concatenate(data)
smin = np.min(dcon, axis=0)                     # lower left corner of shapes
smax = np.max(dcon, axis=0)                     # top right corner
dscaled = (dcon - smin) / max(smax - smin)      # [0,1]^2
dscaled = 2*dscaled - np.max(dscaled, axis=0)   # [-1,-1]^2
dscaled = dscaled * (1-pad)                     # [pad-1,1-pad]^2

# Compute physical domain dimensions and the number of nodes.
dmax = np.max(dscaled, axis=0) + pad
dmax = np.round(dmax, 1)
dmin = np.min(dscaled, axis=0) - pad
dmin = np.round(dmin, 1)
nodes = (dmax - dmin) * min_nodes / 2
print("Domain dimensions: %d x %d" %(nodes[0], nodes[1]))

# Generate mesh nodes.
x = np.linspace(dmin[0], dmax[0], nodes[0])
y = np.linspace(dmin[1], dmax[1], nodes[1])
px, py = np.meshgrid(x, y)
p = np.array([px.flatten(), py.flatten()]).transpose()


print("1/3 - triangulating, computing initial distances...")
t = Delaunay(p)
t = t.simplices
d = cdist(p, dscaled)
phi = np.min(d, axis=1)
print("done.\n")


#
# Given initial distances, find nodes in a neighborhood of the interface & get
# the associated edges from t. For each edge, intersect with the interface; 
# set the intersection point as a temporary (interface) node for new distance
# computations.
#
print("2/3 - finding interface nodes along mesh edges...")

# List of edges as node pairs.
a = np.concatenate( (t[:,[0,1]], t[:,[1,2]], t[:,[2,0]]) )
a.sort()
# Get unique node pairs via byte representation.
atype = np.dtype( (np.void, a.dtype.itemsize * a.shape[1]) )
b = np.ascontiguousarray(a).view(atype)
_, idx, t2e = np.unique(b, return_index = True, return_inverse = True)
edges = a[idx]          # edge -> pair of nodes.

# Find edges that are close to the interface.
d = phi[np.unique(t)]
bwidth = 0.1                                # TODO: no fixed bandwidth!
nodes_close = np.where(d < bwidth)
# np.in1d() returns true/false. Multiply by 1 to get 0/1 for summation:
edges_close = np.in1d(edges[:,0], nodes_close)*1 + np.in1d(edges[:,1], nodes_close)*1
edges_close = np.where(edges_close == 2)[0]

# shape start indices in the full data
idx = [0]
for d in data:
    idx.append(idx[-1] + np.size(d,0))

#
# Interface node positions along mesh edges
#
ifp = []                # interface node positions
i_nodes = []            # interior node indices for which the distance is < 0
for i in range(0, np.size(idx)-1):
    print("* Shape %d/%d" %(i+1, np.size(idx)-1))
    d = dscaled[idx[i]:idx[i+1]]
    if_nodes = [(d[j][0], d[j][1]) for j in range(0, len(d))]
    # Add the first node to the end to create a closed loop for intersecting:
    if_loop = if_nodes
    if_loop.append(if_nodes[0])
    interface = LineString(if_loop) 
    
    # Intersect edges_close with the interface; construct interface nodes.
    for j in range(0, len(edges_close)):
        n1 = edges[edges_close[j], 0]
        n2 = edges[edges_close[j], 1]
        o = LineString( [p[n1,:], p[n2,:]] )
        x = interface.intersection(o)
        if (np.size(x) == 2):
            ifp.append( [x.x, x.y] )
        if (np.size(x) > 2):
            # print("\nWarning: mesh edge crosses multiple interface edges! "
            #     "Using the centroid for interface node position.")
            ifp.append( [x[:].centroid.x, x[:].centroid.y] )
        if ((j+1) % 1000 == 0 or (j+1) == len(edges_close)):
            print("\r- Processing edge %d / %d" %(j+1,len(edges_close)), end='')
    
    print("")
    
    # Determine if a node n is enclosed by the interface: Assuming r = (-1,-1) 
    # is 'outside', then if the number of interface intersections with edge
    # (r, n) is odd, n is enclosed.  
    for j in range(0, len(p)):
        o = LineString( [(-1,-1), (p[j][0],p[j][1])] )
        x = interface.intersection(o)
        ncross = np.size(x) // 2            # number of interface crossings
        if (ncross % 2 == 1):
            i_nodes.append(j)
        if ((j+1) % 1000 == 0 or (j+1) == len(p)):
            print("\r- Update node sign %d / %d" %(j+1,len(p)), end='')
    
    print("")


print("done.\n")


print("3/3 - computing final distances...")
d = cdist(p, ifp)
phi = np.min(d, axis=1)
phi[i_nodes] = -phi[i_nodes]        # negative distance for interior nodes
print("done.\n")


# Discard nodes above y threshold if requested. Requires retriangulation.
if (cropy != []):
    i = np.where(p[:,1] < cropy)
    p = p[i,:][0]
    phi = phi[i]
    t = Delaunay(p)
    t = t.simplices


# Write output
p = p.transpose()                   # 2xN for N nodes
t = t.transpose() + 1               # 3xN, +1 since Matlab indexing begins from 1
phi = np.array([phi]).transpose()   # Nx1
sio.savemat( 'nodes_phi.mat', {'p':p, 't':t, 'phi':phi} )
