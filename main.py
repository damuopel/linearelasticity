import numpy as np
import math

def Mesh(eSize,xElms,yElms):
    # Coordinates
    start = 0
    xStop = eSize*xElms
    yStop = eSize*yElms
    x = np.linspace(start,xStop,xElms+1)
    y = np.linspace(start,yStop,yElms+1)
    X,Y = np.meshgrid(x,y)
    XY = np.array([X.flatten(),Y.flatten()])
    # Topology
    node1 = np.array([iElm+math.ceil(iElm/xElms)-1 for iElm in list(range(1,xElms*yElms+1))])
    node2 = node1 + 1
    node3 = node2 + xElms + 1
    node4 = node1 + xElms + 1
    Topology = np.array([node1,node2,node3,node4])
    
    return Topology,XY

def ShapeFunctions(Xi,Eta,dFlag):
	if dFlag:
		
	else:
		N1 = 0.25*(1-Xi)*(1-Eta)
		N2 = 0.25*(1+Xi)*(1-Eta)
		N3 = 0.25*(1+Xi)*(1+Eta)
		N4 = 0.25*(1-Xi)*(1+Eta)
		N = np.array(N1,N2,N3,N4)

    return N

def Integration(Topology,XY):
    elms = Topology.size/Topology.ndim
    xyGPI = np.array([-0.5774,-0.5774,0.5774,0.5774],[-0.5774,0.5774,-0.5774,0.5774])
    hGPI = np.array([1,1,1,1])
    # Same stiffness matrix per element
        
    return K,F

def Solver(K,F):
    
    return u

if __name__ == '__main__':
    # Input User
    eSize = 1
    xElms = 3
    yElms = 2
	# Create Mesh
    Topology,XY = Mesh(eSize,xElms,yElms)
    K,F = Integration(Topology,XY)
    u = Solver(K,F)