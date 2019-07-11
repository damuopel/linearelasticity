import numpy as np
from numpy.linalg import inv,det
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
    if dFlag==1:
        dN1dXi = -0.25*(1-Eta)
        dN1dEta = -0.25*(1-Xi)
        dN2dXi = 0.25*(1-Eta)
        dN2dEta = -0.25*(1+Xi)
        dN3dXi = 0.25*(1+Eta)
        dN3dEta = 0.25*(1+Xi)
        dN4dXi = -0.25*(1+Eta)
        dN4dEta = 0.25*(1-Xi)
        N = np.array([[dN1dXi,dN2dXi,dN3dXi,dN4dXi],[dN1dEta,dN2dEta,dN3dEta,dN4dEta]])
    else:
        N1 = 0.25*(1-Xi)*(1-Eta)
        N2 = 0.25*(1+Xi)*(1-Eta)
        N3 = 0.25*(1+Xi)*(1+Eta)
        N4 = 0.25*(1-Xi)*(1+Eta)
        N = np.array([N1,N2,N3,N4])
            
    return N

def Integration(Topology,XY):
    xyGPI = np.array([[-0.5774,-0.5774,0.5774,0.5774],[-0.5774,0.5774,-0.5774,0.5774]])
    hGPI = np.array([1,1,1,1])
    refNodes = Topology[:,0]-1 # Pick a reference element
    refVerts = XY[:,refNodes]
    Ke = 0
    for iGP in range(4):
        Xi = xyGPI[0,iGP]
        Eta = xyGPI[1,iGP]
        H = hGPI[iGP]
        dNl = ShapeFunctions(Xi,Eta,1)
        Jacobian = refVerts@dNl.T
        dNg = inv(Jacobian)@dNl
        B = np.array([[dNg[0,0],0,dNg[0,1],0,dNg[0,2],0,dNg[0,3],0],[0,dNg[1,0],0,dNg[1,1],0,dNg[1,2],0,dNg[1,3]],[dNg[1,0],dNg[0,0],dNg[1,1],dNg[0,1],dNg[1,2],dNg[0,2],dNg[1,3],dNg[0,3]]])
        Ke = Ke + B.T@D@B*det(Jacobian)*H     

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