import numpy as np
import matplotlib.pyplot as plt

from piston import SDIRK12step
from movingmesh import Mesh
from movingmesh import Fluid
from movingmesh import step
def couple(tol=10**-4, maxsteps=10000, maxits=10):
    t0=0
    M=maxsteps
    #Mesh size in lfuid
    N=10
    #Final time
    tf=1.
    #Mesh of fluid
    mesh=Mesh(0., 1., N)
    fluid=Fluid(mesh)
    te=tf
    err=2*tol
    it=0
    #Initial guess
    fluid.input(np.array([[t0,te],[0, 0],[0, 0],[0, 00]]))
    while err>tol and it<maxits:
        fluid.solve(te)
        p=fluid.output()
        piston.input(p)
        tdva=piston.output()
        fluid.input(tdva)
        it+=1
    
    return 0;

