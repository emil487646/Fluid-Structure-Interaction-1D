import numpy as np
import matplotlib.pyplot as plt
import time


from matplotlib.animation import FuncAnimation
from matplotlib.animation import FFMpegWriter
from matplotlib.animation import PillowWriter

from piston import Piston
from movingmesh import Mesh
from movingmesh import Fluid
from movingmesh import step

t0=0
tol=10**-5
maxits=10
#Mesh size in lfuid
M=100
#Final time
T=5.
#Mesh of fluid
mesh=Mesh(0., 1., M)
tdva=np.array([[0,T],[1, 1],[0, 0],[0, 0]])
fluid=Fluid(mesh, tdva)

tp=np.array([[0,T],[1, 1]])
piston=Piston(np.array([0., 0.]), tp)

te=T
err=2*tol
it=0
#Initial guess
# fluid.input(np.array([[t0,te],[0, 0],[0, 0],[0, 00]]))
d_prev=0
p_prev=p=0
errvec=np.zeros(maxits)
timepiston=0
timefluid=0

while err>tol/5 and it<maxits:
    t0=time.time()
    piston.SDIRK12int(te, tol)
    timepiston+=time.time()-t0
    d_prev=tdva[1,-1]
    tdva=piston.output()
    # tdva=-tdva
    tdva[1,:]+=1
    mesh=Mesh(0., tdva[1,0], M)
    fluid=Fluid(mesh, tdva)
    
    t0=time.time()
    fluid.EEint(te)
    timefluid+=time.time()-t0
    tp=fluid.output()
    piston=Piston(np.array([0., 0.]), tp)
    
    p_prev=p
    p=tp[1,-1]
    err=abs(p-p_prev)
    errvec[it]=err
    print(err)
    it+=1
errvec=errvec[:it]

w=fluid.w
wtot=fluid.wtot
t=fluid.t
x=fluid.x
p=fluid.p
rho=w[:, :, 0]
u=w[:, :, 1]/rho
rhoE=w[:, :, 2]

# p=(gamma-1)*rhoE-0.5*rho*np.abs(u)**2

#plot values at time t
i=427
ymin=-4
ymax=20
plt.plot(x[i,1:M+1], rho[i, :], 'r-')
plt.plot(x[i,1:M+1], u[i, :], 'g-')
plt.plot(x[i,1:M+1], p[i, :], 'b-')
plt.plot([x[i,-1], x[i, -1]], [ymin, ymax], 'k-')
tt=t[i]
plt.text(1.5,15, f't={tt:.2f}', ha="center", va="top")
plt.xlim(0, 3.2)
plt.legend([r'$\rho$', r'$\rho u$', r'$\rho E$'], loc='upper left')
plt.xlabel('t')

#Plot total values
# plt.plot(t, wtot[:,0], 'r-')
# plt.plot(t, wtot[:,1], 'g-')
# plt.plot(t, wtot[:,2], 'b-')
# plt.legend([r'$\rho$', r'$\rho u$', r'$\rho E$'], loc='upper left')
# plt.xlabel('t')

# fig, ax = plt.subplots()
# ln0, = plt.plot([], [], 'r-')
# ln1, = plt.plot([], [], 'g-')
# ln2, = plt.plot([], [], 'b-')
# ln3, = plt.plot([], [], 'k-')

# plt.legend([r'$\rho$', r'$\rho u$', r'$\rho E$'], loc='upper left')
# #Frames per second in video
# fps=60
# #Time units per second in video
# tps=12
# time = plt.text(1.5,15, str(0), ha="left", va="top")

# ymin=-4
# ymax=20
# def init():
#     ax.set_xlim(0, 3.2)
#     ax.set_ylim(ymin, ymax)
#     return ln0, ln1, ln2, ln3, 

# def update(frame):
#     # plt.plot(mesh.midpoints(), np.ones(N), '*')
#     # print(frame)
    
#     # i=frame
#     i=0
#     while i<len(t)-1 and t[i]<frame/fps/tps:
#         i=i+1
    
#     ln0.set_data(x[i,1:M+1], rho[i, :])
#     ln1.set_data(x[i,1:M+1], u[i, :])
#     ln2.set_data(x[i,1:M+1], p[i, :])
#     ln3.set_data([x[i,-1], x[i, -1]], [ymin, ymax])
#     # plt.title('t='+str(t[frame]))
#     tt=t[i]
#     time.set_text(f't={tt:.2f}')
#     return ln0, ln1, ln2, ln3, time

        
# ani = FuncAnimation(fig, update, frames=np.array(int(tps*fps*T)), init_func=init, blit=True, repeat=True)
    
# plt.show()
# f = "animation4.mp4" 
# writergif = FFMpegWriter(fps=fps, bitrate=-1)
# ani.save(f, writer=writergif)