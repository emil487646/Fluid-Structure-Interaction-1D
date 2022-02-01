import numpy as np
import matplotlib.pyplot as plt
import time


from matplotlib.animation import FuncAnimation
from matplotlib.animation import FFMpegWriter
from matplotlib.animation import PillowWriter

from piston import Piston
from piston import interpolate
from movingmesh import Mesh
from movingmesh import Fluid

t0=0
tol=10**-5
maxits=10
#Mesh size in fluid
M=20
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
tpf=[None]*maxits
timepiston=0
timefluid=0

while err>tol/5 and it<maxits:
    t0=time.time()
    piston.SDIRK12int(te, tol)
    timepiston+=time.time()-t0
    d_prev=tdva[1,-1]
    tdva=piston.output()
    tdva[1:,:]=-tdva[1:,:]
    tdva[1,:]+=1
    mesh=Mesh(0., tdva[1,0], M)
    fluid=Fluid(mesh, tdva)
    
    t0=time.time()
    fluid.EEint(te)
    timefluid+=time.time()-t0
    tp=fluid.output()
    piston=Piston(np.array([0., 0.]), tp)
    
    p_prev=p
    tpf[it]=tp
    p=tp[1,-1]
    if it>0:
        err=abs(p-p_prev)
        errvec[it-1]=err
        print(err)
    it+=1
    
#Plot error over time
# errtvec=[None]*maxits
# errvec=errvec[:it-1]
# errtvec=errtvec[:it-1]
# tpf=tpf[:it]
# t=tpf[it-1][0,:]
# for i in range(it-1):
#     tp1=tpf[i]
#     tp2=tpf[it-1]
#     errtvec[i]=np.array([interpolate(tp1[0], tp1[1], tp2[0, k])-tp2[1, k] for k in range(len(tp2[0, :]))])
#     #     errtvec[i]=np.array([abs(interpolate(tp1[0], tp1[1], tp2[0, s])-tp2[0, s]) for s in len(range(tp2[0, :]))])
# for errt in errtvec[:6:2]:
#     plt.plot(t, errt)
# plt.yscale("log")
# plt.xlabel(r'$t$')
# plt.ylabel(r'$e_j(t)$')
# plt.legend([r'$e_1(t)$', r'$e_3(t)$', r'$e_5(t)$'], loc='upper left')
# errtvec=errtvec[:it-1]



w=fluid.w
wtot=fluid.wtot
t=fluid.t
x=fluid.x
p=fluid.p
rho=w[:, :, 0]
u=w[:, :, 1]/rho
rhoE=w[:, :, 2]

# Plot error
I=range(len(errvec))
plt.plot(I, errvec)
plt.yscale("log")
plt.xlabel('j')
plt.ylabel(r'$e_j^N$')

#plot values at time t
# i=340
# ymin=-1
# ymax=5
# plt.plot(x[i,1:M+1], rho[i, :], 'r-')
# plt.plot(x[i,1:M+1], u[i, :], 'g-')
# plt.plot(x[i,1:M+1], p[i, :], 'b-')
# plt.plot([x[i,-1], x[i, -1]], [ymin, ymax], 'k-')
# tt=t[i]
# plt.text(0.6,2, f't={tt:.2f}', ha="center", va="top")
# plt.xlim(0, 1.2)
# plt.legend([r'$\rho$', r'$\rho u$', r'$\rho E$'], loc='upper left')
# plt.xlabel('t')

# Plot total values
# plt.plot(t, wtot[:,0], 'r-')
# plt.plot(t, wtot[:,1], 'g-')
# plt.plot(t, wtot[:,2], 'b-')
# plt.legend([r'$\rho$', r'$\rho u$', r'$\rho E$'], loc='upper left')
# plt.xlabel('t')

# Animate
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
# time = plt.text(0.6,2, str(0), ha="center", va="top")

# ymin=-1
# ymax=5
# def init():
#     ax.set_xlim(0, 1.2)
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