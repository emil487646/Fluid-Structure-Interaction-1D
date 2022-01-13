import numpy as np
import matplotlib.pyplot as plt


from matplotlib.animation import FuncAnimation
from matplotlib.animation import FFMpegWriter
from matplotlib.animation import PillowWriter

from piston import Piston
from movingmesh import Mesh
from movingmesh import Fluid
from movingmesh import step
def couple(tol=10**-4, maxsteps=10000, maxits=1):
    t0=0
    M=maxsteps
    #Mesh size in lfuid
    M=20
    #Final time
    T=5.
    #Mesh of fluid
    mesh=Mesh(0., 1., M)
    tdva=np.array([[0,T],[1, 1],[0, 0],[0, 0]])
    fluid=Fluid(mesh, tdva)
    te=T
    err=2*tol
    it=0
    #Initial guess
    # fluid.input(np.array([[t0,te],[0, 0],[0, 0],[0, 00]]))
    while err>tol and it<maxits:
        fluid.EEint(te)
        tp=fluid.output()
        piston=Piston(np.array([0., 0.]), tp)
        piston.SDIRK12int(te, 10**-2)
        tdva=piston.output()
        if err>tol and it+1<maxits:
            fluid=Fluid(mesh, tdva)
        it+=1

    w=fluid.w
    t=fluid.t
    x=fluid.x
    p=fluid.p
    
    fig, ax = plt.subplots()
    xdata, ydata = [], []
    ln0, = plt.plot([], [], 'r-')
    ln1, = plt.plot([], [], 'g-')
    ln2, = plt.plot([], [], 'b-')
    ln3, = plt.plot([], [], 'k-')
    rho=w[:, :, 0]
    u=w[:, :, 1]/rho
    rhoE=w[:, :, 2]
    # p=(gamma-1)*rhoE-0.5*rho*np.abs(u)**2
    
    # plt.plot(t, wtot[:])
    plt.legend(['rho', 'u', 'p'], loc='upper left')
    #Frames per second in video
    fps=60
    #Time units per second in video
    tps=12
    time = plt.text(0.5,1.8, str(0), ha="left", va="top")
    
    def init():
        ax.set_xlim(0, 1.1)
        ax.set_ylim(-0.5, 2)
        return ln0, ln1, ln2, ln3,
    
    def update(frame):
        # plt.plot(mesh.midpoints(), np.ones(N), '*')
        # print(frame)
        
        # i=frame
        i=0
        while i<len(t)-1 and t[i]<frame/fps/tps:
            i=i+1
        
        ln0.set_data(x[i,1:M+1], rho[i, :])
        ln1.set_data(x[i,1:M+1], u[i, :])
        ln2.set_data(x[i,1:M+1], p[i, :])
        ln3.set_data([x[i,-1], x[i, -1]], [-0.5, 2])
        # plt.title('t='+str(t[frame]))
        tt=t[i]
        time.set_text(f't={tt:.2f}')
        return ln0, ln1, ln2, ln3, time
    
            
    ani = FuncAnimation(fig, update, frames=np.array(int(tps*fps*T)),
                        init_func=init, blit=True, repeat=True)
        
    plt.show()
    # f = "animation4.mp4" 
    # writergif = FFMpegWriter(fps=fps, bitrate=-1)
    # ani.save(f, writer=writergif)
    
couple()
# M=20
# mesh=Mesh(0., 1., M)
# T=5
# #Initial guess
# tdva=np.array([[0,T],[1, 1],[0, 0],[0, 0]])
# fluid=Fluid(mesh, tdva)
# fluid.EEint(T)

# w=fluid.w
# t=fluid.t
# x=fluid.x
# p=fluid.p

# fig, ax = plt.subplots()
# xdata, ydata = [], []
# ln0, = plt.plot([], [], 'r-')
# ln1, = plt.plot([], [], 'g-')
# ln2, = plt.plot([], [], 'b-')
# ln3, = plt.plot([], [], 'k-')
# rho=w[:, :, 0]
# u=w[:, :, 1]/rho
# rhoE=w[:, :, 2]
# # p=(gamma-1)*rhoE-0.5*rho*np.abs(u)**2

# # plt.plot(t, wtot[:])
# plt.legend(['rho', 'u', 'p'], loc='upper left')
# #Frames per second in video
# fps=60
# #Time units per second in video
# tps=12
# time = plt.text(0.5,1.8, str(0), ha="left", va="top")

# def init():
#     ax.set_xlim(0, 1.1)
#     ax.set_ylim(-0.5, 2)
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
#     ln3.set_data([x[i,-1], x[i, -1]], [-0.5, 2])
#     # plt.title('t='+str(t[frame]))
#     tt=t[i]
#     time.set_text(f't={tt:.2f}')
#     return ln0, ln1, ln2, ln3, time

        
# ani = FuncAnimation(fig, update, frames=np.array(tps*fps*T),
#                     init_func=init, blit=True, repeat=True)
    
# plt.show()
# # f = "animation4.mp4" 
# # writergif = FFMpegWriter(fps=fps, bitrate=-1)
# # ani.save(f, writer=writergif)