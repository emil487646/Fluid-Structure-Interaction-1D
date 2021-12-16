import numpy as np
import matplotlib.pyplot as plt

from matplotlib.animation import FuncAnimation
from matplotlib.animation import FFMpegWriter
from matplotlib.animation import PillowWriter
# from matplotlib import animation

class Cell:
    def __init__(self, xl, xr):
        self.xl=xl
        self.xr=xr
        self.mid=(xl+xr)/2
        self.vl=0
        self.vr=0
        self.volume=xr-xl

class Mesh:
    def __init__(self, x0, xf, N):
        dx=(xf-x0)/N
        self.cells0=[Cell(x0+n*dx, x0+(n+1)*dx) for n in range(N)]
        self.cells=[Cell(x0+n*dx, x0+(n+1)*dx) for n in range(N)]
        self.N=N
        self.x0=x0
        self.xf=xf
        self.z0=xf
    def update(self, xf_new, dt):
        xf_old=self.xf
        # xf_old=self.cells0[0].xl
        x0=self.x0
        factor=(xf_new-x0)/(xf_old-x0)
        for n in range(self.N):
            xl_old=self.cells[n].xl
            xr_old=self.cells[n].xr
            xl_new=(x0+(self.cells0[n].xl-x0)*factor)
            xr_new=(x0+(self.cells0[n].xr-x0)*factor)
            self.cells[n].xl=xl_new
            self.cells[n].xr=xr_new
            self.cells[n].mid=(self.cells[n].xl+self.cells[n].xr)/2
            self.cells[n].volume=(self.cells[n].xr-self.cells[n].xl)
            diffl=xl_new-xl_old
            diffr=xr_new-xr_old
            self.cells[n].vl=diffl/dt
            self.cells[n].vr=diffr/dt
    def midpoints(self):
        return np.array([self.cells[n].mid for n in range(N)])
gamma=7./5.
c=1
def z(t):
    return 1+0.1*np.sin(10*np.pi*t)
def f(w):
    rho=w[0]
    u=w[1]/rho
    rhoE=w[2]
    p=(gamma-1)*rhoE-0.5*rho*np.abs(u)**2
    return w*u+np.array([0, p, p*u])
def ft(w, v):
    rho=w[0]
    u=w[1]/rho
    rhoE=w[2]
    p=(gamma-1)*rhoE-0.5*rho*np.abs(u)**2
    return w*(u-v)+np.array([0, p, p*u])
def f_num(wl, wr, v, n):
    rhol=wl[0]
    ul=wl[1]/rhol
    rhoEl=wl[2]
    pl=(gamma-1)*rhoEl-0.5*rhol*np.abs(ul)**2
    cl=np.sqrt(gamma*pl/rhol)
    
    rhor=wr[0]
    ur=wr[1]/rhor
    rhoEr=wr[2]
    pr=(gamma-1)*rhoEr-0.5*rhor*np.abs(ur)**2
    cr=np.sqrt(gamma*pr/rhor)
    
    A=np.abs(np.array([ul*n-v, ul*n+cl-v, ul*n-cl-v, ur*n-v, ur*n+cr-v, ur*n-cr-v]))
    # print(abs(A))
    lamdamax=max(A)
    # print("lmax=", lamdamax)
    fl=ft(wl, v)
    fr=ft(wr, v)
    return 0.5*(ft(wl, v)+ft(wr, v))*n-0.5*lamdamax*(wr-wl)*n

N=5
M=1000
# mesh=Mesh(0., 1., N)

def solve(N, M, T):
    mesh=Mesh(0., 1., N)
    dt=T/M
    tvec=np.linspace(0, T, M+1)
    u=np.zeros((M+1, N, 3))
    x=np.zeros((M+1, N))
    x[0,:]=mesh.midpoints()
    u[0,:,0]=1
    u[0,:,1]=0
    u[0,:,2]=2.5
    for m in range(1, M+1):
        print("m=", m)
        t=tvec[m]
        #Initilize the new values to the values at the previous time step
        for k in range(N):
            u[m, k, :]=u[m-1, k, :]*mesh.cells[k].volume
        #Calculate the new cell positions, velocities and volumes.
        mesh.update(z(t), dt)
        x[m,:]=mesh.midpoints()
        
        for k in range(N):
            u[m, k, :]/=mesh.cells[k].volume
        
        #Initilize the new values to the values at the previous time step
        # u[m, :, :]=[u[m-1, k, :] for k in range(N)]
        
        #Loop over the edges
        for n in range(N-1):
            celll=mesh.cells[n]
            cellr=mesh.cells[n+1]
            #Approximate the flux by the numeircal flux times the edge length
            #Note that celll.vr=cellr.vl
            flux=f_num(u[m-1, n], u[m-1, n+1], celll.vr, 1)*1
            #Update values
            u[m,n,:]-=dt*flux/celll.volume
            u[m,n+1,:]+=dt*flux/cellr.volume
        u[m, 0, 1]=0
        # u[m, N-1, 1]=0
        u[m, N-1, 1]=mesh.cells[N-1].vr
    return tvec, x, u
        

# t=0
# dt=0.1
# T=5
# M=100
# dt=T/M
# fig, ax = plt.subplots()
# xdata, ydata = [], []
# ln, = plt.plot([], [], 'ro')

# def init():
#     ax.set_xlim(0, 2.)
#     ax.set_ylim(0, 2)
#     return ln,

# def update(frame):
#     # plt.plot(mesh.midpoints(), np.ones(N), '*')
#     ln.set_data(mesh.midpoints(), np.ones(N))
#     t=dt*frame
#     print(z(t))
#     mesh.update(z(t), dt)
#     return ln,

# ani = FuncAnimation(fig, update, frames=np.linspace(0, 5, M),
#                     init_func=init, blit=True)
# plt.show()

N=50
M=10000
# N=10
# M=400
T=5.
t, x, w=solve(N, M, T)
dt=T/M
fig, ax = plt.subplots()
xdata, ydata = [], []
ln0, = plt.plot([], [], 'r-')
ln1, = plt.plot([], [], 'g-')
ln2, = plt.plot([], [], 'b-')
ln3, = plt.plot([], [], 'k-')
rho=w[:, :, 0]
u=w[:, :, 1]/rho
rhoE=w[:, :, 2]
p=(gamma-1)*rhoE-0.5*rho*np.abs(u)**2
plt.legend(['rho', 'u', 'p'], loc='upper left')
time = plt.text(0.5,1.8, str(0), ha="left", va="top")

def init():
    ax.set_xlim(0, 1.1)
    ax.set_ylim(-0.5, 2)
    return ln0, ln1, ln2, ln3,

def update(frame):
    # plt.plot(mesh.midpoints(), np.ones(N), '*')
    # print(frame)
    ln0.set_data(x[frame,:], rho[frame, :])
    ln1.set_data(x[frame,:], u[frame, :])
    ln2.set_data(x[frame,:], p[frame, :])
    ln3.set_data([x[frame,-1], x[frame, -1]], [-0.5, 2])
    # plt.title('t='+str(t[frame]))
    tt=t[frame]
    time.set_text(f't={tt:.2f}')
    return ln0, ln1, ln2, ln3, time

ani = FuncAnimation(fig, update, frames=np.array(range(M)),
                    init_func=init, blit=True, repeat=True)
plt.show()
f = "animation4.mp4" 
writergif = FFMpegWriter(fps=int(M/T/12), bitrate=-1)
ani.save(f, writer=writergif)