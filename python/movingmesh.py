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
        self.vl=xl*zp(0)
        self.vr=xr*zp(0)
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
        print(xf_new)
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
def pf(rhoE, rho, u):
    return (gamma-1)*(rhoE-0.5*rho*u**2)
def z(t):
    return 0.9+0.1*np.cos(10*np.pi*t)
def zp(t):
    return -np.pi*np.sin(10*np.pi*t)
# def z(t):
#     return 1
# def zp(t):
#     return 0
def f(w):
    rho=w[0]
    u=w[1]/rho
    rhou=w[1]
    rhoE=w[2]
    # p=(gamma-1)*rhoE-0.5*(gamma-1)*rho*np.abs(u)**2
    p=pf(rhoE, rho, u)
    return w*u+np.array([0, p, p*u])
def ft(w, v):
    rho=w[0]
    u=w[1]/rho
    rhou=w[1]
    rhoE=w[2]
    # p=(gamma-1)*rhoE-0.5*(gamma-1)*rho*np.abs(u)**2
    p=pf(rhoE, rho, u)
    return w*(u-v)+np.array([0, p, p*u])
def f_num(wl, wr, v, n):
    # print("flux")
    rhol=wl[0]
    ul=wl[1]/rhol
    rhoul=wl[1]
    rhoEl=wl[2]
    # pl=(gamma-1)*rhoEl-0.5*(gamma-1)*rhol*np.abs(ul)**2
    pl=pf(rhoEl, rhol, ul)
    cl=np.sqrt(gamma*pl/rhol)
    # print("pl=", pl)
    
    rhor=wr[0]
    ur=wr[1]/rhor
    rhour=wr[1]
    rhoEr=wr[2]
    # pr=(gamma-1)*rhoEr-0.5*(gamma-1)*rhor*np.abs(ur)**2
    pr=pf(rhoEr, rhor, ur)
    cr=np.sqrt(gamma*pr/rhor)
    # print("pr=", pr)
    
    A=np.abs(np.array([ul*n-v, ul*n+cl-v, ul*n-cl-v, ur*n-v, ur*n+cr-v, ur*n-cr-v]))
    # print(abs(A))
    lamdamax=max(A)
    # print("lmax=", lamdamax)
    fl=ft(wl, v)
    fr=ft(wr, v)
    return 0.5*(ft(wl, v)+ft(wr, v))*n-0.5*lamdamax*(wr-wl)*n, lamdamax

N=5
M=1000
# mesh=Mesh(0., 1., N)

def solve(N, M, T):
    mesh=Mesh(0., 1., N)
    dt=T/M
    tvec=np.linspace(0, T, M+1)
    u=np.zeros((M+1, N, 3))
    p=np.zeros((M+1, N))
    utot=np.zeros((M+1, 3))
    x=np.zeros((M+1, N))
    x[0,:]=mesh.midpoints()
    u[0,:,0]=1
    u[0,:,1]=0
    u[0,:,2]=2.5
    utot[0]=mesh.cells[0].volume*np.sum(u[0], axis=0)
    p[0]=pf(u[0,:,2], u[0,:,0], u[0,:,0]*u[0,:,1])
    for m in range(1, M+1):
        print("m=", m)
        t=tvec[m]
        upd=np.zeros((N, 3))
        
        
        #Initilize the new values to the values at the previous time step
        # u[m, :, :]=[u[m-1, k, :] for k in range(N)]
        
        
        cell0=mesh.cells[0]
        # u[m, 0,:]+=dt*f_num(np.array([u[m-1, 0, 0], -u[m-1, 0, 1], u[m-1, 0, 2]]), u[m-1, 0], cell0.vl,1)/cell0.volume
        # u[m, 0,:]+=dt*f(np.array([u[m-1, 0, 0], 0, u[m-1, 0, 2]]))/cell0.volume
        # flux=np.array([0, p[m-1, 0], 0])
        flux=np.array([0, (gamma-1)*u[m-1, 0, 2], 0])
        upd[0,:]+=flux
        # print((gamma-1)*u[m-1, 0, 2])
        # print(flux)
        #Loop over the edges
        for n in range(N-1):
            celll=mesh.cells[n]
            cellr=mesh.cells[n+1]
            #Approximate the flux by the numerical flux times the edge length
            #Note that celll.vr=cellr.vl
            # print("n=", n)
            flux, lmax=f_num(u[m-1, n], u[m-1, n+1], celll.vr, 1)*1
            #Update values
            upd[n,:]-=flux
            # print(u[m-1, n+1])
            print(flux)
            upd[n+1,:]+=flux
        cellf=mesh.cells[N-1]
        # u[m, N-1,:]-=dt*f_num(u[m-1, N-1], u[m-1, N-1], cellf.vr,1)/cellf.volume
        # u[m, N-1,:]-=dt*f(u[m-1, N-1])/cellf.volume
        # u[m, N-1,:]-=dt*f(np.array([u[m-1, N-1, 0], mesh.cells[N-1].vr, u[m-1, N-1, 2]]))/cellf.volume
        
        # u[m, 0, 1]=0
        # u[m, N-1, 1]=0
        # u[m, N-1,:]-=dt*f_num(u[m-1, N-1]-np.array([0, cellf.vr, 0]), np.array([u[m-1, N-1, 0], -u[m-1, N-1, 1], u[m-1, N-1, 2]]), cellf.vr, 1)/cellf.volume
        # u[m, N-1, 1]=mesh.cells[N-1].vr
        # u[m, N-1,:]-=dt*np.array([0, p[m-1, 0], p[m-1, 0]*cellf.vr])/cellf.volume
        # u[m, N-1,:]-=dt*f(np.array([u[m-1, N-1, 0], u[m-1, N-1, 0]*cellf.vr, u[m-1, N-1, 2]]))/cellf.volume
        # u[m, N-1,:]-=dt*ft(np.array([u[m-1, N-1, 0], cellf.vr*u[m-1, N-1, 0], u[m-1, N-1, 2]]), cellf.vr)/cellf.volume
        
        
        # flux=ft(np.array([u[m-1, N-1, 0], cellf.vr*u[m-1, N-1, 0], u[m-1, N-1, 2]]), cellf.vr)
        # # flux=np.array([0, (gamma-1)*u[m-1, N-1, 2], 0])
        # print(flux)
        p_t=(gamma-1)*u[m-1, N-1, 2]-0.5*(gamma-1)*u[m-1, N-1, 0]*np.abs(cellf.vr)**2
        p_t=pf(u[m-1, N-1, 2], u[m-1, N-1, 0], cellf.vr)
        # p_t=1
        print("pressure=", p_t)
        # u_t=np.array([u[m-1, N-1, 0], cellf.vr*u[m-1, N-1, 0], u[m-1, N-1, 2]])
        # flux=cellf.vr*u_t+np.array([0, p_t, p_t*cellf.vr])
        # print(flux)
        
        flux=np.array([0, p_t, p_t*cellf.vr])
        # flux=np.array([0, (gamma-1)*u[m-1, N-1, 2], 0])
        # flux[1]=15.85
        print(flux)
        upd[N-1,:]-=flux
        
        #Initilize the new values to the values at the previous time step
        for k in range(N):
            u[m, k, :]=u[m-1, k, :]*mesh.cells[k].volume
        #Calculate the new cell positions, velocities and volumes.
        
        mesh.update(z(t), dt)
        x[m,:]=mesh.midpoints()
        
        u[m, :, :]+=dt*upd
        for k in range(N):
            u[m, k, :]/=mesh.cells[k].volume
        
        print("Target velocity: ", zp(t))
        print("Target velocity: ", cellf.vr)
        print("Computed velocity: ", u[m, N-1, 1]/u[m, N-1, 0])
        
        
        #Post computations
        utot[m]=mesh.cells[0].volume*np.sum(u[m], axis=0)
        
        rho=u[m, :, 0]
        u1=u[m, :, 1]/rho
        v=np.array([(mesh.cells[k].vl+mesh.cells[k].vr)/2 for k in range(N)])
        rhoE=u[m, :, 2]
        p[m]=(gamma-1)*rhoE-0.5*(gamma-1)*rho*np.abs(u1)**2
    return tvec, x, u, utot, p
        

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

# N=20
# M=20000
N=5
M=2
T=0.002
t, x, w, wtot, p=solve(N, M, T)
dt=T/M
# fig, ax = plt.subplots()
# xdata, ydata = [], []
# ln0, = plt.plot([], [], 'r-')
# ln1, = plt.plot([], [], 'g-')
# ln2, = plt.plot([], [], 'b-')
# ln3, = plt.plot([], [], 'k-')
rho=w[:, :, 0]
u=w[:, :, 1]/rho
rhoE=w[:, :, 2]
# p=(gamma-1)*rhoE-0.5*rho*np.abs(u)**2

plt.plot(t, wtot[:])
plt.legend(['rho', 'u', 'p'], loc='upper left')
#Frames per second in video
fps=60
#Time units per second in video
tps=3
# time = plt.text(0.5,1.8, str(0), ha="left", va="top")

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
    
    ln0.set_data(x[i,:], rho[i, :])
    ln1.set_data(x[i,:], u[i, :])
    ln2.set_data(x[i,:], p[i, :])
    ln3.set_data([x[i,-1], x[i, -1]], [-0.5, 2])
    # plt.title('t='+str(t[frame]))
    tt=t[i]
    time.set_text(f't={tt:.2f}')
    return ln0, ln1, ln2, ln3, time

# ani = FuncAnimation(fig, update, frames=np.array(range(M)),
#                     init_func=init, blit=True, repeat=True)
# plt.show()
# f = "animation4.mp4" 
# writergif = FFMpegWriter(fps=int(M/T/12), bitrate=-1)
# ani.save(f, writer=writergif)

# ani = FuncAnimation(fig, update, frames=np.array(tps*fps*T),
#                     init_func=init, blit=True, repeat=True)
    
# plt.show()
# f = "animation4.mp4" 
# writergif = FFMpegWriter(fps=fps, bitrate=-1)
# ani.save(f, writer=writergif)