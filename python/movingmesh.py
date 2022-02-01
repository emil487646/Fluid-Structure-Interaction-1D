import numpy as np
import matplotlib.pyplot as plt
import bisect

from matplotlib.animation import FuncAnimation
from matplotlib.animation import FFMpegWriter
from matplotlib.animation import PillowWriter
# from matplotlib import animation
def interpolate(t, f, s):
    index=max(bisect.bisect_left(t, s)-1,0)
    h=s-t[index]
    return f[index]+(f[index+1]-f[index])/(t[index+1]-t[index])*(s-t[index])
class Cell:
    def __init__(self, xl, xr):
        self.xl=xl
        self.xr=xr
        self.mid=(xl+xr)/2
        self.vl=0#xl*zp(0)
        self.vr=0#xr*zp(0)
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
    def update2(self, xf_new, vf_new):
        xf_old=self.xf
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
            
            vl_new=(vf_new*xl_new/xf_new)
            vr_new=(vf_new*xr_new/xf_new)
            self.cells[n].vl=vl_new
            self.cells[n].vr=vr_new
    def midpoints(self):
        return np.array([self.cells[n].mid for n in range(self.N)])
class Fluid:
    def __init__(self, mesh, tdva, maxsteps=10000):
        self.gamma=7./5.
        
        self.N=maxsteps
        N=self.N
        self.M=mesh.N
        M=self.M
        
        self.mesh=Mesh(0., 1., M)
        self.x=np.zeros((N+1, M+2))
        self.x[0,1:M+1]=self.mesh.midpoints()
        self.x[:,0]=0
        self.x[0,M+1]=1
        
        self.dt=np.zeros(N)
        
        t0=0
        self.t=np.zeros(N+1)
        self.dt=np.zeros(N)
        self.t[0]=t0
        
        self.w=np.zeros((N+1, M, 3))
        self.w[0,:,0]=1
        self.w[0,:,1]=0#np.sin(np.pi*self.x[0])
        self.w[0,:,2]=2.5
        
        
        self.wtot=np.zeros((N+1, 3))
        self.wtot[0]=self.mesh.cells[0].volume*np.sum(self.w[0], axis=0)
        
        self.pr=np.zeros(N+1)
        self.p=np.zeros((N+1, M))
        self.p[0]=pf(self.w[0,:,2], self.w[0,:,0], self.w[0,:,0]*self.w[0,:,1])
        self.pr[0]=pf(self.w[0, M-1, 2], self.w[0, M-1, 0], tdva[2,0])
        
        #Boundary data (t, d, v, a)
        self.tdva=tdva
        
        #Current step
        self.n=0
    def interpolateinput(self, s):
        tdva=self.tdva
        d=interpolate(tdva[0], tdva[1], s)
        v=interpolate(tdva[0], tdva[2], s)
        a=interpolate(tdva[0], tdva[3], s)
        return np.array([d, v, a])
    def step(self, dva_n):
        t_n=self.t[self.n]
        w_n=self.w[self.n]
        mesh=self.mesh
        M=mesh.N
        
        upd=np.zeros((M, 3))
        lmaxmax=0
        
        #Left boundary
        cell0=mesh.cells[0]
        flux=np.array([0, (gamma-1)*w_n[0, 2], 0])
        upd[0,:]+=flux
        
        #Internal cells
        # upd[0:-1]-=np.array([f_num(u[m-1, n], u[m-1, n+1], mesh.cells[n].vr, 1)[0] for n in range(N-1)])
        # upd[1:]+=np.array([f_num(u[m-1, n], u[m-1, n+1], mesh.cells[n].vr, 1)[0] for n in range(N-1)])
        for m in range(M-1):
            celll=mesh.cells[m]
            cellr=mesh.cells[m+1]
            #Approximate the flux by the numerical flux times the edge length
            #Note that celll.vr=cellr.vl
            flux, lmax=f_num(w_n[m], w_n[m+1], celll.vr, 1)*1
            lmaxmax=max(lmaxmax, lmax)
            #Update values
            upd[m,:]-=flux
            upd[m+1,:]+=flux
            
        cellf=mesh.cells[M-1]
        
        #Interpolate the input variables
        d=dva_n[0]
        v=dva_n[1]
        a=dva_n[2]
        
        p=pf(w_n[M-1, 2], w_n[M-1, 0], v)
        if p<0:
            print("Negative pressure, p=", p)
        
        rho=w_n[M-1,0]
        c=np.sqrt(gamma*p/rho)
        lmaxmax=max(lmaxmax, c)
        
        #Right boundary
        forcedensity=w_n[M-1, 0]*cellf.volume*a
        F=np.array([0, forcedensity, v*forcedensity])
        flux=np.array([0, p, p*v])
        upd[M-1,:]-=flux
        upd[M-1,:]+=F
        
        return upd, lmaxmax
    def EEstep(self, dt):
        return 0
    #Integrate until te. Time step is chosen adaptively to ensure stability
    def EEint(self, te):
        n=0
        t0=self.t[self.n]
        M=self.M
        N=self.N

        while n<N and self.t[n]<te-10**(-15):
            dva_n=self.interpolateinput(self.t[n])
            upd, lmax=self.step(dva_n)
            
            cellf=self.mesh.cells[M-1]
            dx=cellf.volume
            self.dt[n]=dx/lmax
            
            self.t[n+1]=self.t[n]+self.dt[n]
            if self.t[n+1]>te:
                self.dt=self.dt[:n+1]
                self.dt[n]=te-self.t[n]
                self.t=self.t[:n+2]
                self.t[n+1]=te
                self.w=self.w[:n+2,:,:]
                self.wtot=self.wtot[:n+2,:]
                self.p=self.p[:n+2,:]
                self.pr=self.pr[:n+2]
                
            #Initilize the new values to the values at the previous time step
            # for m in range(M):
            #     self.w[n+1, m, :]=self.w[n, m, :]*self.mesh.cells[m].volume
            self.w[n+1, :, :]=[self.w[n, m, :]*self.mesh.cells[m].volume for m in range(M)]
                
            #Calculate the new cell positions, velocities and volumes.
            dva_np1=self.interpolateinput(self.t[n+1])
            self.mesh.update(dva_np1[0], self.dt[n])
            self.x[n+1,1:M+1]=self.mesh.midpoints()
            self.x[n+1,M+1]=self.mesh.cells[M-1].xr
            
            self.w[n+1, :, :]+=self.dt[n]*upd
            # for m in range(M):
            #     self.w[n+1, m, :]/=self.mesh.cells[m].volume
            self.w[n+1, :, :]=[self.w[n+1, m, :]/self.mesh.cells[m].volume for m in range(M)]
            
            
            #Post computations
            self.wtot[n+1]=self.mesh.cells[0].volume*np.sum(self.w[n+1], axis=0)
            
            rho=self.w[n+1, :, 0]
            u1=self.w[n+1, :, 1]/rho
            rhoE=self.w[n+1, :, 2]
            self.p[n+1]=(gamma-1)*rhoE-0.5*(gamma-1)*rho*np.abs(u1)**2
            self.pr[n+1]=pf(self.w[n+1, M-1, 2], self.w[n+1, M-1, 0], dva_np1[1])
            
            n=n+1
            self.n=n
    def output(self):
        return np.stack((self.t, self.pr))
gamma=7./5.
def pf(rhoE, rho, u):
    return (gamma-1)*(rhoE-0.5*rho*u**2)
def f(w):
    rho=w[0]
    u=w[1]/rho
    rhou=w[1]
    rhoE=w[2]
    p=pf(rhoE, rho, u)
    return w*u+np.array([0, p, p*u])
def ft(w, v):
    rho=w[0]
    u=w[1]/rho
    rhou=w[1]
    rhoE=w[2]
    p=pf(rhoE, rho, u)
    return w*(u-v)+np.array([0, p, p*u])
def f_num(wl, wr, v, n):
    rhol=wl[0]
    ul=wl[1]/rhol
    rhoul=wl[1]
    rhoEl=wl[2]
    pl=pf(rhoEl, rhol, ul)
    cl=np.sqrt(gamma*pl/rhol)
    
    rhor=wr[0]
    ur=wr[1]/rhor
    rhour=wr[1]
    rhoEr=wr[2]
    pr=pf(rhoEr, rhor, ur)
    cr=np.sqrt(gamma*pr/rhor)
    
    A=np.abs(np.array([ul*n-v, ul*n+cl-v, ul*n-cl-v, ur*n-v, ur*n+cr-v, ur*n-cr-v]))
    lamdamax=max(A)
    return 0.5*(ft(wl, v)+ft(wr, v))*n-0.5*lamdamax*(wr-wl)*n, lamdamax
