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
        # print(xf_new)
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
    def update2(self, xf_new, vf_new):
        # print(xf_new)
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
        # self.tdva=np.zeros((N+1, 4))
        self.tdva=tdva
        
        #Current step
        self.n=0
    def interpolateinput(self, s):
        # d=z(s)
        # v=zp(s)
        # a=zpp(s)
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
            # print("n=", n)
            flux, lmax=f_num(w_n[m], w_n[m+1], celll.vr, 1)*1
            lmaxmax=max(lmaxmax, lmax)
            #Update values
            upd[m,:]-=flux
            # print(u[m-1, n+1])
            # print(flux)
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
            #print(dva_n)
            # print(dva_n[1])
            # print(dva_n[2])
            upd, lmax=self.step(dva_n)
            
            cellf=self.mesh.cells[M-1]
            dx=cellf.volume
            # print(lmax)
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
            # print("n=", n+1, ", t=", self.t[n+1])
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
            self.wtot[n]=self.mesh.cells[0].volume*np.sum(self.w[n], axis=0)
            
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
def z(t):
    return 0.9+0.1*np.cos(10*np.pi*t)
def zp(t):
    return -np.pi*np.sin(10*np.pi*t)
def zpp(t):
    return -10*np.pi**2*np.cos(10*np.pi*t)
# def z(t):
#     return 1
# def zp(t):
#     return 0
# def zpp(t):
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
    p=pf(rhoE, rho, u)
    # c=np.sqrt(gamma*p/rho)
    # A=np.abs(np.array([u-v, u+c-v, u-c-v]))
    # lmax=max(A)
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
    # fl, lmaxl = ft(wl, v)
    # fr, lmaxr =ft(wr, v)
    return 0.5*(ft(wl, v)+ft(wr, v))*n-0.5*lamdamax*(wr-wl)*n, lamdamax

# N=5
# M=1000
# mesh=Mesh(0., 1., N)
def step(mesh, w_m, t_m, v_m, a_m):
    N=mesh.N
    
    upd=np.zeros((N, 3))
    lmaxmax=0
        
    cell0=mesh.cells[0]
    flux=np.array([0, (gamma-1)*w_m[0, 2], 0])
    upd[0,:]+=flux
    # print((gamma-1)*u[m-1, 0, 2])
    # print(flux)
    #Loop over the edges
    # upd+=np.array([])
    
    # upd[0:-1]-=np.array([f_num(u[m-1, n], u[m-1, n+1], mesh.cells[n].vr, 1)[0] for n in range(N-1)])
    # upd[1:]+=np.array([f_num(u[m-1, n], u[m-1, n+1], mesh.cells[n].vr, 1)[0] for n in range(N-1)])
    for n in range(N-1):
        celll=mesh.cells[n]
        cellr=mesh.cells[n+1]
        #Approximate the flux by the numerical flux times the edge length
        #Note that celll.vr=cellr.vl
        # print("n=", n)
        flux, lmax=f_num(w_m[n], w_m[n+1], celll.vr, 1)*1
        lmaxmax=max(lmaxmax, lmax)
        #Update values
        upd[n,:]-=flux
        # print(u[m-1, n+1])
        # print(flux)
        upd[n+1,:]+=flux
        
    cellf=mesh.cells[N-1]
    v=zp(t_m)
    # v=u[m-1, N-1, 1]/u[m-1, N-1, 0]
    # v=cellf.vr
    p=pf(w_m[N-1, 2], w_m[N-1, 0], v)
    if p<0:
        print("Negative pressure, p=", p)
    
    rho=w_m[N-1,0]
    c=np.sqrt(gamma*p/rho)
    lmaxmax=max(lmaxmax, c)
    
    forcedensity=w_m[N-1, 0]*cellf.volume*zpp(t_m)
    F=np.array([0, forcedensity, v*forcedensity])
    # print(F)
    flux=np.array([0, p, p*v])
    # flux, lmax=ft(w_m[N-1], v)
    # lmaxmax=max(lmaxmax, lmax)
    # print("flux=", flux)
    upd[N-1,:]-=flux
    upd[N-1,:]+=F
    return upd, lmaxmax
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
    u[0,:,1]=np.sin(x[0])
    # u[0,:,1]=0
    u[0,:,2]=2.5
    utot[0]=mesh.cells[0].volume*np.sum(u[0], axis=0)
    p[0]=pf(u[0,:,2], u[0,:,0], u[0,:,0]*u[0,:,1])
    for m in range(1, M+1):
        lmaxmax=0
        if m % 100==1:
            print("m=", m)
        t=tvec[m]
        
        
        #Old Code
        
        # upd=np.zeros((N, 3))
        
        
        # cell0=mesh.cells[0]
        # # u[m, 0,:]+=dt*f_num(np.array([u[m-1, 0, 0], -u[m-1, 0, 1], u[m-1, 0, 2]]), u[m-1, 0], cell0.vl,1)/cell0.volume
        # # u[m, 0,:]+=dt*f(np.array([u[m-1, 0, 0], 0, u[m-1, 0, 2]]))/cell0.volume
        # # flux=np.array([0, p[m-1, 0], 0])
        # flux=np.array([0, (gamma-1)*u[m-1, 0, 2], 0])
        # upd[0,:]+=flux
        # # print((gamma-1)*u[m-1, 0, 2])
        # # print(flux)
        # #Loop over the edges
        # # upd+=np.array([])
        
        # # upd[0:-1]-=np.array([f_num(u[m-1, n], u[m-1, n+1], mesh.cells[n].vr, 1)[0] for n in range(N-1)])
        # # upd[1:]+=np.array([f_num(u[m-1, n], u[m-1, n+1], mesh.cells[n].vr, 1)[0] for n in range(N-1)])
        # for n in range(N-1):
        #     celll=mesh.cells[n]
        #     cellr=mesh.cells[n+1]
        #     #Approximate the flux by the numerical flux times the edge length
        #     #Note that celll.vr=cellr.vl
        #     # print("n=", n)
        #     flux, lmax=f_num(u[m-1, n], u[m-1, n+1], celll.vr, 1)*1
        #     lmaxmax=max(lmaxmax, lmax)
        #     #Update values
        #     upd[n,:]-=flux
        #     # print(u[m-1, n+1])
        #     # print(flux)
        #     upd[n+1,:]+=flux
            
        # cellf=mesh.cells[N-1]
        # v=zp(tvec[m-1])
        # # v=u[m-1, N-1, 1]/u[m-1, N-1, 0]
        # # v=cellf.vr
        # p_t=pf(u[m-1, N-1, 2], u[m-1, N-1, 0], v)
        # # print("pressure=", p_t)
        
        # force=u[m-1, N-1, 0]**2*cellf.volume*zpp(tvec[m-1])
        # flux=np.array([0, p_t-force, p_t*v])
        # # flux=np.array([0, -force, 0])
        # # print(flux)
        # upd[N-1,:]-=flux
        
        
        
        #New Code
        v=zp(tvec[m-1])
        a=zpp(tvec[m-1])
        upd, lmax=step(mesh, u[m-1], tvec[m-1], v, a)
        cellf=mesh.cells[N-1]
        
        dx=cellf.volume
        cfl=dt/dx
        if cfl*lmaxmax>1:
            print("CFL condition not satisfied")
            
        
        #Initilize the new values to the values at the previous time step
        for k in range(N):
            u[m, k, :]=u[m-1, k, :]*mesh.cells[k].volume
        #Calculate the new cell positions, velocities and volumes.
        
        mesh.update(z(t), dt)
        x[m,:]=mesh.midpoints()
        
        u[m, :, :]+=dt*upd
        for k in range(N):
            u[m, k, :]/=mesh.cells[k].volume
        
        # print("Target velocity: ", zp(t))
        # print("Target velocity: ", cellf.vr)
        # print("Computed velocity: ", u[m, N-1, 1]/u[m, N-1, 0])
        
        #Post computations
        utot[m]=mesh.cells[0].volume*np.sum(u[m], axis=0)
        
        rho=u[m, :, 0]
        u1=u[m, :, 1]/rho
        v=np.array([(mesh.cells[k].vl+mesh.cells[k].vr)/2 for k in range(N)])
        rhoE=u[m, :, 2]
        p[m]=(gamma-1)*rhoE-0.5*(gamma-1)*rho*np.abs(u1)**2
    return tvec, x, u, utot, p
        
def adaptivesolve(M, tf, maxsteps=10000):
    t0=0
    N=maxsteps
    
    mesh=Mesh(0., 1., M)
    x=np.zeros((N+1, M))
    x[0,:]=mesh.midpoints()
    
    dt=np.zeros(N)
    
    t=np.zeros(N+1)
    t[0]=t0
    
    w=np.zeros((N+1, M, 3))
    w[0,:,0]=1
    # w[0,:,1]=np.sin(np.pi*x[0])
    w[0,:,1]=0
    w[0,:,2]=2.5
    
    p=np.zeros((N+1, M))
    p[0]=pf(w[0,:,2], w[0,:,0], w[0,:,0]*w[0,:,1])
    
    wtot=np.zeros((N+1, 3))
    wtot[0]=mesh.cells[0].volume*np.sum(w[0], axis=0)
    n=1
    while n<N+1 and t[n-1]<tf-10**(-15):
        v=zp(t[n-1])
        a=zpp(t[n-1])
        upd, lmax=step(mesh, w[n-1], t[n-1], v, a)
        
        cellf=mesh.cells[M-1]
        dx=cellf.volume
        # dt[n-1]=tf/10000
        dt[n-1]=dx/lmax
        
        t[n]=t[n-1]+dt[n-1]
        if t[n]>tf:
            dt=dt[:n]
            dt[n-1]=tf-t[n-1]
            t=t[:n+1]
            t[n]=tf
            w=w[:n+1,:,:]
            p=p[:n+1,]
            wtot=wtot[:n+1,:]
        print("n=", n, ", t=", t[n])
        #Initilize the new values to the values at the previous time step
        for m in range(M):
            w[n, m, :]=w[n-1, m, :]*mesh.cells[m].volume
        #Calculate the new cell positions, velocities and volumes.
        
        # mesh.update(z(t[n]), dt[n-1])
        mesh.update2(z(t[n]), zp(t[n]))
        x[n,:]=mesh.midpoints()
        
        w[n, :, :]+=dt[n-1]*upd
        for m in range(M):
            w[n, m, :]/=mesh.cells[m].volume
        
        
        #Post computations
        wtot[n]=mesh.cells[0].volume*np.sum(w[n], axis=0)
        
        rho=w[n, :, 0]
        u1=w[n, :, 1]/rho
        # v=np.array([(mesh.cells[k].vl+mesh.cells[k].vr)/2 for k in range(N)])
        rhoE=w[n, :, 2]
        p[n]=(gamma-1)*rhoE-0.5*(gamma-1)*rho*np.abs(u1)**2
        
        n=n+1
    return w, x, t, dt, wtot, p


def convtest():
    #space
    refNvec=np.array([10, 20])
    for refN in refNvec:
        T=1.
        #time
        refM=16000
        t, x, refw, _, _=solve(refN, refM, T)
        Mvec=np.array([125, 250, 500, 1000])
        errvec=np.zeros(4)
        for i in range(len(Mvec)):
            M=Mvec[i]
            t, x, w, _, _=solve(refN, M, T)
            errvec[i]=np.linalg.norm(w[-1,:,:]-refw[-1,:,:])/np.sqrt(refN)
        print("M=", refN)
        print(errvec)
        plt.plot(T/Mvec, errvec)
        plt.xscale("log")
        plt.yscale("log")
    plt.legend(['M=10', 'M=20'])
def convtest2():
    T=1.
    #space
    refN=40
    refMvec=[2000, 4000]
    for refM in refMvec:
        Nvec=np.array([5, 10, 20])
        t, x, refw, _, _=solve(refN, refM, T)
        errvec=np.zeros(len(Nvec))
        for j in range(len(Nvec)):
            N=Nvec[j]
            #time
            # t, x, refw, _, _=solve(refN, refM, T)
            t, x, w, _, _=solve(N, refM, T)
            errvec[j]=np.linalg.norm(w[-1,:,:]-refw[-1,::2**(len(Nvec)-j),:])/np.sqrt(N)
        print("M=", N)
        print(errvec)
        plt.plot(T/Nvec, errvec)
        plt.xscale("log")
        plt.yscale("log")
    plt.legend(['N=2000', 'N=4000'])
# convtest2()
N = 120
# M=2
# # N=100
# # M=40
T=5

# t, x, w, wtot, p=solve(N, M, T)
# dt=T/M

# w, x, t, dt, wtot, p=adaptivesolve(N, T)


# # mesh=Mesh(0., 1., N)
# # fluid=Fluid(mesh)
# # fluid.EEint(T)
# # w2=fluid.w
# # t2=fluid.t


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
# tps=1
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
    
#     ln0.set_data(x[i,:], rho[i, :])
#     ln1.set_data(x[i,:], u[i, :])
#     ln2.set_data(x[i,:], p[i, :])
#     ln3.set_data([x[i,-1], x[i, -1]], [-0.5, 2])
#     # plt.title('t='+str(t[frame]))
#     tt=t[i]
#     time.set_text(f't={tt:.2f}')
#     return ln0, ln1, ln2, ln3, time

        
# ani = FuncAnimation(fig, update, frames=np.array(tps*fps*T),
#                     init_func=init, blit=True, repeat=True)
    
# plt.show()
# f = f"C:\\animation4.mp4" 
# writergif = FFMpegWriter(fps=fps, bitrate=-1)
# ani.save(f, writer=writergif)