import numpy as np
import matplotlib.pyplot as plt
import bisect

# from scipy.optimize import newton

def newton(g, dg, y0, eps, max_iter=50):
    yn = y0
    for n in range(0,max_iter):
        gyn = g(yn)
        nor=np.linalg.norm(gyn)
        if nor < eps:
            # print('Found solution after',n,'iterations.')
            return yn
        
        dgyn=dg(yn)
        if np.linalg.norm(dgyn, 'fro')< eps:
            print('Zero derivative. No solution found.')
            return None
        yn = yn - np.linalg.solve(dgyn, gyn)
    print('Exceeded maximum iterations. No solution found.')
    return yn

alpha=1-1./np.sqrt(2)
alphat=2-5./4.*np.sqrt(2)
b1=1-alpha
b2=alpha
bt1=1-alphat
bt2=alphat

def newstep(tol, err, hold):
    return hold*(tol/err)**(1./2.)
def interpolate(t, f, s):
    index=max(bisect.bisect_left(t, s)-1,0)
    h=s-t[index]
    return f[index]+(f[index+1]-f[index])/(t[index+1]-t[index])*(s-t[index])

class Piston:
    def __init__(self, y0, tp, maxsteps=50000):
        N=self.N=maxsteps
        self.h=np.zeros(N)
        self.t=np.zeros(N)
        self.y=np.zeros((N, 2))
        self.y[0]=y0
        self.tp=tp
        
        self.p0=50
        self.ms=100
        self.A=0.1
        self.k=5000
    def interpolateinput(self, t):
        return interpolate(self.tp[0], self.tp[1], t)
    def F(self, t):
        t=min(t, np.pi)
        return 100*np.sin(t)
    def f(self, u, t): 
        pt=self.interpolateinput(t)
        # print(pt)
        return np.array([u[1], -self.k/self.ms*u[0]+self.A*(self.p0-pt)])
    def df(self, u, t):
        return np.array([[0, 1], [-self.k/self.ms, 0]])
    def IEstep(self, told, uold, dt):
        def g(u):
            return u-dt*self.f(u, told+dt)-uold
        def dg(u):
            return np.eye(np.size   (uold))-dt*self.df(u, told+dt)
        u= newton(g, dg, uold, 10**-14)
        return u
    def SDIRK12step(self, told, yold, dt):
        S1=yold
        Y1=self.IEstep(told, S1, dt*alpha)
        K1=self.f(Y1, told+dt*alpha)
        S2=yold+dt*(1-alpha)*K1
        Y2=self.IEstep(told+dt*(1-alpha), S2, dt*alpha)
        K2=self.f(Y2, told+dt)
        ynew=S2+dt*alpha*K2
        ytnew=yold+dt*bt1*K1+dt*bt2*K2
        err=ynew-ytnew
        return ynew, np.linalg.norm(err)

    def SDIRK12int(self, te, tol):
        N=self.N
        self.h[0]=tol#(te-self.t[0])*tol**(1/2)/(100*(1+np.linalg.norm(self.f(self.y[0], self.t[0]))))
        self.t[0]=self.t[0]
        errold=tol
        err=errold
        k=1
        while k<N and self.t[k-1]<te-10**(-15):
            errold=err
            self.t[k]=self.t[k-1]+self.h[k-1]
            if self.t[k]>te:
                self.h=self.h[:k+1]
                self.h[k-1]=te-self.t[k-1]
                self.t=self.t[:k+1]
                self.t[k]=te
                self.y=self.y[:k+1]
            unew, err = self.SDIRK12step(self.t[k-1], self.y[k-1], self.h[k-1])
            self.y[k]=unew
            # h[k]=newstep(tol, err, errold, h[k-1], 2)
            self.h[k]=newstep(tol, err, self.h[k-1])
            k=k+1
    def SDIRK2int(self, te, N):
        self.N=N
        self.h[0]=te/N
        for k in range(1, N):
            self.t[k]=self.t[k-1]+self.h[k-1]
            self.y[k], _ = self.SDIRK12step(self.t[k-1], self.y[k-1], self.h[k-1])
            self.h[k]=self.h[k-1]
    
    def output(self):
        a=np.array([self.f(self.y[n], self.t[n])[1] for n in range(len(self.t))])
        return np.stack((self.t, self.y[:,0], self.y[:,1], a))
    
