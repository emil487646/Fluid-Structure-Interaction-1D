import numpy as np
import matplotlib.pyplot as plt
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


def p(t):
    return 10**5+100*np.sin(10*np.pi*t)
p0=p(0)
ms=100
A=0.1
k=50000
P=0.2
k=ms/(P/2/np.pi)**2
# k=10*np.pi**2*ms
def f(u, t):
    return np.array([u[1], -k/ms*u[0]+A*(p0-p(t))])
def df(u, t):
    return np.array([[0, 1], [-k/ms, 0]])
# def f(u, t):
#     return np.array([u[0]])
# def df(u, t):
#     return np.array([[1.]])
# def exact(t):
#     return np.e**t
def EEstep(f, told, uold, dt):
    return uold+dt*f(uold, told)
def IEstep(f, told, uold, dt):
    def g(u):
        return u-dt*f(u, told+dt)-uold
    def dg(u):
        return np.eye(np.size(uold))-dt*df(u, told+dt)
    u= newton(g, dg, uold, 10**-14)
    # print("tnew=",told+dt)
    # print("u=", u)
    # print("f(u)=", f(u, told+dt))
    # print("g(u)=", g(u))
    # print("dg(u)=", dg(u))
    return u
alpha=1-1./np.sqrt(2)
alphat=2-5./4.*np.sqrt(2)
b1=1-alpha
b2=alpha
bt1=1-alphat
bt2=alphat
def SDIRK2step(f, told, yold, dt):
    S1=yold
    Y1=IEstep(f, told+dt*alpha, S1, dt*alpha)
    K1=f(Y1, told+dt*alpha)
    S2=yold+dt*(1-alpha)*K1
    Y2=IEstep(f, told+dt*alpha, S2, dt*alpha)
    K2=f(Y2, told+dt)
    ynew=S2+dt*alpha*K2
    return ynew
def SDIRK12step(f, told, yold, dt):
    S1=yold
    Y1=IEstep(f, told+dt*alpha, S1, dt*alpha)
    K1=f(Y1, told+dt*alpha)
    S2=yold+dt*(1-alpha)*K1
    Y2=IEstep(f, told+dt*alpha, S2, dt*alpha)
    K2=f(Y2, told+dt)
    ynew=S2+dt*alpha*K2
    ytnew=yold+dt*bt1*K1+dt*bt2*K2
    err=ynew-ytnew
    return ynew, np.linalg.norm(err)
#One step of RK4
# def RK4step(f, told, uold, h):
#     yp1=f(told, uold)
#     yp2=f(told+h/2, uold+h*yp1/2)
#     yp3=f(told+h/2, uold+h*yp2/2)
#     yp4=f(told+h, uold+h*yp3)
#     return uold+h/6*(yp1+2*yp2+2*yp3+yp4)
# #One step of RK34, returns new value and norm of the local error estimate
# def RK34step(f, told, uold, h):
#     yp1=f(told, uold)
#     yp2=f(told+h/2, uold+h*yp1/2)
#     yp3=f(told+h/2, uold+h*yp2/2)
#     zp3=f(told+h, uold-h*yp1+2*h*yp2)
#     yp4=f(told+h, uold+h*yp3)
#     l=h/6*(2*yp2+zp3-2*yp3-yp4)
#     return uold+h/6*(yp1+2*yp2+2*yp3+yp4), norm(l)


#EE integration
def EEint(f, y0, t0, tf, N):
    h=tf/N
    t=t0
    y=np.zeros((N+1, 2))
    tvec=np.zeros((N+1))
    y[0]=y0
    tvec[0]=t0
    for k in range(1, N+1):
        y[k] = EEstep(f, t, y[k-1], h)
        t+=h
        tvec[k]=t
    return tvec, y

#IE integration
def IEint(f, y0, t0, tf, N):
    h=tf/N
    t=t0
    y=np.zeros((N+1, 2))
    tvec=np.zeros((N+1))
    y[0]=y0
    tvec[0]=t0
    for k in range(1, N+1):
        y[k] = IEstep(f, t, y[k-1], h)
        t+=h
        tvec[k]=t
    return tvec, y

#SDIRK2 integration
def SDIRK2int(f, y0, t0, tf, N):
    h=tf/N
    t=t0
    y=np.zeros((N+1, 2))
    tvec=np.zeros((N+1))
    y[0]=y0
    tvec[0]=t0
    for k in range(1, N+1):
        y[k] = SDIRK2step(f, t, y[k-1], h)
        t+=h
        tvec[k]=t
    return tvec, y
def newstep(tol, err, hold):
    return hold*(tol/err)**(1./2.)
#Calculate the new step size
# def newstep(tol, err, errold, hold, k):
#     return (tol/err)**(2/(3*k))*(tol/errold)**(-1/(3*k))*hold
#Adaptive SDIRK12 integration
def SDIRK12int(f, y0, t0, tf, tol, maxsteps=500000):
    N=maxsteps
    h=np.zeros(N)
    t=np.zeros(N)
    y=np.zeros((N, 2))
    h[0]=(tf-t0)*tol**(1/2)/(100*(1+np.linalg.norm(f(y0, t0))))
    t[0]=t0
    y[0]=y0
    errold=tol
    err=errold
    k=1
    while k<N and t[k-1]<tf-10**(-15):
        errold=err
        t[k]=t[k-1]+h[k-1]
        if t[k]>tf:
            h=h[:k+1]
            h[k-1]=tf-t[k-1]
            t=t[:k+1]
            t[k]=tf
            y=y[:k+1]
        unew, err = SDIRK12step(f, t[k-1], y[k-1], h[k-1])
        y[k]=unew
        # h[k]=newstep(tol, err, errold, h[k-1], 2)
        h[k]=newstep(tol, err, h[k-1])
        k=k+1
        
    return t, y, h

# def stifftest():
#     mus=np.array([10, 15, 22, 33, 47, 68, 100, 150, 220, 330])
# #    mus=np.array([10, 15, 22, 33, 47, 68, 100, 150, 220, 330, 470, 680, 1000])
#     its1=np.zeros(len(mus))
#     its2=np.zeros(len(mus))
#     k=0
#     for mu in mus:
#         def f(t, u):
#             return np.array([u[1], mu*(1-u[0]**2)*u[1]-u[0]])
#         u, t=adaptiveRK34(f, [2, 0], 0.0, 0.7*mu, 10**(-6))
#         o=solve_ivp(f, [0, 0.7*mu], [2, 0], method='BDF')
#         its1[k]=len(u)
#         its2[k]=len(o.t)
#         k=k+1
#     plt.plot(mus, its1)
#     plt.plot(mus, its2)
#     plt.xscale('log')
#     plt.yscale('log')
#     plt.show()

# mu=100
# def f(t, u):
#     return np.array([u[1], mu*(1-u[0]**2)*u[1]-u[0]])
# #fig = plt.figure()
# #ax1 = fig.add_subplot(111)
# #ax1.set_xlabel('$t$')
# #ax1.grid()
# #ax1.set_ylabel('$y_2$')
# #t, u=adaptiveRK34(f, [2, 0], 0.0, 2*mu, 10**(-5))
# ##ax1.plot(t, u[:, 1])
# #ax1.plot(u[:, 0], u[:, 1])

# t,u=IEint(f, np.array([0, 0.0]), 0.0, 10., 1000)
# plt.plot(t, u[:,0])

# t,u=SDIRK2int(f, np.array([0., 0.]), 0.0, 10., 1000)
# plt.plot(t, u[:,0])

t,u, h=SDIRK12int(f, np.array([0., 0.]), 0.0, 10., 10**-4)
plt.plot(t, u[:,0], label='u')
# plt.plot(t, u[:,1], label='v')
# plt.plot(t, 50*h, label='50h')
plt.legend()


def errortest():
    intevector=[EEint, IEint, SDIRK2int]
    for inte in intevector:
        M=9
        err=np.zeros(M-2)
        h=np.zeros(M-2)
        for m in range(2, M):
            n=2**m
            t, u=inte(f, np.array([1.]), 0.0, 1., n)
            h[m-2]=1./n
            err[m-2]=np.abs(u[-1,0]-exact(t[-1]))
        plt.plot(h, err)
        plt.xscale("log")
        plt.yscale("log")
        print("Convergence order: ", (np.log(err[-1])-np.log(err[-2]))/(np.log(h[-1])-np.log(h[-2])))
# errortest()
