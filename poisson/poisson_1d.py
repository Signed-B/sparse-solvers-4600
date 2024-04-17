import numpy as np
import matplotlib as plt


a=0;
b=5;
n=10;


f= lambda x: x**2 #forcing
g= lambda x: x #exact



def poisson_solve(a,b,n,f,g):
    #a-b interval
    #number discrete point on interval
    #f forcing function
    #exact for endpoints
    h=(b-a)/(n-1)
    x=np.linspace(a,b,n)
    A=np.zeros([n,n])
    rhs=np.zeros([n,1])
    i=0
    A[i,i]=1
    rhs[i]=g(x[i])
    for i in range(1,n-1):
        A[i,i-1]=-1.0/(h**2);
        A[i,i]=2/(h**2);
        A[i,i+1]=-1/(h**2);
        rhs[i]=f(x[i]);
    i=n-1
    A[i,i]=1
    rhs[i]=g(x[i])
    u=np.linalg.solve(A, rhs)
    return x,u,A
        
x,u,A=poisson_solve(a, b, n, f, g)