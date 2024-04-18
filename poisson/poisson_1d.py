import numpy as np
import matplotlib.pyplot as plt


a=-2
b=2
n=50



f= lambda x: x**2 #forcing
g= lambda x: x #exact for endpoints


#u''(x)=f(x)
def poisson_solve(a,b,n,f,g):
    #a-b interval
    #n-number discrete point on interval
    #f forcing function
    #exact function to calc endpoints
    #c_1, c_2 are coefficients on diff eq
    h=(b-a)/(n-1)
    x=np.linspace(a,b,n)
    A=np.zeros([n,n])
    rhs=np.zeros([n,1])
    i=0
    A[i,i]=1
    rhs[i]=g(x[i])
    for i in range(1,n-1):
        A[i,i-1]=1.0/(h**2)
        A[i,i]=-2/(h**2)
        A[i,i+1]=1/(h**2)
        rhs[i]=f(x[i])
    i=n-1
    A[i,i]=1
    rhs[i]=g(x[i])
    u=np.linalg.solve(A, rhs)
    return x,u,A
        
x,u,A=poisson_solve(a, b, n, f, g)
plt.plot(x, u)
plt.plot(x,g(x))
plt.show()