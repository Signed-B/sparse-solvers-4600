import numpy as np
import matplotlib.pyplot as plt
from time import process_time
from thomas import thomas

a=-2
b=2
n=750



f= lambda x: x**2 #forcing
g= lambda x: x #exact for endpoints


#u''(x)=f(x)
def poisson(a,b,n,f,g):
    #a-b interval
    #n-number discrete point on interval
    #f forcing function
    #exact function to calc endpoints
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
    return x,rhs,A
        
x,rhs,A=poisson(a, b, n, f, g)

#If we use linalg.solve, we need to understand how it works. Currently I think it uses an LU factorization
start_npsolve = process_time()
unpsolve = np.linalg.solve(A,rhs)
time_npsolve = process_time() - start_npsolve
print('Np.linalg.solve time: ', time_npsolve)


start_thomas = process_time()
uthomas = np.linalg.solve(A,rhs)
time_thomas = process_time() - start_thomas
print('Thomas time: ', time_thomas)

#Probably want to implement a way to solve for multiple n's and then visualize the difference in times as n increases.
#Some n's will be too small (e.g. n<500) but I don't know when you start seeing differences
'''
plt.plot(x,u)
plt.plot(x,g(x))
plt.show()
'''