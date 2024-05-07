import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter_ns
# import scipy.linalg.lu as lu
import scipy as sp

import matplotlib

matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['figure.dpi'] = 300

#Thomas algo stuff because it's in a different folder :/

def thomas(A: np.array, d: np.array, *args) -> np.array:
    """
    Solves a tridiagonal system of equations using the Thomas algorithm by
    parsing the matrix A into its three diagonals.
    """
    return abc_thomas(
        np.concatenate([np.array([0]), np.diag(A, -1)]),
        np.diag(A),
        np.concatenate([np.diag(A, 1), np.array([0])]),
        d
    )

def abc_thomas(a: np.array, b: np.array, c: np.array, d: np.array) -> np.array:
    """
    Solves a tridiagonal system of equations using the Thomas algorithm 
    given the three diagonals of the matrix A.
    """
    # TODO rewrite to make more verbose.
    # TODO comment the code.
    n = len(d)
    c_ = np.zeros(n)
    d_ = np.zeros(n)
    x = np.zeros(n)
    
    c_[0] = c[0] / b[0]
    d_[0] = d[0] / b[0]
    
    for i in range(1, n):
        c_[i] = c[i] / (b[i] - a[i] * c_[i - 1])
        d_[i] = (d[i] - a[i] * d_[i - 1]) / (b[i] - a[i] * c_[i - 1])
        
    x[n - 1] = d_[n - 1]
    
    for i in range(n - 2, -1, -1):
        x[i] = d_[i] - c_[i] * x[i + 1]
        
    return x

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

a=0
b=1

f= lambda x: 12*x**2+8 #forcing
g= lambda x: x**4+4*x**2-6*x+2 #exact for endpoints

n = np.arange(10,1001,10)
time_lusolve = np.zeros([n.size,1])
time_thomas = np.zeros([n.size,1])
for i in n:
    #u''(x)=f(x)
            
    x,rhs,A=poisson(a, b, i, f, g)
    
    #If we use linalg.solve, we need to understand how it works. Currently I think it uses an LU factorization
    start_lusolve = perf_counter_ns()
    p, l, u = sp.linalg.lu(A)
    y = sp.linalg.solve_triangular(l, p@rhs, lower=True)
    lusoln = sp.linalg.solve_triangular(u, y)
    time_lusolve[int(i/10)] = perf_counter_ns() - start_lusolve
#    print('LU decomposition time: ', time_lusolve, 'ns')
    
    
    start_thomas = perf_counter_ns()
    thomassoln = thomas(A,rhs) #fix to correct solver
    #uthomas = uthomas[:,0]
    time_thomas[int(i/10)] = perf_counter_ns() - start_thomas
#    print('Thomas time: ', time_thomas, 'ns')

plt.title('Thomas vs LU Decomposition')
plt.xlabel('n')
plt.ylabel('Time (ns)')
plt.plot(n, time_lusolve[:,0], label = 'LU Decomposition')
plt.plot(n, time_thomas[:,0], label = 'Thomas')
plt.legend()
plt.show()

#Probably want to implement a way to solve for multiple n's and then visualize the difference in times as n increases.
#Some n's will be too small (e.g. n<500) but I don't know when you start seeing differences
n_test = 500
x_test = np.linspace(a,b,n_test)

interp_thomas = np.zeros([n_test,1])
for i in range(0,n_test):
    interp_thomas[i] = np.interp(x_test[i],x,thomassoln)
    
g = (g(x_test)) 
error = np.zeros([n_test,1])   
for i in range(0,n_test):
    error[i] = abs(interp_thomas[i] - g[i])

#plot solutions
plt.plot(x,thomassoln,'--') #calculated
plt.plot(x_test,g)
plt.title('Actual Function vs Finite Diff')
plt.legend(["Finite","Actual"])
plt.show()

#plot log error
plt.plot(x_test,np.log10(error))
plt.title('Log10 Error')
plt.show()
