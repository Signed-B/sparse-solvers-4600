import numpy as np
import matplotlib.pyplot as plt
from time import process_time

n=10 #number points per direction
x1=0
x2=1
y1=0
y2=1



def f(x,y): #forcing
    return x*y

def g(x,y): #end conditions function
    return 2

def poisson_solve(A,F): #solve linear system given banded matrix A and vecotr F
    U=np.linalg.solve(A,F)
    return U

def poisson_setup(x1,x2,y1,y2,n,f,g): #2d finite difference setup
    #square bounds only
    xrange=np.linspace(x1,x2,n)
    yrange=np.linspace(y1,y2,n)
    n2=n
    n=n**2
    A=np.zeros([n,n])
    F=np.zeros(n)
    h=(x2-x1)/(n2)
    A[0,0]=4
    A[1,1]=4
    A[2,2]=4
    A[n-1,n-1]=4
    A[n-2,n-2]=4
    A[n-3,n-3]=4
    A[1,0]=-1
    A[0,1]=-1
    A[2,1]=-1
    A[1,2]=-1
    A[2,3]=-1
    A[0,3]=-1
    A[1,4]=-1
    A[2,5]=-1
    A[n-1,n-2]=-1
    A[n-2,n-1]=-1
    A[n-2,n-3]=-1
    A[n-3,n-2]=-1
    A[n-3,n-4]=-1
    A[n-1,n-4]=-1
    A[n-2,n-5]=-1
    A[n-3,n-6]=-1
    for i in range(3,n-3):
        A[i,i]=4
        A[i,i+1]=-1
        A[i,i+3]=-1
        A[i,i-1]=-1
        A[i,i-3]=-1
    A=np.transpose(A)*(1/(h**2))
    for i in range(n2):
        for j in range(n2):
            F[i*n2+j]=g(xrange[i],yrange[j])
    for i in range(1,n2-1):
        for j in range(1,n2-1):
            F[i*n2+j]=f(xrange[i],yrange[j])


    return A,F
        
A,F=poisson_setup(x1,x2,y1,y2,n,f,g)
start_time=process_time()
U=poisson_solve(A, F)
end_time=process_time()
run_time=end_time-start_time

solution_field=np.zeros([n,n])
for i in range(n):
    for j in range(n):
        solution_field[i,j]=U[i*n+j]
    
    
    
    
    
    
    
    
    
