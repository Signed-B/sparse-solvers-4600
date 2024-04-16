import numpy as np
# import scipy as sp

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
