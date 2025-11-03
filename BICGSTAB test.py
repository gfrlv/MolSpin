#BICGSTAB test

import numpy as np
import scipy.sparse as sp
import scipy.linalg as spl

from scipy.sparse import identity
from scipy.sparse import diags
from scipy.sparse.linalg import onenormest

def bicgstab(A, b, K = None, x0=None, tol=1e-10, max_iter=1000):
    '''
    Preconditioned Biconjugate Gradient Stabilized Method for solving Ax = b.
    
    Parameters:
    A : 2D array-like
        Coefficient matrix.
    B : 1D array-like
        Right-hand side vector.
    K : 2D array-like
        Preconditioner matrix (optional).
        If not provided, diagonal of A is used as preconditioner.
    x0 : 1D array-like, optional
        Initial guess for the solution.
    tol : float, optional
        Tolerance for convergence.
    max_iter : int, optional
        Maximum number of iterations.
    
    Returns:
    x : 1D array
        Approximate solution vector.
    '''
    
    A = np.array(A)
    b = np.array(b)

    if K is None:
        K = np.diag(np.diag(A))
    
    K_1, K_2 = spl.lu(K, permute_l=True)
    K_1_inv = spl.solve(K_1, np.eye(K.shape[0]))
 
    if x0 is None:
        x0 = np.zeros_like(b)
    
    x = x0.copy()
    r_naught = b - A @ x
    r_naught_hat = r_naught.copy()
    r = r_naught.copy()
    rho_naught = r_naught.copy()
    rho_prev = rho_naught.copy()
    rho_k1 = np.dot(rho_prev,rho_prev)

    for k in range(1,max_iter+1):
        y = np.linalg.solve(K,rho_prev)
        v = A @ y
        alpha =  rho_k_1 / np.dot(r_naught_hat, v)
        h = x + alpha * y
        s = r - alpha * v
        if np.linalg.norm(s) < tol:
            print(f"Converged in {k} iterations.")
            return h
        z = np.linalg.solve(K, s)
        t = A @ z
        omega = np.dot(K_1_inv @ t, K_1_inv @ s) / np.dot(K_1_inv @ t, K_1_inv @ t)
        x = h + omega * z
        r = s - omega * t
        if np.linalg.norm(r) < tol:
            print(f"Converged in {k} iterations.")
            return x
        rho_k = np.dot(r_naught_hat, r)
        beta = (rho_k / rho_k_1) * (alpha / omega)
        rho_prev = r + beta * (rho_prev - omega * v)
        rho_k_1 = rho_k
    
    return x

def spai(A,m):
    n = A.shape[0]
    condition_numbers = []

    if A is not sp.csr_matrix:
        A = sp.csr_matrix(A)

    identity_matrix = identity(n, format='csr')
    alpha = 2 / onenormest(A @ A.T)
    M = alpha * A
    #condition_number = np.linalg.cond(A.todense() @ M.todense())
    for i in range(m):
        C = A @ M
        G = identity_matrix - C
        AG = A @ G
        trace = (G.T @ AG).diagonal().sum()
        alpha = trace / np.linalg.norm(AG.data)**2
        M = M + alpha * G
        #do a check to see if the condition number is lower
        #condition_number_2 = np.linalg.cond(A.todense() @ M.todense())
        #if condition_number_2 < condition_number:
        #    condition_number = condition_number_2
        #else:
        #    print(f"Condition number did not improve at iteration {i+1}. Stopping early.")
        #    #break
    return M

def IncompleteBICGSTAB(A, max_iter = 5):
    n = A.shape[0]
    identity_matrix = identity(n, format='csr')
    x = []
    for i in range(0,n):
        col = identity_matrix.getcol(i)
        colarray = []
        for a in col:
            colarray.append(a.toarray()[0][0])
        col = np.array(colarray)
        x0 = np.zeros(n)
        returnvec = bicgstab(A,col,K = identity_matrix.todense(), x0=x0, tol=1e-10, max_iter=max_iter)
        x.append(returnvec)
    x = np.array(x)
    M = sp.csr_matrix(x.T)
    print(M)
    return M


# Example usage
if __name__ == "__main__":
    #A = np.array([[1, 1, 1], [6, -4, 5], [5, 2, 2]], dtype=float)
    A = np.array([[2,0,4,0,1,0,0,0],
                  [0,3,0,0,0,1,0,0],
                  [2,0,1,3,0,0,1,0],
                  [0,4,0,1,0,0,0,1],
                  [1,0,0,0,3,0,1,0],
                  [0,1,0,0,2,1,0,0],
                  [0,0,1,0,0,3,2,0],
                  [0,0,0,1,0,0,2,4]],dtype=float)
    print(np.linalg.cond(A))
    #M = spai(A, 100)

    from matplotlib import pyplot as plt

    #fig = plt.figure(figsize=(8 ,8))
    #ax = fig.add_subplot(111)
    #ax.spy(M, markersize=1)
    #plt.show()
    
    #print(np.linalg.cond(A @ M))

    n = 1000

    data = [2.001 * np.ones(n),
        -1. * np.ones(n - 1),
        -1. * np.ones(n - 1)]

    offsets = [0, 1, -1]

    A2  = diags(data, offsets=offsets, shape=(n, n), format='csr')

    print(np.linalg.cond(A2.todense()))
    M = spai(A2, 50)
    fig = plt.figure(figsize=(8 ,8))
    ax = fig.add_subplot(111)
    ax.spy(M, markersize=1)
    plt.show()
    print(np.linalg.cond(A2.todense() @ M.todense()))

    M = IncompleteBICGSTAB(A,10)
    print(np.linalg.cond(A@M))


    b = np.array([4, 2, 0,-1,0,1,0.5,-2], dtype=float)
    x0 = np.zeros_like(b)
    #K = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)  # Identity preconditioner
    x = bicgstab(A, b, x0=x0, tol=1e-10, max_iter=1000)
    print("Solution:", x)
    print("Residual:", np.linalg.norm(A @ x - b))
    print("Expected:", np.linalg.solve(A, b))
    print("Expected Residual:", np.linalg.norm(np.linalg.solve(A, b) - x))
    print("Convergence:", np.linalg.norm(A @ x - b) < 1e-10)
    print("Iterations:", "Converged" if np.linalg.norm(A @ x - b) < 1e-10 else "Did not converge")
                 


