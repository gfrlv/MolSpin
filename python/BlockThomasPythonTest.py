from typing import Union, Tuple
import numpy as np
import scipy.sparse as sp
from scipy.sparse import linalg as spl


def GetBlockMatrix(A : sp.csc_matrix, r0 : int, r1 : int, c0 : int, c1 : int) -> sp.csc_matrix:
    blk = A[r0:r1+1, c0:c1+1]
    return blk

def GetBlockVector(A : np.ndarray, r0 : int, r1 : int) -> np.ndarray:
    blk = A[r0:r1+1]
    return blk

def AugmentedMatrix(Mat: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Return [Mat | b] as a dense array."""
    b = np.asarray(b).reshape(-1, 1)
    return np.hstack([Mat, b])

def UndoAugmentedMatrix(AugMat: np.ndarray) -> Tuple[sp.csc_matrix, np.ndarray]:
    """Split [Mat | b] into (Mat, b) where b is a column vector returned as 1D array."""
    rows, cols = AugMat.shape
    Mat = AugMat[:, : cols - 1]
    Mat = sp.csc_matrix(Mat)
    b = AugMat[:, cols - 1]
    return Mat, b
    

def ThomasBlockSolver(A : sp.csc_matrix, b : np.ndarray, block_size : int) -> np.ndarray:

    N = A.shape[0]
    if N % block_size != 0:
        raise ValueError("block_size must divide matrix size")
    n_blocks = N // block_size

    A_blocks = []
    B_blocks = []

    for i in range(0, n_blocks):
        B_blocks.append(GetBlockVector(b, i*block_size, (i+1)*block_size - 1))

        #Off-diagonal LEFT
        if i > 0:
            A_blocks.append(GetBlockMatrix(A, i*block_size, (i+1)*block_size - 1,
                                     (i-1)*block_size, i*block_size - 1))
        #Diagonal block
        A_blocks.append(GetBlockMatrix(A, i*block_size, (i+1)*block_size - 1,
                                 i*block_size, (i+1)*block_size - 1))
        #Off-diagonal RIGHT
        if i < n_blocks - 1:
            A_blocks.append(GetBlockMatrix(A, i*block_size, (i+1)*block_size - 1,
                                     (i+1)*block_size, (i+2)*block_size - 1))
    
    #if one of the diagonal blocks is zero (A_blocks[3*i]) is zero, pivot the matrix
    #Permutation matrices
    
    #print(A.todense())
    rowswaplst = []
    colswaplst = []
    for i in range(n_blocks):
        D_block = A_blocks[3*i]
        if spl.norm(D_block) == 0:
            #FULL PIVOT
            #find a non-zero block to swap with, ideally from either the left or right, as this will maintain the tridiagonal form
            swap_found = False
            #For a column swap, check either a left or right block first
            offset = 0
            operations = []
            for offsets in [-1,1]:
                if swap_found:
                    break
                if not (0 <= i + offsets < n_blocks):
                    operations.append(0)
                    continue
                if CheckColSwap(A_blocks, i, offsets):
                    operations.append(1)
                else:
                    operations.append(0)
                if CheckRowSwap(A_blocks, i, offsets):
                    operations.append(1)
                else:
                    operations.append(0)
            op = 0
            for valid in operations:
                if not valid:
                    op += 1
                if valid:
                    break
            colswap = 1 if op % 2 == 0 else 0
            offset = -1 if op - 2 < 0 else 1
            BlockIndicies = []
            if colswap:
                BlockIndicies = [-2,0,2]
                colswaplst.append([i,i-offset])
            else:
                BlockIndicies = [-1,0,1]
                rowswaplst.append[i,i-offset]

            #swapblocks
            BlocksToSwap = []
            if (3*i)-2 > 0:
                BlocksToSwap.append([(3*i - 2), -1 if offset == -1 else (3*i - 3)])
            if (3*i + 2) < len(A_blocks):
                BlocksToSwap.append([(3*i + 2), -1 if offset == 1 else (3*i + 3)])
            BlocksToSwap.append([(3*i), (3*i + offset)])

            print(BlocksToSwap)
            for Blocks in BlocksToSwap:
                if Blocks[1] != -1:
                    A_blocks[Blocks[0]], A_blocks[Blocks[1]] = A_blocks[Blocks[1]], A_blocks[Blocks[0]]
                else:
                    A_blocks[Blocks[0]] = A_blocks[Blocks[0]]#] sp.csc_matrix(A_blocks[0].shape, dtype=np.complex128)
            
    
    for rowswap in rowswaplst:
        B_blocks[rowswap[0]], B_blocks[rowswap[1]] = B_blocks[rowswap[1]], B_blocks[rowswap[0]]

    # Forward elimination
    for i in range(1, n_blocks):
        #Get Blocks
        D_prev = A_blocks[3*(i-1)]
        U_prev = A_blocks[(3*(i-1)) + 1]
        #L_prev = A_blocks[(3*(i-1)) - 1]
        B_prev = B_blocks[i - 1]

        D_curr = A_blocks[3*i]
        #U_curr = A_blocks[(3*i) + 1]
        L_curr = A_blocks[(3*i) - 1]
        B_curr = B_blocks[i]

        #Step 1: Form Augmented Matrix and solve
        AugMat = AugmentedMatrix(U_prev.toarray(), B_prev)
        X = np.linalg.solve(D_prev.toarray(), AugMat)
        #Step 2: Update current blocks
        #B_blocks[i-1] = UndoAugmentedMatrix(X)[1]
        AugMat = AugmentedMatrix(D_curr.toarray(), B_curr)
        AugMatNew = AugMat - L_curr.toarray() @ X
        D_new, B_new = UndoAugmentedMatrix(AugMatNew)
        A_blocks[3*i] = D_new
        B_blocks[i] = B_new
        #A_blocks[3*(i-1)] = sp.csc_matrix(sp.identity(block_size))
        A_blocks[(3*(i) - 1)] = sp.csc_matrix(np.zeros((block_size, block_size)))

    # Back substitution
    X_blocks = []
    for i in range(n_blocks - 1, -1, -1):
        D_curr = A_blocks[3*i]
        B_curr = B_blocks[i]

        if i == n_blocks - 1:
            #Last block
            X_curr = np.linalg.solve(D_curr.toarray(), B_curr)
            X_blocks.insert(0, X_curr)
        else:
            U_curr = A_blocks[(3*i) + 1]
            X_next = X_blocks[0]

            value = B_curr - U_curr.toarray() @ X_next
            X_curr = np.linalg.solve(D_curr.toarray(), value)
            #X_curr = value
            X_blocks.insert(0, X_curr)
    #Combine X_blocks into single solution vector
    X_full = np.zeros(N, dtype=b.dtype)

    for colswap in colswaplst:
        X_blocks[colswap[0]], X_blocks[colswap[1]] = X_blocks[colswap[1]], X_blocks[colswap[0]]

    for i in range(n_blocks):
        X_full[i*block_size:(i+1)*block_size] = X_blocks[i]
    return X_full

def CheckColSwap(Blocks, i, offset) -> bool:
    #check block either above or below
    if 0 <= (3*i + 2 * offset) < len(Blocks):
        A1 = Blocks[3*i + 2 * offset]
        if spl.norm(A1) != 0:
            return False
    if 0 <= (3*i - 5 * offset) < len(Blocks):
        A2 = Blocks[3*i - 5 * offset]
        if spl.norm(A2) != 0:
            return False
    return True

def CheckRowSwap(Blocks, i, offset) -> bool:
    #check block either left or right
    A1 = Blocks[3*i + offset]
    if spl.norm(A1) != 0:
        return False
    A2 = Blocks[3*i + 4 * offset]
    if spl.norm(A2) != 0:
        return False
    return True

if __name__ == '__main__':
    #np.random.seed(0)
    bs = 3
    nblk = 5
    N = bs * nblk
    Ad = np.zeros((N, N), dtype=np.complex128)
    for i in range(nblk):
        Ad[i * bs:(i + 1) * bs, i * bs:(i + 1) * bs] = np.eye(bs) * (5.0 + 0.1j) + 0.1 * (np.random.randn(bs, bs) + 1j * np.random.randn(bs, bs))
        if i < nblk - 1:
            Ad[i * bs:(i + 1) * bs, (i + 1) * bs:(i + 2) * bs] = 0.1 * (np.random.randn(bs, bs) + 1j * np.random.randn(bs, bs))
            Ad[(i + 1) * bs:(i + 2) * bs, i * bs:(i + 1) * bs] = 0.1 * (np.random.randn(bs, bs) + 1j * np.random.randn(bs, bs))

    A = sp.csc_matrix(Ad)
    b = np.random.randn(N) + 1j * np.random.randn(N)

    #make one of the diagonal blocks a zero
    print(A)
    A[2*bs:3*bs, 2*bs:3*bs] = sp.csc_matrix(np.zeros((bs, bs)))
    print(A)
    A[1*bs:2*bs, 2*bs:3*bs] = sp.csc_matrix(np.zeros((bs, bs)))
    A[4*bs:5*bs, 3*bs:4*bs] = sp.csc_matrix(np.zeros((bs, bs)))
    #make a entire row of blocks zero to test pivoting
    #A[2*bs:3*bs, :] = sp.csc_matrix(np.zeros((bs, N)))
    print(A)
    x_cpp_like = ThomasBlockSolver(A, b, bs)
    print(x_cpp_like)
    x_direct = np.linalg.solve(Ad, b)
    print(x_direct)
    print('norm diff (cpp-like):', np.linalg.norm(x_cpp_like - x_direct))