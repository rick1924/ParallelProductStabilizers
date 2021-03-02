import numpy as np
from numpy import genfromtxt
import os
import time

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
my_mat_path = os.path.join(THIS_DIR, "Inputs/starting_tableau.csv")
my_vec_path = os.path.join(THIS_DIR, "Inputs/starting_phases.csv")


def alg1(mat, phases):
    """Sequential algorithm 1. Takes as input a rectangular matrix nx2n and a vector of phases of size n. Computes the canonical form of the tableau."""

    i = 0
    # Algorithm applied to the right-hand side of the matrix
    for j in range(num_cols // 2):
        k_exists = False
        # Find the first "k" that has a 1 in the "j+n"th column
        for k in range(i, num_rows):
            if mat[k][j + n] == 1:
                k_exists = True
                break
        # If "k" exists we swap rows and row-reduce
        if k_exists == True:
            mat[[i, k]] = mat[[k, i]]
            phases[i], phases[k] = phases[k], phases[i]
            z_i, x_i = mat[i, :n], mat[i, n:]
            # Row reduce
            for m in range(num_rows):
                if mat[m][j + n] == 1 and m != i:
                    z_m, x_m = mat[m, :n], mat[m, n:]
                    partial_symplectic = z_m @ x_i - x_m @ z_i
                    unmodded_row = mat[m] + mat[i]
                    mat[m] = (mat[m] + mat[i]) % 2
                    exp_mod_4 = (
                        partial_symplectic
                        + (
                            mat[m, :n] @ mat[m, n:]
                            - unmodded_row[:n] @ unmodded_row[n:]
                        )
                    ) % 4
                    phases[m] = (exp_mod_4 // 2 + phases[m] + phases[i]) % 2
            i += 1
    # Algorithm applied to the left-hand side of the matrix. Similar thing as above
    for j in range(num_cols // 2):
        k_exists = False
        for k in range(i, num_rows):
            if mat[k][j] == 1:
                k_exists = True
                break
        if k_exists == True:
            mat[[i, k]] = mat[[k, i]]
            phases[i], phases[k] = phases[k], phases[i]
            z_i, x_i = mat[i, :n], mat[i, n:]
            for m in range(num_rows):
                if mat[m][j] == 1 and m != i:
                    z_m, x_m = mat[m, :n], mat[m, n:]
                    partial_symplectic = z_m @ x_i - x_m @ z_i
                    unmodded_row = mat[m] + mat[i]
                    mat[m] = (mat[m] + mat[i]) % 2
                    exp_mod_4 = (
                        partial_symplectic
                        + (
                            mat[m, :n] @ mat[m, n:]
                            - unmodded_row[:n] @ unmodded_row[n:]
                        )
                    ) % 4
                    phases[m] = (exp_mod_4 // 2 + phases[m] + phases[i]) % 2
            i += 1
    return mat, phases


def alg2(mat, phases):
    """Sequential algorithm 2. Takes as input the canonical reduced form of the stabilizer tableau computed in algorithm 1. Outputs a quantum circuit of the form H-CX-CZ-S-H, the matrix in basis state form, and the vector of phases."""

    circuit = []
    i = 0
    # First H BLOCK
    for j in range(n, 2 * n):
        k_exists = False
        for k in range(i, n):
            if mat[k, j] == 1:  # If entry is an X, Y, or Z literal
                mat[[k, i]] = mat[[i, k]]
                phases[i], phases[k] = phases[k], phases[i]
                k_exists = True
                break
        if k_exists != True:
            # Search from the back. On the left hand side (look for Z literals)
            for k2 in range(i, n)[::-1]:
                if mat[k2, j - n] == 1 and mat[k2, j] == 0:
                    mat[[i, k2]] = mat[[k2, i]]
                    phases[i], phases[k2] = phases[k2], phases[i]
                    for jj in range(j + 1, 2 * n):
                        if mat[i, jj] == 1 or mat[i, jj - n] == 1:
                            phases = (phases + mat[:, jj] * mat[:, jj - n]) % 2
                            mat[:, (jj, jj - n)] = mat[:, (jj - n, jj)]
                            circuit.append(("H", jj - n))
                    break
        i += 1
    # CX BLOCK
    for j in range(n, 2 * n):
        for k in range(j + 1, 2 * n):
            if mat[j - n, k] == 1:  # If entry is an X, Y, or Z literal
                phases = phases + mat[:, k - n] * mat[:, j] * (
                    mat[:, j - n] + mat[:, k] + 1
                )
                mat[:, j - n] = (mat[:, j - n] + mat[:, k - n]) % 2
                mat[:, k] = (mat[:, k] + mat[:, j]) % 2
                circuit.append(("CX", j - n, k - n))
    # CZ BLOCK
    for j in range(0, n):
        for k in range(j + 1, n):
            if mat[j, k] == 1 and mat[j, k + n] == 0:  # If entry is a Z literal
                phases = phases + mat[:, j + n] * mat[:, k + n]
                mat[:, j] = (mat[:, j] + mat[:, k + n]) % 2
                mat[:, k] = (mat[:, k] + mat[:, j + n]) % 2
                circuit.append(("CZ", j, k))
    # S BLOCK
    for j in range(n, 2 * n):
        if mat[j - n, j] == 1 and mat[j - n, j - n] == 1:  # If entry is a Y literal
            phases = (phases + mat[:, j] * mat[:, j - n]) % 2
            mat[:, j - n] = (mat[:, j - n] + mat[:, j]) % 2
            circuit.append(("S", j - n))
    # H BLOCK 2
    for j in range(n, 2 * n):
        if mat[j - n, j] == 1 and mat[j - n, j - n] == 0:  # If entry is a Z literal
            phases = (phases + mat[:, j] * mat[:, j - n]) % 2
            mat[:, [j - n, j]] = mat[:, [j, j - n]]
            circuit.append(("H", j - n))
    # Eliminate trailing Z literals
    for j in range(0, n):
        z_j, x_j = mat[j, :n], mat[j, n:]
        for k in range(j + 1, n):
            if mat[k, j] == 1 and mat[k, j + n] == 0:  # If entry is a Z literal
                z_k, x_k = mat[k, :n], mat[k, n:]
                partial_symplectic = z_k @ x_j - x_k @ z_j
                unmodded_row = mat[j, :] + mat[k, :]
                mat[k, :] = (mat[j, :] + mat[k, :]) % 2
                exp_mod_4 = (
                    partial_symplectic
                    + (mat[k, :n] @ mat[k, n:] - unmodded_row[:n] @ unmodded_row[n:])
                ) % 4
                phases[k] = (exp_mod_4 // 2 + phases[k] + phases[j]) % 2
    return circuit, mat, phases


# Main
n = 8
mat = genfromtxt(my_mat_path, dtype="int", delimiter=",")
mat = mat.reshape((n, 2 * n))
num_rows, num_cols = mat.shape
phases = genfromtxt(my_vec_path, dtype="int", delimiter=",")

start_time = time.time()
mat, phases = alg1(mat, phases)
circuit, mat, phases = alg2(mat, phases)
end_time = time.time()

print(end_time - start_time)
