import numpy as np


def check_C_diagnolize_S(S,C):
    S_diag = {}
    for key in C:
        S_diag[key] = np.dot(C[key].T,np.dot(S[key],C[key]))
    max_error = 0
    for key in S_diag:
        for i in range(0, len(S_diag[key])):
            max_error = np.max((np.abs(S_diag[key][i,i] - 1), max_error))
            for j in range(0, i):
                max_error = np.max((np.abs(S_diag[key][i,i] - 1), max_error))
    return max_error