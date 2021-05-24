import numpy as np


def punchin_mo_coefficients(CMO, outfile_name):
    outfile = open(outfile_name, "w+")
    outfile.write("**MOLORB")
    for sym in range(1, len(CMO) + 1):
        if isinstance(CMO[sym], bool):
            continue
        number_orbitals = len(CMO[sym])
        for i in range(0, number_orbitals):
            outfile.write("\n")
            for j in range(0, number_orbitals):
                if j != 0 and j % 4 == 0:
                    outfile.write("\n")
                if np.max(np.abs(CMO[sym][:, i])) < 100:
                    outfile.write("{0:18.14f}".format(CMO[sym][j, i]))
                elif np.max(np.abs(CMO[sym][:, i])) < 1000:
                    outfile.write("{0:18.13f}".format(CMO[sym][j, i]))
                elif np.max(np.abs(CMO[sym][:, i])) < 10000:
                    outfile.write("{0:18.12f}".format(CMO[sym][j, i]))
                elif np.max(np.abs(CMO[sym][:, i])) < 100000:
                    outfile.write("{0:18.11f}".format(CMO[sym][j, i]))
                elif np.max(np.abs(CMO[sym][:, i])) < 1000000:
                    outfile.write("{0:18.10f}".format(CMO[sym][j, i]))
                else:
                    print("MO coefficient to large, punchin file will be broken")
    outfile.write("\n")
    outfile.close()


def return_punchin_mo_coefficients(CMO):
    mo_out = ""
    mo_out += "**MOLORB"
    for sym in range(1, len(CMO) + 1):
        if isinstance(CMO[sym], bool):
            continue
        number_orbitals = len(CMO[sym])
        for i in range(0, number_orbitals):
            mo_out += "\n"
            for j in range(0, number_orbitals):
                if j != 0 and j % 4 == 0:
                    mo_out += "\n"
                if np.max(np.abs(CMO[sym][:, i])) < 100:
                    mo_out += "{0:18.14f}".format(CMO[sym][j, i])
                elif np.max(np.abs(CMO[sym][:, i])) < 1000:
                    mo_out += "{0:18.13f}".format(CMO[sym][j, i])
                elif np.max(np.abs(CMO[sym][:, i])) < 10000:
                    mo_out += "{0:18.12f}".format(CMO[sym][j, i])
                elif np.max(np.abs(CMO[sym][:, i])) < 100000:
                    mo_out += "{0:18.11f}".format(CMO[sym][j, i])
                elif np.max(np.abs(CMO[sym][:, i])) < 1000000:
                    mo_out += "{0:18.10f}".format(CMO[sym][j, i])
                else:
                    print("MO coefficient to large, punchin file will be broken")
    mo_out += "\n"
    return mo_out

def format_AVAS_eigenvalues(occ_eigenvalues, virt_eigenvalues, print_threshold):
    occ_eigenvalues_all = np.zeros((len(occ_eigenvalues[1]), 2))
    occ_eigenvalues_all[:, 0] = occ_eigenvalues[1]
    occ_eigenvalues_all[:, 1] = 1
    virt_eigenvalues_all = np.zeros((len(virt_eigenvalues[1]), 2))
    virt_eigenvalues_all[:, 0] = virt_eigenvalues[1]
    virt_eigenvalues_all[:, 1] = 1
    for key in occ_eigenvalues:
        if key == 1:
            continue
        occ_tmp = np.zeros((len(occ_eigenvalues[key]), 2))
        occ_tmp[:, 0] = occ_eigenvalues[key]
        occ_tmp[:, 1] = key
        occ_eigenvalues_all = np.vstack((occ_eigenvalues_all, occ_tmp))
        virt_tmp = np.zeros((len(virt_eigenvalues[key]), 2))
        virt_tmp[:, 0] = virt_eigenvalues[key]
        virt_tmp[:, 1] = key
        virt_eigenvalues_all = np.vstack((virt_eigenvalues_all, virt_tmp))
    idx = np.argsort(occ_eigenvalues_all[:, 0])[::-1]
    occ_eigenvalues_all = occ_eigenvalues_all[idx]
    idx = np.argsort(virt_eigenvalues_all[:, 0])[::-1]
    virt_eigenvalues_all = virt_eigenvalues_all[idx]
    a = 13
    print("Eig.Val.".ljust(a) + "Symm. #")
    for i in range(len(occ_eigenvalues_all)):
        if occ_eigenvalues_all[i, 0] > print_threshold:
            print(
                " {:02.3f}".format(occ_eigenvalues_all[i, 0]).ljust(a)
                + " {:d}".format(int(occ_eigenvalues_all[i, 1]))
            )
        else:
            break
    print("======================")
    for i in range(len(virt_eigenvalues_all)):
        if virt_eigenvalues_all[i, 0] > print_threshold:
            print(
                " {:02.3f}".format(virt_eigenvalues_all[i, 0]).ljust(a)
                + " {:d}".format(int(virt_eigenvalues_all[i, 1]))
            )
        else:
            break
