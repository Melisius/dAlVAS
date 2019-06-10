import numpy as np
from daltools import one


def load_coordinates(load_file):
    found_coordinates = False
    Nxyz = []
    for line in load_file:
        if "Total number of coordinates:" in line:
            found_coordinates = True
        elif found_coordinates:
            # Could be shorter than five,
            # just need to fine empty line after coordinates.
            if len(line) < 5:
                break
            N = line[0:5].strip()
            x = float(line[20:36])
            y = float(line[43:59])
            z = float(line[66:82])
            Nxyz.append([N,x,y,z])
    if found_coordinates == False:
        print("Could not find coordinates in Dalton output file")
        Nxyz = False
    return Nxyz
    
    
def load_charge(load_file):
    found_charge = False
    for line in load_file:
        if "Total charge of the molecule" in line:
            charge = int(line.split("molecule")[1])
            found_charge = True
            break
    if found_charge == False:
        print("Could not find charge in Dalton output file")
        charge = False
    return charge
    
    
def load_number_orbitals(load_file):
    found_number_orbitals = False
    for line in load_file:
        if "Total number of orbitals" in line:
            number_orbitals = [int(i) for i in line.split("|")[1].split()]
            found_number_orbitals = True
            break
    if found_number_orbitals == False:
        print("Could not find number of orbitals in Dalton output file")
        number_orbitals = False
    return number_orbitals
    
    
def load_number_basis_functions(load_file):
    found_number_basis_functions = False
    for line in load_file:
        if "Number of basis functions" in line:
            number_basis_functions = [int(i) for i in line.split("|")[1].split()]
            found_number_basis_functions = True
            break
    if number_basis_functions == False:
        print("Could not find number of basis functions in Dalton output file")
        found_number_basis_functions = False
    return number_basis_functions
    
    
def load_number_occupied_orbitals(load_file):
    found_occupations = False
    for line in load_file:
        if "Orbital occupations :" in line:
            orbital_occupations = [int(i) for i in line.split(":")[1].split()]
            found_occupations = True
            break
    if found_occupations == False:
        print("Could not find orbital occupations in Dalton output file")
        orbital_occupations = False
    return orbital_occupations
    
    
def load_symmetry(load_file):
    found_symmetry = False
    for line in load_file:
        if "Total number of symmetries" in line:
            symmetry = line.split(":")[1].split(")")[0].strip()
            found_symmetry = True
            break
    if found_symmetry == False:
        print("Could not find symmetry in Dalton output file")
        symmetry = False
    return symmetry
    
    
def load_punchout_mo_coefficients(load_file, number_basis_functions):
    mo_coefficients = {}
    for i,orbitals in enumerate(number_basis_functions, 1):
        mo_coefficients[i] = np.zeros(orbitals*orbitals)
    idx = 0
    current_symmetry = 1
    for line in load_file:
        if "MOLORB" in line:
            continue
        for i in range(0, len(line)-17, 18):
            mo_coefficients[current_symmetry][idx] = float(line[i:i+18])
            idx += 1
        if number_basis_functions[current_symmetry-1]**2 == idx:
            idx = 0
            current_symmetry += 1
    for i,orbitals in enumerate(number_basis_functions, 1):
        mo_coefficients[i] = mo_coefficients[i].reshape((orbitals,orbitals)).T
    return mo_coefficients
    
    
def load_aolabels(load_file, symmetry):
    found_aolabels = False
    aolabels = {}
    current_symmetry = 0
    if symmetry.lower() == "c1":
        for line in load_file:
            if "Contracted Orbitals" in line:
                found_aolabels = True
                aolabels[1] = []
            elif "Orbital exponents and normalized contraction coefficients" in line:
                break
            elif found_aolabels and len(line) > 20 and "------" not in line:
                atom = line.split()[1]
                label = line.split()[2]
                aolabels[1].append([atom, label])
    else:
        for line in load_file:
            if "Symmetry Orbitals" in line:
                found_aolabels = True
            elif "Symmetries of electric field" in line or "Symmetry pointer indices" in line:
                break
            elif found_aolabels:
                if "Symmetry" in line:
                    current_symmetry += 1
                    aolabels[current_symmetry] = []
                elif len(line) > 25 and "orbitals" not in line:
                    atom = line.split()[1]
                    label = line.split()[2]
                    aolabels[current_symmetry].append([atom, label])
    if found_aolabels == False:
        print("Could not find AO labels in Dalton output file")
        print("Have the calculation been run with:")
        print("*MOLBAS")
        print(".PRINT")
        print(" 10")
        aolabels = False
    return aolabels
    
    
def load_overlap_matrix(AOONEINT_file, number_basis_functions):
    S_total = one.read(label='OVERLAP', filename=AOONEINT_file).unpack()
    S_calculation = {}
    S_reference = {}
    S_mixed = {}
    for i in range(len(number_basis_functions)):
        numb_bf = number_basis_functions[i]
        S_calculation[i+1] = np.array(S_total[i])[:numb_bf,:numb_bf:]
        S_reference[i+1] = np.array(S_total[i])[numb_bf:,numb_bf:]
        S_mixed[i+1] = np.array(S_total[i])[:numb_bf,numb_bf:]
    return S_calculation, S_reference, S_mixed






