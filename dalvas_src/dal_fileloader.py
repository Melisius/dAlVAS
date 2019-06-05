import numpy as np


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
     
     
def load_symmetry_orbital_definitions(load_file):
    found_symmetry_orbitals = False
    current_symmetry = 0
    symmetry_orbital_definition = {}
    for line in load_file:
        if "Symmetry Orbitals" in line:
            found_symmetry_orbitals = True
        elif "Symmetries of electric field" in line or "Symmetry pointer indices" in line:
            break
        elif found_symmetry_orbitals:
            if "Symmetry" in line:
                current_symmetry += 1
                symmetry_orbital_definition[current_symmetry] = []
            elif len(line) > 25 and "orbitals" not in line:
                symmetry_orbital_definition[current_symmetry].append([line[0:6].strip(),line[26:].strip()])
    if found_symmetry_orbitals == False:
        print("Could not find symmetry orbital in Dalton output file")
        symmetry_orbital_definition = False
    return symmetry_orbital_definition
    
    
def load_punchout_mo_coefficients(load_file, number_orbitals):
    mo_coefficients = {}
    for i,orbitals in enumerate(number_orbitals, 1):
        mo_coefficients[i] = np.zeros(orbitals*orbitals)
    idx = 0
    current_symmetry = 1
    for line in load_file:
        if "MOLORB" in line:
            continue
        for i in range(0, len(line)-17, 18):
            mo_coefficients[current_symmetry][idx] = float(line[i:i+18])
            idx += 1
        if number_orbitals[current_symmetry-1]**2 == idx:
            idx = 0
            current_symmetry += 1
    for i,orbitals in enumerate(number_orbitals, 1):
        mo_coefficients[i] = mo_coefficients[i].reshape((orbitals,orbitals)).T
    return mo_coefficients
    
    
def load_symorder_to_aoorder(load_file):
    found_orderings = False
    ordering_converter = {}
    for line in load_file:
        if "i ->" in line:
            found_orderings = True
        elif "ICNTAO" in line or "Starting in" in line:
            break
        elif found_orderings:
            if len(line) > 22 and "----" not in line:
                orderings = line.split()
                ordering_converter[int(orderings[0])] = int(orderings[1])
    if found_orderings == False:
        print("Could not find orbital orderings in Dalton output file")
        print("Have the calculation been run with:")
        print("*MOLBAS")
        print(".PRINT")
        print(" 10")
        ordering_converter = False
    return ordering_converter
