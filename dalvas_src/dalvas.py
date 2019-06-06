import numpy as np
import scipy.linalg
from dalvas_src import dal_fileloader
from dalvas_src import result_printer
from dalvas_src import dalvas_utility
from pyscf import gto


class DALVAS():
    def __init__(self, dalton_output, mo_punched_out):
        self.S_zeroing_threshold = 10**-12
    
        with open(dalton_output, "r", encoding="utf-8") as f:
            self.__dalton_output_file = list(f)
        with open(mo_punched_out, "r", encoding="utf-8") as f:
            self.__mo_punched_out_file = list(f)
        
        self.coordinates = dal_fileloader.load_coordinates(self.__dalton_output_file)
        self.charge = dal_fileloader.load_charge(self.__dalton_output_file)
        self.number_orbitals = dal_fileloader.load_number_orbitals(self.__dalton_output_file)
        self.number_basisfunctions = dal_fileloader.load_number_basis_functions(self.__dalton_output_file)
        self.number_occupied_orbitals = dal_fileloader.load_number_occupied_orbitals(self.__dalton_output_file)
        self.number_symmetries = len(self.number_occupied_orbitals)
        self.symmetry = dal_fileloader.load_symmetry(self.__dalton_output_file)
        if self.symmetry != "C1":
            self._symmetry_orbital_definition = dal_fileloader.load_symmetry_orbital_definitions(self.__dalton_output_file)
            self._symmetry_to_ao_ordering = dal_fileloader.load_symorder_to_aoorder(self.__dalton_output_file)
        else:
            self._symmetry_orbital_definition = {1: []}
            self._symmetry_to_ao_ordering = {}
            for i in range(1, self.number_orbitals[0]+1):
                self._symmetry_orbital_definition[1].append([str(i),str(i)])
                self._symmetry_to_ao_ordering[i] = i
        if np.sum(np.array(self.number_orbitals) - np.array(self.number_basisfunctions)) != 0:
            print("Number of basis functions is different from number of orbitals")
            print("Run the calculation with:")
            print("*ORBITAL INPUT")
            print(".AO DELETE")
            print(" 1.0D-20")
            print(".CMOMAX")
            print(" 1.0D+20")
        else:
            self.mo_coefficients = dal_fileloader.load_punchout_mo_coefficients(self.__mo_punched_out_file, self.number_orbitals)
        
        self.atom_to_idx = {}
        for idx,atom in enumerate(self.coordinates):
            if atom[0] not in self.atom_to_idx:
                self.atom_to_idx[atom[0]] = [idx]
            else:
                self.atom_to_idx[atom[0]].append(idx)
           
           
    def AVAS(self):
        self.AVAS_eigenvalues_occupied = {}
        self.AVAS_eigenvectors_occupied = {}
        self.AVAS_eigenvalues_virtuel = {}
        self.AVAS_eigenvectors_virtuel = {}
        self.AVAS_mo_coefficients = {}
        for symmetry in range(1, self.number_symmetries+1):
            occupied_idx = self.number_occupied_orbitals[symmetry-1]
            S_reference_inv = np.linalg.inv(self.S_reference)
            Projection_matrix = np.dot(self.S_mixed[symmetry].T,np.dot(S_reference_inv,self.S_mixed[symmetry]))
            S_occupied = np.dot(self.mo_coefficients[symmetry][:,:occupied_idx].T,np.dot(Projection_matrix,self.mo_coefficients[symmetry][:,:occupied_idx]))
            eig_val, eig_vec = np.linalg.eigh(S_occupied)
            self.AVAS_eigenvalues_occupied[symmetry] = eig_val
            self.AVAS_eigenvectors_occupied[symmetry] = eig_vec
        for symmetry in range(1, self.number_symmetries+1):
            occupied_idx = self.number_occupied_orbitals[symmetry-1]
            S_reference_inv = np.linalg.inv(self.S_reference)
            Projection_matrix = np.dot(self.S_mixed[symmetry].T,np.dot(S_reference_inv,self.S_mixed[symmetry]))
            S_virtuel = np.dot(self.mo_coefficients[symmetry][:,occupied_idx:].T,np.dot(Projection_matrix,self.mo_coefficients[symmetry][:,occupied_idx:]))
            eig_val, eig_vec = np.linalg.eigh(S_virtuel)
            self.AVAS_eigenvalues_virtuel[symmetry] = eig_val[::-1]
            self.AVAS_eigenvectors_virtuel[symmetry] = eig_vec[:,::-1]
        for symmetry in range(1, self.number_symmetries+1):
            occupied_idx = self.number_occupied_orbitals[symmetry-1]
            C_occupied = np.dot(self.mo_coefficients[symmetry][:,:occupied_idx],self.AVAS_eigenvectors_occupied[symmetry])
            C_virtuel = np.dot(self.mo_coefficients[symmetry][:,occupied_idx:],self.AVAS_eigenvectors_virtuel[symmetry])
            C_new = np.hstack((C_occupied, C_virtuel))
            eigVal, eigVec = np.linalg.eigh(np.dot(C_new.T,np.dot(self.S_calculation[symmetry],C_new)))
            diaVal = np.diag(eigVal)
            if len(diaVal) > 0:
                M = np.dot(np.dot(eigVec,scipy.linalg.sqrtm(np.linalg.inv(np.diag(eigVal)))),np.matrix.transpose(eigVec))
                C_orthonormal = np.dot(C_new,M)
                self.AVAS_mo_coefficients[symmetry] = C_orthonormal
            else:
                self.AVAS_mo_coefficients[symmetry] = False
    
    
    def punch_mo_coefficients(self, output_name="MO_PUNCHIN"):
        result_printer.punchin_mo_coefficients(self.AVAS_mo_coefficients,outfile_name=output_name)
        
        
    def print_AVAS_eigenvalues(self, print_eigenvalue_threshold=10**-3, cas_eigenvalue_threshold=0.01):
        result_printer.format_AVAS_eigenvalues(self.AVAS_eigenvalues_occupied, self.AVAS_eigenvalues_virtuel, print_threshold=print_eigenvalue_threshold)
        print("")
        print("Partial Dalton input; with eigenvalue threshold of: {:02.3f}".format(cas_eigenvalue_threshold))
        active_occupied = np.zeros(self.number_symmetries)
        active_virtuel = np.zeros(self.number_symmetries)
        for key in self.AVAS_eigenvalues_occupied:
            active_occupied[key-1] = np.sum(self.AVAS_eigenvalues_occupied[key]>cas_eigenvalue_threshold)
            active_virtuel[key-1] = np.sum(self.AVAS_eigenvalues_virtuel[key]>cas_eigenvalue_threshold)
        print(".ELECTRONS")
        print(" "+str(int(np.sum(active_occupied))*2))
        print(".INACTIVE")
        print(" ",end="")
        for i in range(0, len(active_occupied)):
            print(str(int(self.number_occupied_orbitals[i]-active_occupied[i])),end=" ")
        print("")
        print(".CAS SPACE")
        print(" ",end="")
        for i in range(0, len(active_occupied)):
            print(str(int(active_occupied[i]+active_virtuel[i])),end=" ")
        print("")
    
    def set_atomic_valence_selection(self, label_dict, calculation_basis, reference_basis="minao"):
        atom_string = ""
        basis_dict = {}
        for atom in self.coordinates:
            atom_string += atom[0] + " " + str(atom[1]) + " " + str(atom[2]) + " " + str(atom[3]) + "; "
            if atom[0] not in basis_dict:
                basis_dict[atom[0]] = calculation_basis
        atom_string_without_targets = atom_string
        for key in label_dict:
            if isinstance(key, int):
                atom_string += self.coordinates[key][0] + "00 " + str(self.coordinates[key][1]) + " " + str(self.coordinates[key][2]) + " " + str(self.coordinates[key][3]) + "; "
                if self.coordinates[key][0] + "00" not in basis_dict:
                    basis_dict[self.coordinates[key][0] + "00"] = reference_basis
            else:
                indices = self.atom_to_idx[key]
                for idx in indices:
                    atom_string += self.coordinates[idx][0] + "00 " + str(self.coordinates[idx][1]) + " " + str(self.coordinates[idx][2]) + " " + str(self.coordinates[idx][3]) + "; "
                    if self.coordinates[idx][0] + "00" not in basis_dict:
                        basis_dict[self.coordinates[idx][0] + "00"] = reference_basis
        mol = gto.Mole(atom=atom_string, unit="bohr")
        mol.basis = basis_dict
        mol.charge = self.charge
        try:
            mol.build()
        except:
            mol.charge = self.charge - 1
            mol.build()
        S_total = mol.intor('int1e_ovlp_sph')
        S_total[np.abs(S_total)<self.S_zeroing_threshold] = 0
        S_calculation_temp = S_total[:np.sum(self.number_orbitals),:np.sum(self.number_orbitals)]
        self.S_reference = S_total[np.sum(self.number_orbitals):,np.sum(self.number_orbitals):]
        S_mixed_temp = S_total[:np.sum(self.number_orbitals),np.sum(self.number_orbitals):]
        picked_array = np.full(len(S_mixed_temp[0]), False)
        idx_shift = np.sum(self.number_orbitals) # due to only wanting the idx of dummy placed atoms.
        for key in label_dict:
            if isinstance(key, int):
                atom = self.coordinates[key][0] + "00"
                for orbital in label_dict[key]:
                    for idx,label in enumerate(mol.ao_labels()):
                        if atom in label and orbital in label:
                            picked_array[idx-idx_shift] = True
            else:
                atom = self.coordinates[self.atom_to_idx[key][0]][0] + "00"
                for orbital in label_dict[key]:
                    for idx,label in enumerate(mol.ao_labels()):
                        if atom in label and orbital in label:
                            picked_array[idx-idx_shift] = True
        S_mixed_temp = S_mixed_temp.T[picked_array].T
        self.S_reference = self.S_reference[picked_array]
        self.S_reference = self.S_reference.T[picked_array].T
        self.S_calculation = {}
        self.S_mixed = {}
        for idx,orbitals in enumerate(self.number_orbitals, 1):
            self.S_calculation[idx] = np.zeros((orbitals,orbitals))
            self.S_mixed[idx] = np.zeros((orbitals,len(S_mixed_temp[0])))   
        for key in self._symmetry_orbital_definition:
            orbital_definition = self._symmetry_orbital_definition[key]
            for idx_i in range(len(orbital_definition)):
                for idx_j in range(len(orbital_definition)):
                    symmetry_idx_shift = int(orbital_definition[0][0]) - 1 # due to the way indices are printed from Dalton
                    symmetry_idx_i = int(orbital_definition[idx_i][0]) - 1 - symmetry_idx_shift
                    nosymmetry_idx_i = orbital_definition[idx_i][1]
                    symmetry_idx_j = int(orbital_definition[idx_j][0]) - 1 - symmetry_idx_shift
                    nosymmetry_idx_j = orbital_definition[idx_j][1]
                    if "+" in nosymmetry_idx_i:
                        ai = self._symmetry_to_ao_ordering[int(nosymmetry_idx_i.split("+")[0])] - 1
                        bi = self._symmetry_to_ao_ordering[int(nosymmetry_idx_i.split("+")[1])] - 1
                        if "+" in nosymmetry_idx_j:
                            aj = self._symmetry_to_ao_ordering[int(nosymmetry_idx_j.split("+")[0])] - 1
                            bj = self._symmetry_to_ao_ordering[int(nosymmetry_idx_j.split("+")[1])] - 1
                            self.S_calculation[key][symmetry_idx_i,symmetry_idx_j] = S_calculation_temp[ai,aj] + S_calculation_temp[ai,bj] + S_calculation_temp[bi,aj] + S_calculation_temp[bi,bj]
                        elif "-" in nosymmetry_idx_j:
                            aj = self._symmetry_to_ao_ordering[int(nosymmetry_idx_j.split("-")[0])] - 1
                            bj = self._symmetry_to_ao_ordering[int(nosymmetry_idx_j.split("-")[1])] - 1
                            self.S_calculation[key][symmetry_idx_i,symmetry_idx_j] = S_calculation_temp[ai,aj] - S_calculation_temp[ai,bj] + S_calculation_temp[bi,aj] - S_calculation_temp[bi,bj]
                        else:
                            aj = self._symmetry_to_ao_ordering[int(nosymmetry_idx_j)] - 1
                            self.S_calculation[key][symmetry_idx_i,symmetry_idx_j] = S_calculation_temp[ai,aj] + S_calculation_temp[bi,aj]
                    elif "-" in nosymmetry_idx_i:
                        ai = self._symmetry_to_ao_ordering[int(nosymmetry_idx_i.split("-")[0])] - 1
                        bi = self._symmetry_to_ao_ordering[int(nosymmetry_idx_i.split("-")[1])] - 1
                        if "+" in nosymmetry_idx_j:
                            aj = self._symmetry_to_ao_ordering[int(nosymmetry_idx_j.split("+")[0])] - 1
                            bj = self._symmetry_to_ao_ordering[int(nosymmetry_idx_j.split("+")[1])] - 1
                            self.S_calculation[key][symmetry_idx_i,symmetry_idx_j] = S_calculation_temp[ai,aj] + S_calculation_temp[ai,bj] - S_calculation_temp[bi,aj] - S_calculation_temp[bi,bj]
                        elif "-" in nosymmetry_idx_j:
                            aj = self._symmetry_to_ao_ordering[int(nosymmetry_idx_j.split("-")[0])] - 1
                            bj = self._symmetry_to_ao_ordering[int(nosymmetry_idx_j.split("-")[1])] - 1
                            self.S_calculation[key][symmetry_idx_i,symmetry_idx_j] = S_calculation_temp[ai,aj] - S_calculation_temp[ai,bj] - S_calculation_temp[bi,aj] + S_calculation_temp[bi,bj]
                        else:
                            aj = self._symmetry_to_ao_ordering[int(nosymmetry_idx_j)] - 1
                            self.S_calculation[key][symmetry_idx_i,symmetry_idx_j] = S_calculation_temp[ai,aj] - S_calculation_temp[bi,aj]
                    else:
                        ai = self._symmetry_to_ao_ordering[int(nosymmetry_idx_i)] - 1
                        if "+" in nosymmetry_idx_j:
                            aj = self._symmetry_to_ao_ordering[int(nosymmetry_idx_j.split("+")[0])] - 1
                            bj = self._symmetry_to_ao_ordering[int(nosymmetry_idx_j.split("+")[1])] - 1
                            self.S_calculation[key][symmetry_idx_i,symmetry_idx_j] = S_calculation_temp[ai,aj] + S_calculation_temp[ai,bj]
                        elif "-" in nosymmetry_idx_j:
                            aj = self._symmetry_to_ao_ordering[int(nosymmetry_idx_j.split("-")[0])] - 1
                            bj = self._symmetry_to_ao_ordering[int(nosymmetry_idx_j.split("-")[1])] - 1
                            self.S_calculation[key][symmetry_idx_i,symmetry_idx_j] = S_calculation_temp[ai,aj] - S_calculation_temp[ai,bj]
                        else:
                            aj = self._symmetry_to_ao_ordering[int(nosymmetry_idx_j)] - 1
                            self.S_calculation[key][symmetry_idx_i,symmetry_idx_j] = S_calculation_temp[ai,aj]
        for key in self._symmetry_orbital_definition:
            orbital_definition = self._symmetry_orbital_definition[key]
            for idx_i in range(len(orbital_definition)):
                for idx_j in range(len(S_mixed_temp[0])):
                    symmetry_idx_shift = int(orbital_definition[0][0]) - 1 # due to the way indices are printed from Dalton
                    symmetry_idx_i = int(orbital_definition[idx_i][0]) - 1 - symmetry_idx_shift
                    nosymmetry_idx_i = orbital_definition[idx_i][1]
                    if "+" in nosymmetry_idx_i:
                        ai = self._symmetry_to_ao_ordering[int(nosymmetry_idx_i.split("+")[0])] - 1
                        bi = self._symmetry_to_ao_ordering[int(nosymmetry_idx_i.split("+")[1])] - 1
                        self.S_mixed[key][symmetry_idx_i,idx_j] = S_mixed_temp[ai,idx_j] + S_mixed_temp[bi,idx_j]
                    elif "-" in nosymmetry_idx_i:
                        ai = self._symmetry_to_ao_ordering[int(nosymmetry_idx_i.split("-")[0])] - 1
                        bi = self._symmetry_to_ao_ordering[int(nosymmetry_idx_i.split("-")[1])] - 1
                        self.S_mixed[key][symmetry_idx_i,idx_j] = S_mixed_temp[ai,idx_j] - S_mixed_temp[bi,idx_j]
                    else:
                        ai = self._symmetry_to_ao_ordering[int(nosymmetry_idx_i)] - 1
                        self.S_mixed[key][symmetry_idx_i,idx_j] = S_mixed_temp[ai,idx_j]
        for key in self.S_mixed:
            self.S_mixed[key] = self.S_mixed[key].T
            
        # Check sanity
        Input_C_diag = dalvas_utility.check_C_diagnolize_S(self.S_calculation,self.mo_coefficients)
        print("Max deviation from input coefficient diagalization of Overlap:",Input_C_diag)
        if Input_C_diag > 10**-5:
            print(atom_string_without_targets)
            print("Input MO coefficients does not diagonalize the overlap matrix.")
            mol2 = gto.Mole(atom=atom_string_without_targets, unit="bohr")
            mol2.basis = basis_dict
            mol2.charge = self.charge
            try:
                mol2.build()
            except:
                mol2.charge = self.charge - 1
                mol2.build()
            if mol2.nao != np.sum(self.number_basisfunctions):
                print("Number of basis functions in PySCF different from the number of basis functions in DALTON")
                print("PySCF:",mol2.nao,"basis functions")
                print("DALTON:",self.number_basisfunctions,"->",np.sum(self.number_basisfunctions),"basis functions")
                print("PySCF and DALTON is not using the same basis set")
            else:
                print("PySCF and DALTON might have different basis set defintions for the chosen basis set")
            print("")
            
        