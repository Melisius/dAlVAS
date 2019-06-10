import numpy as np
import scipy.linalg
from dalvas_src import dal_fileloader
from dalvas_src import result_printer
from dalvas_src import dalvas_utility


class DALVAS():
    def __init__(self, dalton_output="None", mo_punched_out="None"):
        self.coordinates = False
        self.charge = False
        self.number_orbitals = False
        self.number_basisfunctions = False
        self.number_occupied_orbitals = False
        self.number_symmetries = False
        self.symmetry = False
        self.mo_coefficients = False
        if dalton_output != "None":
            with open(dalton_output, "r", encoding="utf-8") as f:
                self.__dalton_output_file = list(f)
            self.coordinates = dal_fileloader.load_coordinates(self.__dalton_output_file)
            self.charge = dal_fileloader.load_charge(self.__dalton_output_file)
            self.number_orbitals = dal_fileloader.load_number_orbitals(self.__dalton_output_file)
            self.number_basisfunctions = dal_fileloader.load_number_basis_functions(self.__dalton_output_file)
            self.number_occupied_orbitals = dal_fileloader.load_number_occupied_orbitals(self.__dalton_output_file)
            self.number_symmetries = len(self.number_occupied_orbitals)
            self.symmetry = dal_fileloader.load_symmetry(self.__dalton_output_file)
            
        if mo_punched_out != "None":
            with open(mo_punched_out, "r", encoding="utf-8") as f:
                self.__mo_punched_out_file = list(f)
            self.mo_coefficients = dal_fileloader.load_punchout_mo_coefficients(self.__mo_punched_out_file, self.number_basisfunctions)
        
           
    def AVAS(self):
        self.AVAS_eigenvalues_occupied = {}
        self.AVAS_eigenvectors_occupied = {}
        self.AVAS_eigenvalues_virtuel = {}
        self.AVAS_eigenvectors_virtuel = {}
        self.AVAS_mo_coefficients = {}
        for symmetry in range(1, self.number_symmetries+1):
            occupied_idx = self.number_occupied_orbitals[symmetry-1]
            S_reference_inv = np.linalg.inv(self.S_reference[symmetry])
            Projection_matrix = np.dot(self.S_mixed[symmetry].T,np.dot(S_reference_inv,self.S_mixed[symmetry]))
            S_occupied = np.dot(self.mo_coefficients[symmetry][:,:occupied_idx].T,np.dot(Projection_matrix,self.mo_coefficients[symmetry][:,:occupied_idx]))
            eig_val, eig_vec = np.linalg.eigh(S_occupied)
            self.AVAS_eigenvalues_occupied[symmetry] = eig_val
            self.AVAS_eigenvectors_occupied[symmetry] = eig_vec
        for symmetry in range(1, self.number_symmetries+1):
            occupied_idx = self.number_occupied_orbitals[symmetry-1]
            S_reference_inv = np.linalg.inv(self.S_reference[symmetry])
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
    
    
    def set_atomic_valence_selection(self, reference_dalton_outfile, AOONEINT_file, label_dict):
        self.S_calculation, self.S_reference, self.S_mixed = dal_fileloader.load_overlap_matrix(AOONEINT_file, self.number_basisfunctions)
        with open(reference_dalton_outfile, "r", encoding="utf-8") as f:
            refrence_output = list(f)
        self._ao_labels = dal_fileloader.load_aolabels(refrence_output, self.symmetry)

        picked_array = {}
        for key in self.S_mixed:
            picked_array[key] = np.full(len(self.S_mixed[key][0]), False)
        for atom in label_dict:
            for symmetry in self._ao_labels:
                idx_shift = self.number_basisfunctions[symmetry-1] # due to only wanting the idx of dummy placed atoms.
                for bf_idx in range(self.number_basisfunctions[symmetry-1], len(self._ao_labels[symmetry])):
                    if atom in self._ao_labels[symmetry][bf_idx][0] and label_dict[atom] in self._ao_labels[symmetry][bf_idx][1]:
                        picked_array[symmetry][bf_idx-idx_shift] = True
        for key in picked_array:
            self.S_mixed[key] = self.S_mixed[key].T[picked_array[key]]
            self.S_reference[key] = self.S_reference[key][picked_array[key]]
            self.S_reference[key] = self.S_reference[key].T[picked_array[key]]
        
        # Check sanity
        Input_C_diag = dalvas_utility.check_C_diagnolize_S(self.S_calculation,self.mo_coefficients)
        if Input_C_diag > 10**-4:
            print("Input MO coefficients does not diagonalize the overlap matrix.")
            print("Max deviation from input coefficient diagalization of Overlap:",Input_C_diag)
            
        