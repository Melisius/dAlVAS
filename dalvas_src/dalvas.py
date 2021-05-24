import tarfile

import numpy as np
import scipy.linalg

from dalvas_src import dal_fileloader, dalvas_utility, result_printer


class DALVAS:
    def __init__(
        self, dalton_hf_out: str, dalton_hf_targz: str, dalton_integral_out: str, AOONEINT_file: str
    ) -> None:
        self.coordinates = False
        self.charge = False
        self.number_orbitals = False
        self.number_basisfunctions = False
        self.number_occupied_orbitals = False
        self.number_symmetries = False
        self.symmetry = False
        self.mo_coefficients = False

        with open(dalton_hf_out, "r", encoding="utf-8") as f:
            self.__dalton_output_file = list(f)
        self.coordinates = dal_fileloader.load_coordinates(self.__dalton_output_file)
        self.charge = dal_fileloader.load_charge(self.__dalton_output_file)
        self.number_orbitals = dal_fileloader.load_number_orbitals(self.__dalton_output_file)
        self.number_basisfunctions = dal_fileloader.load_number_basis_functions(self.__dalton_output_file)
        self.number_occupied_orbitals = dal_fileloader.load_number_occupied_orbitals(
            self.__dalton_output_file
        )
        self.number_symmetries = len(self.number_occupied_orbitals)
        self.symmetry = dal_fileloader.load_symmetry(self.__dalton_output_file)

        tar = tarfile.open(dalton_hf_targz, "r:gz")
        for member in tar.getmembers():
            if member.name == "DALTON.MOPUN":
                f = tar.extractfile(member)
                content = f.read().decode("UTF-8")
                break
        self.__mo_punched_out_file = content.split("\n")[:-1]
        self.mo_coefficients = dal_fileloader.load_punchout_mo_coefficients(
            self.__mo_punched_out_file, self.number_basisfunctions
        )

        with open(dalton_integral_out, "r", encoding="utf-8") as f:
            refrence_output = list(f)
        self._ao_labels = dal_fileloader.load_aolabels(refrence_output, self.symmetry)

        self.S_calculation, self.S_reference, self.S_mixed = dal_fileloader.load_overlap_matrix(
            AOONEINT_file, self.number_basisfunctions
        )

    def AVAS(self) -> None:
        self.AVAS_eigenvalues_occupied = {}
        self.AVAS_eigenvectors_occupied = {}
        self.AVAS_eigenvalues_virtuel = {}
        self.AVAS_eigenvectors_virtuel = {}
        self.AVAS_mo_coefficients = {}
        for symmetry in range(1, self.number_symmetries + 1):
            occupied_idx = self.number_occupied_orbitals[symmetry - 1]
            S_reference_inv = np.linalg.inv(self.S_reference[symmetry])
            Projection_matrix = np.dot(
                self.S_mixed[symmetry].T, np.dot(S_reference_inv, self.S_mixed[symmetry])
            )
            self.S_occupied = np.dot(
                self.mo_coefficients[symmetry][:, :occupied_idx].T,
                np.dot(Projection_matrix, self.mo_coefficients[symmetry][:, :occupied_idx]),
            )
            eig_val, eig_vec = np.linalg.eigh(self.S_occupied)
            self.AVAS_eigenvalues_occupied[symmetry] = eig_val
            self.AVAS_eigenvectors_occupied[symmetry] = eig_vec
        for symmetry in range(1, self.number_symmetries + 1):
            occupied_idx = self.number_occupied_orbitals[symmetry - 1]
            S_reference_inv = np.linalg.inv(self.S_reference[symmetry])
            Projection_matrix = np.dot(
                self.S_mixed[symmetry].T, np.dot(S_reference_inv, self.S_mixed[symmetry])
            )
            self.S_virtuel = np.dot(
                self.mo_coefficients[symmetry][:, occupied_idx:].T,
                np.dot(Projection_matrix, self.mo_coefficients[symmetry][:, occupied_idx:]),
            )
            eig_val, eig_vec = np.linalg.eigh(self.S_virtuel)
            self.AVAS_eigenvalues_virtuel[symmetry] = eig_val[::-1]
            self.AVAS_eigenvectors_virtuel[symmetry] = eig_vec[:, ::-1]
        for symmetry in range(1, self.number_symmetries + 1):
            occupied_idx = self.number_occupied_orbitals[symmetry - 1]
            C_occupied = np.dot(
                self.mo_coefficients[symmetry][:, :occupied_idx], self.AVAS_eigenvectors_occupied[symmetry]
            )
            C_virtuel = np.dot(
                self.mo_coefficients[symmetry][:, occupied_idx:], self.AVAS_eigenvectors_virtuel[symmetry]
            )
            C_new = np.hstack((C_occupied, C_virtuel))
            eigVal, eigVec = np.linalg.eigh(np.dot(C_new.T, np.dot(self.S_calculation[symmetry], C_new)))
            diaVal = np.diag(eigVal)
            if len(diaVal) > 0:
                M = np.dot(
                    np.dot(eigVec, scipy.linalg.sqrtm(np.linalg.inv(np.diag(eigVal)))),
                    np.matrix.transpose(eigVec),
                )
                C_orthonormal = np.dot(C_new, M)
                self.AVAS_mo_coefficients[symmetry] = C_orthonormal
            else:
                self.AVAS_mo_coefficients[symmetry] = False

    def punch_mo_coefficients(self, output_name: str = "MO_PUNCHIN") -> None:
        result_printer.punchin_mo_coefficients(self.AVAS_mo_coefficients, outfile_name=output_name)
        
    def return_punchin_mo_coefficients(self) -> str:
        return result_printer.return_punchin_mo_coefficients(self.AVAS_mo_coefficients)

    def print_AVAS_eigenvalues(
        self,
        print_eigenvalue_threshold: float = 10 ** -3,
        space_eigenvalue_threshold: float = 0.1,
        silent: bool = False,
    ) -> None:
        if not silent:
            result_printer.format_AVAS_eigenvalues(
                self.AVAS_eigenvalues_occupied,
                self.AVAS_eigenvalues_virtuel,
                print_threshold=print_eigenvalue_threshold,
            )
        self._active_occupied = np.zeros(self.number_symmetries)
        self._active_virtuel = np.zeros(self.number_symmetries)
        for key in self.AVAS_eigenvalues_occupied:
            self._active_occupied[key - 1] = np.sum(
                self.AVAS_eigenvalues_occupied[key] > space_eigenvalue_threshold
            )
            self._active_virtuel[key - 1] = np.sum(
                self.AVAS_eigenvalues_virtuel[key] > space_eigenvalue_threshold
            )

    def print_CAS(self, eigenvalue_threshold: float = 0.1) -> None:
        self.print_AVAS_eigenvalues(space_eigenvalue_threshold=eigenvalue_threshold, silent=True)
        print("Partial Dalton input; Eigenvalue threshold of: {:02.3f}".format(eigenvalue_threshold))
        print("")
        print(".ELECTRONS")
        print(" " + str(int(np.sum(self._active_occupied)) * 2))
        print(".INACTIVE")
        print("", end=" ")
        for i in range(0, len(self._active_occupied)):
            print(str(int(self.number_occupied_orbitals[i] - self._active_occupied[i])), end=" ")
        print("")
        print(".CAS SPACE")
        print("", end=" ")
        for i in range(0, len(self._active_occupied)):
            print(str(int(self._active_occupied[i] + self._active_virtuel[i])), end=" ")
        print("")

    def print_RAS(self, eigenvalue_threshold: float = 0.1) -> None:
        print("Partial Dalton input; Eigenvalue threshold of: {:02.3f}".format(eigenvalue_threshold))
        print("")
        print(self.return_RAS(eigenvalue_threshold=eigenvalue_threshold))

    def set_atomic_valence_selection(self, label_dict):
        picked_array = {}
        for key in self.S_mixed:
            picked_array[key] = np.full(len(self.S_mixed[key][0]), False)
        for atom in label_dict:
            for symmetry in self._ao_labels:
                idx_shift = self.number_basisfunctions[
                    symmetry - 1
                ]  # due to only wanting the idx of dummy placed atoms.
                for bf_idx in range(self.number_basisfunctions[symmetry - 1], len(self._ao_labels[symmetry])):
                    if (
                        atom in self._ao_labels[symmetry][bf_idx][0]
                        and label_dict[atom] in self._ao_labels[symmetry][bf_idx][1]
                    ):
                        picked_array[symmetry][bf_idx - idx_shift] = True
        for key in picked_array:
            self.S_mixed[key] = self.S_mixed[key].T[picked_array[key]]
            self.S_reference[key] = self.S_reference[key][picked_array[key]]
            self.S_reference[key] = self.S_reference[key].T[picked_array[key]]
        # Check sanity
        Input_C_diag = dalvas_utility.check_C_diagnolize_S(self.S_calculation, self.mo_coefficients)
        if Input_C_diag > 10 ** -4:
            print("Input MO coefficients does not diagonalize the overlap matrix.")
            print("Max deviation from input coefficient diagalization of Overlap:", Input_C_diag)

    def return_RAS(self, eigenvalue_threshold: float = 0.1) -> str:
        self.print_AVAS_eigenvalues(space_eigenvalue_threshold=eigenvalue_threshold, silent=True)
        partial_ras_input = ""
        partial_ras_input += ".ELECTRONS\n"
        partial_ras_input += f" {int(np.sum(self._active_occupied)) * 2}\n"
        partial_ras_input += ".INACTIVE\n"
        for i in range(0, len(self._active_occupied)):
            partial_ras_input += f" {int(self.number_occupied_orbitals[i] - self._active_occupied[i])}"
        partial_ras_input += "\n"
        partial_ras_input += ".RAS1 SPACE\n"
        for i in range(0, len(self._active_occupied)):
            partial_ras_input += f" {int(self._active_occupied[i])}"
        partial_ras_input += "\n"
        partial_ras_input += ".RAS2 SPACE\n"
        for i in range(0, len(self._active_occupied)):
           partial_ras_input += " 0"
        partial_ras_input += "\n"
        partial_ras_input += ".RAS3 SPACE\n"
        for i in range(0, len(self._active_occupied)):
            partial_ras_input += f" {int(self._active_virtuel[i])}"
        partial_ras_input += "\n"
        return partial_ras_input