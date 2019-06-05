from dalvas_src import dalvas
import numpy as np


def test_punchin_punchout():
    # Load some orbitals in and punch them out unmodified.
    A = dalvas.DALVAS("data/testfiles/hf_ScF6_aug-cc-pVTZ-J.out","data/testfiles/hf_ScF6_aug-cc-pVTZ-J.MOPUN")
    A.AVAS_mo_coefficients = A.mo_coefficients
    A.punch_mo_coefficients("data/temp/test_punchin_punchout_punchin.MOPUN")
    
    with open("data/testfiles/hf_ScF6_aug-cc-pVTZ-J.MOPUN", "r", encoding="utf-8") as f:
        reference_file = list(f)
    with open("data/temp/test_punchin_punchout_punchin.MOPUN", "r", encoding="utf-8") as f:
        punchin_file = list(f)
    
    for i in range(1, len(reference_file)):
        assert reference_file[i] == punchin_file[i]
        
        
def test_AVAS_nosymmetry():
    A = dalvas.DALVAS("data/testfiles/input_punchout_FeNOCO3.out","data/testfiles/input_punchout_FeNOCO3.MOPUN")
    A.set_atomic_valence_selection({"Fe":"d"}, "cc-pVDZ-DK")
    A.AVAS()
    
    reference = np.zeros(42)
    reference[37] = 0.5282566095943796
    reference[38] = 0.528647966424764
    reference[39] = 0.8642205015483625
    reference[40] = 0.9443794440304968
    reference[41] = 0.9444059478593652
    for i,eigenvalue in enumerate(A.AVAS_eigenvalues_occupied[1]):
        assert np.abs(eigenvalue - reference[i]) < 10**-6
            
    reference = np.zeros(133)
    reference[0] = 0.4717433416616607
    reference[1] = 0.4713519848343533
    reference[2] = 0.13577944815203144
    reference[3] = 0.055620501241933934
    reference[4] = 0.05559399741335239
    for i,eigenvalue in enumerate(A.AVAS_eigenvalues_virtuel[1]):
        assert np.abs(eigenvalue - reference[i]) < 10**-6
        
        
def test_AVAS_symmetry():
    A = dalvas.DALVAS("data/testfiles/input_punchout_FeNOCO3_sym.out","data/testfiles/input_punchout_FeNOCO3_sym.MOPUN")
    A.set_atomic_valence_selection({"Fe":"d"}, "cc-pVDZ-DK")
    A.AVAS()
    
    reference = np.zeros(29)
    reference[26] = 0.5282551186932969
    reference[27] = 0.864221716045548
    reference[28] = 0.9443789299054083
    for i,eigenvalue in enumerate(A.AVAS_eigenvalues_occupied[1]):
        assert np.abs(eigenvalue - reference[i]) < 10**-6
    reference = np.zeros(13)
    reference[11] = 0.5286494441180509
    reference[12] = 0.9444058557370661
    for i,eigenvalue in enumerate(A.AVAS_eigenvalues_occupied[2]):
        assert np.abs(eigenvalue - reference[i]) < 10**-6
        
    reference = np.zeros(68)
    reference[0] = 0.4717448325627313
    reference[1] = 0.13577823365483857
    reference[2] = 0.05562101536704647
    for i,eigenvalue in enumerate(A.AVAS_eigenvalues_virtuel[1]):
        assert np.abs(eigenvalue - reference[i]) < 10**-6
    reference = np.zeros(45)
    reference[0] = 0.4713505071410456
    reference[1] = 0.055594089535672644
    for i,eigenvalue in enumerate(A.AVAS_eigenvalues_virtuel[2]):
        assert np.abs(eigenvalue - reference[i]) < 10**-6
        