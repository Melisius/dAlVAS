from dalvas_src import dalvas


def test_punchin_punchout():
    # Load some orbitals in and punch them out unmodified.
    A = dalvas.DALVAS("data/testfiles/hf_ScF6_aug-cc-pVTZ-J.out","data/testfiles/hf_ScF6_aug-cc-pVTZ-J.MOPUN")
    A.AVAS_mo_coefficients = A.mo_coefficients
    A.punch_mo_coefficients("test_punchin_punchout_punchin.MOPUN")
    
    with open("data/testfiles/hf_ScF6_aug-cc-pVTZ-J.MOPUN", "r", encoding="utf-8") as f:
        reference_file = list(f)
    with open("data/test_punchin_punchout_punchin.MOPUN", "r", encoding="utf-8") as f:
        punchin_file = list(f)
    
    for i in range(1, len(reference_file):
        assert reference_file[i] == punchin_file[i]