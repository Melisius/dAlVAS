
# dAlVAS

AVAS for [Dalton](https://daltonprogram.org/).

AVAS is **A**tomic **V**alence **A**ctive **S**pace.
 Sayfutyarova E. R., Sun, Q., Chan, G. K.-L., & Knizia, G. (2017). [Automated Construction of Molecular Active Spaces from Atomic Valence Orbitals](https://doi.org/10.1021/acs.jctc.7b00128). Journal of Chemical Theory and Computation, 13(9), 4063–4078.

Using AVAS requires to first run two Dalton calculations, here let us use water as an example.

Before being able to use the following description the Dalton code need to be modified (three lines need to be commented out).
In **/DALTON/abacus/herrdn.F** the following three lines needs to be commented out.
```FORTRAN
              IF (.NOT. NOATMD) THEN
               IF (DIST.LT.0.1.AND.CHARGE(M)*CHARGE(N).NE.D0) GOTO 5000
              END IF
```

### Step one; Hartree-Fock
First, a Hartree-Fock calculation using the following dal-file:
```
**DALTON INPUT
.RUN WAVE FUNCTIONS
*MOLBAS
.PRINT
 10
**WAVE FUNCTIONS
.HF
*ORBITAL INPUT
.PUNCHOUTORBITALS
**END OF DALTON INPUT
```
With the following mol-file:
```
BASIS
6-31G**
Title
Description
Atomtypes=2
Charge=8.0 Atoms=1
O      .0000000000        -.2249058930         .0000000000
Charge=1.0 Atoms=2
H     1.4523499293         .8996235720         .0000000000
H    -1.4523499293         .8996235720         .0000000000
```
This can be run as:
```
datlon input.dal molecule.mol
```
### Step two; Integrals
Secondly, a calculation getting the needed overlap integrals is needed, with the following input:
```
**DALTON INPUT
.INTEGRALS
*MOLBAS
.PRINT
 10
**END OF DALTON INPUT
```
Say that we want to target the p-orbitals on oxygen, then the mol-file would look like:
```
ATOMBASIS
Title
Description
Atomtypes=3
Charge=8.0 Atoms=1 Basis=6-31G**
O      .0000000000        -.2249058930         .0000000000
Charge=1.0 Atoms=2 Basis=6-31G**
H     1.4523499293         .8996235720         .0000000000
H    -1.4523499293         .8996235720         .0000000000
Charge=8.0 Atoms=1 Basis=ANO-RCC 2 1
O00      .0000000000        -.2249058930         .0000000000
```
Note here that **Atomtypes** have been increased, and a new oxygen atom has been added on top of the original one. The new atom has been tagged with **00**. The basis set has also been assigned using **ATOMBASIS** instead of **BASIS**.
This can be run as:
```
datlon -get AOONEINT input.dal molecule.mol
```
Note here the additional **-get AOONEINT**

### Step three; dAlVAS
Now that all the required files have been created, we can use AVAS to pick our active space.
```python
from dalvas_src import dalvas

A = dalvas.DALVAS("HARTREEFOCK.out","HARTREEFOCK.tar.gz","INTEGRAL.out","INTEGRAL.AOONEINT")
A.set_atomic_valence_selection({"O00":"p"})
A.AVAS()
```
Names in capital letters are placeholder names.
The eigenvalues can be inspected by running.
```python
A.print_AVAS_eigenvalues()
```
Partial Dalton input for a CAS calculation can be acquired by running.
```python
A.print_CAS()
```
Partial Dalton input for a RAS calculation can be acquired by running.
```python
A.print_RAS()
```
### Step four; New input with AVAS orbitals
The AVAS orbitals can be written to file by running.
```python
A.punch_mo_coefficients("AVASORBITALS.MOPUN")
```
Now an MCSCF or CI calculation can be constructed based on the AVAS orbitals by inserting the following into the file.
```
**DALTON INPUT
.RUN WAVE FUNCTIONS
**WAVE FUNCTIONS
.MCSCF or .CI
*CONFIGURATION INPUT
... !CONFIGURTION INTPUT
*ORBITAL INPUT
.MOSTART
FORM18
.PUNCHINORBITALS
.PUNCHOUTORBITALS
**MOLORB
 -0.00744466432882 -0.04142041072879  0.03712725172283  0.00122516460384
 -0.00209975963736  ... !REST OF MOPUN FILE
 **END OF DALTON INPUT
 ```
Note to **NOT** run a Hartree-Fock before MCSCF or CI in the input file, where the AVAS orbitals are used.

### Step three and a half (OPTIONAL); Inspect AVAS orbitals
The AVAS orbitals can be written to a molden-file by running the following.
```
**DALTON INPUT
.RUN WAVE FUNCTIONS
**WAVE FUNCTIONS
.HF
*SCF INPUT
.NOQCSCF
.NODIIS
.NONCANONICAL
*ORBITAL INPUT
.MOSTART
FORM18
.PUNCHINORBITALS
.PUNCHOUTORBITALS
**MOLORB
 -0.00744466432882 -0.04142041072879  0.03712725172283  0.00122516460384
 -0.00209975963736  ... REST OF MOPUN FILE
**END OF DALTON INPUT
 ```
 The orbitals will then be in *molden.inp* found in the *OUTPUT.tar.gz*
 These orbitals can be viewed in [MOLDEN](https://www3.cmbi.umcn.nl/molden/) or [Jmol](http://jmol.sourceforge.net/).