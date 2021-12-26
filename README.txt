This folder contains a simple python code to compute
manybody XXZ Hamiltonian in the spin basis.
The program uses various symmetry sectors to split the 
specturm of the system into individual blocks and then
diagonalizes in each sector.
The program is heavily based on the notes of Prof. Sandvik
see attached paper.

Sometimes constructing the Hamiltonian takes up a lot 
of time. This module contains methods to construct it
quickly. 

After the Hamiltonian is constructed it is diagonalized
and the respective eigenvalues and eigenvectors are obtained.

Given an initial state, the program evolves the
corrosponding density matrix in time.

Sparse Matrix methods are not used in this implementation.

--2021,Dec 25
Aritra Kundu, Trieste
