The code is based on PHYSICAL REVIEW B 97, 045145 (2018)
and SciPost Phys. Lect. Notes 7 (2019), and follows their
notation. One should cite these papers when using this code. 

imps.h and imps.cc give a C++ realization of an infinite matrix product state.  
They are based on iTensor version 3.1.3

The main member functions are: 
void mixed_gauge(): bring an iMPS into the canonical form iteratively 
T get_lw() and T get_rw(): obtain the quasi-fixed points of the MPO transfer matrix
void iTDVP(): carry out the time-dependent variational principle of an iMPS 

See sample for the source code and an example on the dynamics of the Ising model.
