This example gives the dynamics of an infinite transverse-field Ising model (H = -\sum_i Z_i Z_{i+1} + g_x \sum_i X_i) after a quantum quench.
g0 in run.sh is the pre-quenched g_x, and g1 in run.sh is the post-quenched g_x. 

To compile, fix the Makefile so that the itensor3 directory is correct for your machine. 
To run, edit run.sh, and run ./run.sh
To check the result, run compare.sh or plot_itdvp.sh, which use gnuplot. 
The compare.sh compares the time-dependence of the transverse magnetization after the quench with the exact solution. 
(Note that the dynamics is bound to be inaccurate in long time due to the entanglement growth after an extended quench.)
