# Data-Analysis
Atomistic simulations: Data analysis codes

---------------------------cluster-analysis.py-----------------------------------------------------------------

This data processing module finds clusters in a system of molecules (arbitrary number of molecules consisting of 
arbitrary number of atoms) based upon cutoff distance and test each candidate on the basis of force magnitude to 
exclude potentially unstable structures. The input parameters are all included in the input_cluster.dat file. An 
example of input file and explanations are below.

1: Enter the charcateristic distance (LAMMPS units; for the cluster analysis):

2.400
2: Enter the number of atom types [from the LAMMPS data file]:

8

3: Atom type IDs (different from atom IDs) and corresponding chemical element symbol:

1  2  3  4  5  6  7  8

O  C  H  N  S  O  F  Mg

4: Enter the ion type (e.g., Li+, Mg+2, etc.) type ID [1,2,3,...] from the LAMMPS data-file:

8

5: Enter the number of molecule types:

3

6: Enter the number of molecules of each type:

591 66 33

7: Enter the number of atoms in a single molecule of each type:

16 15 1

8: Enter the maximum force for clusterization criterion

0.95

9: Enter the LAMMPS coordinate file (dump file) name:

lammps_dump_coordinate_file

10: Choose Working directory:

/Users/some_name/Folder1/Folder2

11: Enter the number of files and frame number (1 N if single file contains N frames)

1 15

12: Output mode

verbose

