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


---------------------------rdf_partial.py-----------------------------------------------------------------
This code computes the full rdf and partial rdfs for selected pairs of atoms. The purpose of the code is to 
compute partial rdfs for systems with large number of atom types. The input file provides all the further 
necessary information.

1: Enetr cut-off distance in LAMMPS units:

4.800

2: Enter bin size for RDF calculations in the LAMMPS units:

0.020

3:Enter total number of atoms

10479

4: Enter number of atomic types (from LAMMPS data file)

8

5: Eneter atomic masses for each atomic type in LAMMPS file (i.e., the number
of entries must be equal to the above number of types:

16.000 12.010 1.008 14.010 32.060 16.000 19.000 24.305

6: Enter the number of partial RDFs to be computed:

6

7: Enter atomic types for which partial RDF should be
computed as a matrix consisting of two rows:

8  8  8  4  8  5
1  2  5  7  7  6

8: Enter "multi" and a unique (for the set) part of LAMMPS coordinate file (dump file) names:
or "single" and a name of your file with multiple frames

single file_name

9: Select your working directory:

/Users/user_name/Folder1/Folder2/


Coordination number uses the following input file (i_coord.dat)

1:[RDF] Enetr cut-off distance in LAMMPS units:
4.800
2:[RDF] Enter bin size for RDF calculations in the LAMMPS units:
0.020
3:Enter total number of atoms
10479
4:Enter number of atomic types (from LAMMPS data file)
8 
5:Eneter atomic masses for each atomic type in LAMMPS file (i.e., the number 
of entries must be equal to the above number of types:
16.000 12.010 1.008 14.010 32.060 16.000 19.000 24.305 
6:Enter the number of partial RDFs to be computed:
6
7:Enter atomic types for which partial RDF should be
computed as a matrix consisting of two rows:
8  8  8  4  8  5
1  2  5  7  7  6
8:Enter "multi" and a unique (for the set) part of LAMMPS coordinate file (dump file) names:
or "single" and a name of your file with multiple frames 
single Mg_2TFSI_G1.lammpstrj
9:Select your working directory:
/Users/username/FolderOne/Folder2/Folder3/
