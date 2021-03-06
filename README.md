rdtables
========

Generate tables of values of disc scale length for use in galaxy formation 
models. This program is based in the disc formation model by Mo, Mao & 
White (1998).

RUNNING THE PROGRAM
===================
This program is written in C and has been tested with gcc version 4.2.4. To
compile, just run 'make'.

To run the program:

   ./rdtables parameterfile

The program generates tables for selected values of redshift, halo
concentration, and bulge mass fraction m_b = M_bulge / M_vir. The virial
mass is defined as the mass enclosed by r_200, the radius within which the
enclosed mass density is 200 times the universal critical density.

Two different sets of tables can created: rdtable_NNN_NNN stores values of
disc scale length, and vctable_NNN_NNN stores the value of the circular
velocity at the flat portion of the rotation curve (here identified with
the circular velocity at a specified multiple of the scale radius). Both
sets of data are always generated as part of the calculation, but you can
turn off the generation of one or the other set of files to save time.

We have chosen to use the velocity at the flat portion of the rotation
curve because bulges are modeled as point masses, and thus the maximum
circular velocity is not a good choice. Another option would be to select
the velocity at the optical radius, defined e.g. as the radius which
encloses a given percentage of the total mass, but this also causes
problems whenever the bulge mass alone is equal to or greater than the
selected fraction.

Compilation options (enable in the makefile):
- TESTMODE runs a single calculation in interactive mode and exits, without
  generating any table files. This mode requires the user to provide a set
  of tables that relate halo concentration to virial mass for the selected
  cosmological model (generate these tables with program conctables).
- DEBUG prints extra information while running.

Each table file generated by this program is a collection of double-entry
tables for a given redshift and halo concentration: table_m(i,j) is the value
corresponding to row i in disc mass fraction m_d = M_disc / M_vir and
column j in halo spin parameter. There are Nbulge such tables in each file,
one for each of the Nbulge considered values of bulge mass fraction. The
format of the files is the following:

(Redshift) (Concentration)
(Bulge mass fraction[0])
(Value[0][0] Value[0][1] ... Value[0][Ncols-2] Value[0][Ncols-1]
 Value[1][0] Value[1][1] ... Value[1][Ncols-2] Value[1][Ncols-1]
 .
 .
 .
 Value[Nrows-1][0] Value[Nrows-1][1] ... Value[Nrows-1][Ncols-1])
(Bulge mass fraction[1])
(Value[0][0] Value[0][1] ... Value[0][Ncols-2] Value[0][Ncols-1]
 Value[1][0] Value[1][1] ... Value[1][Ncols-2] Value[1][Ncols-1]
 .
 .
 .
 Value[Nrows-1][0] Value[Nrows-1][1] ... Value[Nrows-1][Ncols-1])
.
.
.

An additional file, tables_header.dat, is created in the output directory.
This file lists all the selected values of redshift, concentration, number
of rows, columns and values of m_b, then label values of rows, columns and
m_b.


READING/TESTING THE OUTPUT
==========================
Directory test/ contains a sample program, rdtest.c, that reads the tables
generated with rdtables and interpolates to find the scale radius and
circular velocity for a given set of values entered by the user. Compile it
by running 'make' inside the test/ directory (it has its own makefile) and
to run it type

   ./rdtest path_to_table_files

You can use this program as a template to include these tables in another
code.

The values obtained from this test program should not differ from those
obtained by running rdtables in test mode by more than 10^-3. Please email
me if you find a case where the differences are higher than this.


DISTRIBUTION
===========
This package is to be publicly available from me. Feel free to distribute
the package as it is (unchanged). Please do not pass modified versions to
others, as the resulting variants may cause errors that I can no longer
sort out. Instead of modifying yourself, suggest the modification to me (if
it is a good idea, others can then use it as well). It is presumably always
a good idea to ask me for another copy before distributing one yourself,
since in the meantime I might well have updated the package to a new
version.

