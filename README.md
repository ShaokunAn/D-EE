# D-EE: a distributed software for visualizing intrinsic structure of large-scale single-cell data

D-EE is a distributed software for visualization of large-scale single cell data written in the C language, and is based on Portable 
Extensible Toolkit Scientific Computation (PETSc) package and Message Passing Interface (MPI) standard. It employs distributed 
storage and distributed computation technique to achieve memory efficiency and high-performance computation simultaneously. 

In our LSSC-IV super cmputer, which is comprised of two 18-core Intel Xeon Gold CPUs with 192 Gb local memory and are interconnected via a proprietary high performance network, the largest sample size that is feasibly analyzed of D-EE is about 850k. 

Before using D-EE, please install PETSc and MPI in advance. The installation of MPI can be done when install PETSc. The version 
of PETSc used in D-EE is 3.11.4 (https://www.mcs.anl.gov/petsc/download/index.html)
The version of PETSc used in D-EE is 3.11.4 with external package 'elemental'.

### Required
+ MPI i.e., openMPI or IntelMPI
+ PETSc v3.11.4 or higher with external package Elemental
+ R, Python
+ operating systems: Linux
+ wget

### Installation 
 **install MPI**
 
 If MPI wasn't not installed in advance, you can download MPI from https://www.mpich.org/ 
 and install it according to its corresponding guidelines or the README file in the download files.
 
**install PETSc**

The PETSc package used in our software needs external Elemental package. You can install them as follows:


```
mkdir PETSc
cd PETSc
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.11.4.tar.gz
tar -xvf petsc-3.11.4.tar.gz
cd petsc-3.11.4/
# put our configure.py into the current directory and add your MPI directory
'--with-mpi-dir=/home/<USERNAME>/<MPI_INSTALL_DIR>'
python configure.py
# according to the final prompt set PETSC_DIR and PETSC_ARCH
make PETSC_DIR=/home/<USERNAME>/PETSC-3.11.4/petsc-3.11.4 PETSC_ARCH=arch-linux-c-opt all
# check if libraries are working
make PETSC_DIR=/home/<USERNAME>/PETSC-3.11.4/petsc-3.11.4 PETSC_ARCH=arch-linux-c-opt check
# set PETSC_DIR and PETSC_ARCH in your ~/.bash_profile 
export PETSC_DIR=/home/<USERNAME>/PETSC-3.11.4/petsc-3.11.4
export PETSC_ARCH=arch-linux-c-opt
# exit .bash_profile and source it in command
source ~/.bash_profile
```


### use D-EE

Here we use the data in D-EE/testdata/TestData.csv as input to show how to use D-EE software

``` 
mkdir ~/D-EE
cd ~/D-EE
git clone https://github.com/ShaokunAn/D-EE.git \\ download D-EE to your computer
cd D-EE
cd src
# make the input file name as TestData.csv, and the output file name as out_testdata
make main
make run
```

**parameters**

In "make run" command, there exist some parameters in the command, the common format
of the command in makefile is:

```
mpirun -np n ./main  -pc_type lu -pc_factor_mat_solver_type elemental -inputfile InputFile.csv -lambda l -outfile OutputFile
```
and the parameters in the command are:
+ n: the number of process you want to use
+ l: the lambda value in objective function, trades off the weights of attractive and repulsive terms, the defaults is 10
+ InputFile.csv: the name of input file in csv format. It's a matrix where rows represent samples and columns represent features. Note the matrix cannot have 
colnames or rownames.
+ OutputFile: the name of output file in csv format. It's a matrix where rows represents samples and columns represent dimensions. 
The default dimension of low-dimensional space is 2.

And in "make run" of D-TSEE, the common format of the command in makefile is:
```
mpirun -np n ./main -pc_type lu -pc_factor_mat_solver_type elemental -lambda l -beta b -inputfile InputFile.csv -timefile TimeFile.csv -outfile OutputFile
```
Besides the parameters the same as in D-EE, the additional parameters are
+ TimeFile.csv: the time file in csv format containing the time point of each sample. It's a one-column csv file and each row represents the time of one sample.
+ b: beta, the beta value in objective function of D-TSEE, trades off weights of expression-based repulsive force and time-based repulsive force. The defaults is 10.

### License

This package is distributed under the GNU GPLv3 license. 
Please see the https://github.com/ShaokunAn/D-EE/blob/master/License file for the complete version of the license.






 



