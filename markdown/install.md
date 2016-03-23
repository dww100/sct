---
author:
- name: David W. Wright
  orcid: 0000-0002-5124-8044
  affiliation: Structural Immunology Group, University College London, London
---

# Installing SCT

The SCT suite contains both the modern python based and classic fortran based tools for constrained modelling.
Both are installed at the same time using cmake.

## Requirements

The following software is needed to install SCT:

*  Linux
*  Fortran compiler (gfortran verison 4.3 or newer have been tested)
*  cmake
*  Python 2.7 - with [NumPy](http://www.numpy.org/), [SciPy](http://www.scipy.org/scipylib/index.html), [matplotlib](http://matplotlib.org/) and [PyYAML](http://pyyaml.org/)

[Gfortran](http://gcc.gnu.org/wiki/GFortran) is the name of the GNU Fortran compiler and is available as part of the GNU Compiler Collection (GCC) and [cmake](http://www.cmake.org/) is a cross-platform, open-source build system.
Both are available from the package repositories of most major Linux distributions (i.e. can be installed via apt-get, yum, etc.).

All of the Python packages should also be available via standard package managers.
However, the easiest way to acquire them is to use the free scientific Python distribution Anaconda.
Download and installation instructions for Anaconda can be found [here](https://store.continuum.io/cshop/anaconda/).

## Getting the Source Code

Links to archives of the code are provided at the top of this page in two formats zip and compressed tarball.
To uncompress the zip archive use the command:

    unzip sct.zip

To uncompress the compressed tarball use:

    tar xvfz sct.tar.gz

The code can also be cloned from the GitHub repository using git:

    git clone https://github.com/dww100/sct.git

## Compiling the Code

The compilation and installation is handled by cmake.
To install in user space:

```
cmake -DCMAKE_INSTALL_PREFIX:PATH=/path/to/install . && make install
```

where /path/to/install is the path into which you want the classic SCT executables and modern SCT python scripts to be installed.
NOTE: You must have permission to write to this directory.
Remember to place this directory into your PATH environment variable.

System wide installation is achieved using:

```
sudo cmake -DCMAKE_INSTALL_PREFIX:PATH=/path/to/install . && sudo make install
```
