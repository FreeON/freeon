---
layout: default
title: Users Guide
---

-   [Input File Syntax](Input File Syntax "wikilink")
-   [Example Inputs](Example Inputs "wikilink")

Downloading and Building the code
---------------------------------

The location of the source code releases are described in [Main\_Page\#Download](Main_Page#Download "wikilink"). Download the latest tar file and unpack it with

    tar xvf freeon*.tar.bz2

Change into the newly created directory (it will be named just like the tar file without the .tar.bz2 suffix), and configure the sources.

    ./configure [options]

The options are described in more detail using the `--help` command line option for configure. The only external library dependencies are BLAS/LAPACK and HDF5. We currently include [netlib/LAPACK](http://www.netlib.org/lapack/) and [HDF5](http://www.hdfgroup.org/HDF5/) 1.8.3 with the tar releases. To use these internal libraries, run

    ./configure --enable-internal-lapack --enable-internal-hdf5

Please note that by linking in a vendor tuned LAPACK library, substantial performance gains in the linear algebra parts of the code can be achieved. See [ATLAS](http://math-atlas.sourceforge.net/), [MKL](http://software.intel.com/en-us/articles/intel-mkl/), or [ACML](http://developer.amd.com/libraries/acml/pages/default.aspx) for instance. Once configured, the sources are built using

    make

and installed with

    make install

[Category:Users Guide](Category:Users Guide "wikilink")
