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

The options are described in more detail using the `--help` command line option for configure. Once configured, the sources are built using

    make

and installed with

    make install
