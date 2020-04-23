How to install {#install}
==============

Installation requires at least **CMake** of version 3.2.3. To build STARS-H,
follow these instructions:

1.  Get STARS-H from git repository

        git clone git@github.com:ecrc/stars-h

    or

        git clone https://github.com/ecrc/stars-h

2.  Go into STARS-H folder

        cd stars-h

3.  Get submodules

        git submodule update --init

4.  Create build directory and go there

        mkdir build && cd build

5.  Use CMake to get all the dependencies. For a list of CMake options,
    click [here](@ref cmake_opts)

        cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install/

6.  Build STARS-H

        make -j

7.  Run tests (optional)

        make test

8.  Build local documentation (optional)

        make docs

9.  Install STARS-H

        make install

10. Add line

        export PKG_CONFIG_PATH=/path/to/install/lib/pkgconfig:$PKG_CONFIG_PATH

    to your .bashrc file.

Now you can use pkg-config executable to collect compiler and linker flags for
STARS-H.


