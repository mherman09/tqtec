# TQTec

## Overview

**TQTec** and **TQClim** calculate the one-dimensional transient thermal field for an area that experiences geologic, tectonic, and/or climatic events. The package also includes programs for estimating thermochronologic observables corresponding to the temperature history.

Example input model files (which include comments describing their format) are available in the `examples` directory.

### Programs

- **tqtec**: calculate thermal histories for crustal scale models (units: km, Ma)
- **tqclim**: calculate thermal histories for surface scale models (units: m, yr)
- **readtqtec**: post-process **tqtec**/**tqclim** outputs
- **minage**: model mineral ages for temperature histories
- **petro**: model petroleum-related observables for temperature histories


## Scripts

Bash scripts showing how to run models and plot the model results are available in the `scripts` directory. These scripts are included in the `build` and `bin` directories when **TQTec** is built and installed. The plotting scripts require [Generic Mapping Tools
(GMT)](www.generic-mapping-tools.org).


## Installation

### Dependencies

We recommend installing these programs with a package manager, either [Homebrew](https://brew.sh) on macOS or something like APT on Linux OS.

- GNU Fortran and C compilers
- CMake
- Make
- Git
- Generic Mapping Tools (required for plotting scripts)

### Instructions

1. Download **TQTec** from Github. Navigate to the directory where you want to download the **TQTec** code, then type the following commands:

    ```
    git clone https://github.com/mherman09/tqtec ./tqtec
    cd tqtec
    mkdir build
    cd build
    ```

2. From the `build` directory, build **TQTec** with:

    ```
    cmake ..
    cmake --build .
    ```

3. Install the programs in `/opt/tqtec` (you may specify a different installation directory if you prefer) with the following command:

    ```
    cmake --install . --prefix /opt/tqtec
    ```

    You may need to run this command with `sudo` in front of `cmake` to give permission to install in this location.

4. Add `/opt/tqtec/bin` (or wherever you installed **TQTec**) to your PATH, so your computer can find the programs. For example, if your shell is `zsh` and you installed **TQTec** in `/opt/tqtec`, add the following to the `.zprofile` file in your home directory:

    ```
    export PATH=$PATH:/opt/tqtec/bin
    ```


### Updating **TQTec**

Updating **TQTec** is as simple as going to the `build` directory and typing:

```
git pull
```

Then repeat steps 2 and 3 from the installation instructions above.


## Authors

- *Kevin Furlong*: Thermal physicist extraordinaire; created original **TQTec**; advises **TQTec** development
- *Matthew Herman*: Thermal physicist in training; slightly obsessive; updated **TQTec** to Modern Fortran; maintains and develops **TQTec**
- *Chris Guzofski*: Bulk thickening and metamorphism
- *Matthew Legg*: Apatite fission track and petroleum observables
- *Rachel Piotraschke*: (U-Th)/He dating and helium diffusion in apatite
