# Superancillary equations for PC-SAFT

This repository provides runnable, reproducible code for fitting superancillary equations for the PC-SAFT equation of state

By Ian Bell, NIST (ian.bell@nist.gov)

The outputs of running the fitting are included in the repository, but the build files (needed for Visual Studio/ninja/make, etc.) and binaries are not.

## Building and running

On windows, standard cmake commands executed from this a shell in the folder containing this file will do the job:

```
cd bld
cmake ..
cmake --build . --config Release
Release\fitmain
Release\fittest
Release\fitbench
```

On WSL2 in vscode, running on an ubuntu 20.04 image, the C++ build tools will pick up the correct building arguments if you open this folder in the WSL2 image and have installed the C++ build tools in the WSL2 image. Then:

* Select the Build configuration (Debug, Release)
* Select the target (fitmain)
* Build & run
* The outputs will be in the ``output`` folder at completion

On linux, similar to windows:

```
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
./fitmain
./fittest
./fitbench
```

## Python package

To build the python package:
```
python setup.py bdist_wheel
```
to build the binary wheel, which can then be installed with pip

## Plotting

To ensure complete reproducibility, the script used to generate the figures is also included in this repository: ``make_figures.py``

## License

See LICENSE file