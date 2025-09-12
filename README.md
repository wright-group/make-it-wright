# make-it-wright
Jin Group Tools for handling a variety of common data.
Builds upon the WrightTools Data object.

## Features

- a module for each instrument featured
  - AFM (Gwiddion)
  - Andor Neo Camera (Solis)
  - Becker and Hickl SPCM
  - Horiba LabRAM
  - Generic Images
  - ion TOF 
  - XRD (Bruker)
- preset styles and routines for making quick figures


## Installation

### Basic

First download this repository and unzip it. 
Navigate to the package directory (the folder containing the "pyproject.toml" file). 
Create a conda environment for this package:
```
conda create -n makeitwright pip
conda activate makeitwright
```
Still in this directory, install the package:
```
pip install .
```
If you wish to install in editable mode, use the same steps of above, but include the editable flag when installing:
```
pip install --editable .
```

For whatever IDE you use (pyCharm, Spyder, VSCode), be sure to configure the editor so that you use the correct Conda environment (consult the documentation for the IDE). 

### IonTOF

support for iontof data is optional; if you need to use iontof data, specify additional imports using:

`pip install .[iontof]`

Note that at the time of this writing, iontof must be used on python version <3.13.

## Examples

TODO
