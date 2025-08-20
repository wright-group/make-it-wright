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

`pip install git+https://github.com/wright-group/make-it-wright.git`

If you use conda, consider making an environment for this package prior to installing:
`conda create -n makeitwright pip`

### IonTOF

support for iontof data is optional; if you need to use iontof data, specify additional imports using:

`pip install git+https://github.com/wright-group/make-it-wright.git[iontof]`

Note that at the time of this writing, iontof must be used on python version <3.13.

## Examples

TODO
