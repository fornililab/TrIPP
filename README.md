# TrIPP: Trajectory Iterative pKa Predictor


TrIPP (Trajectory Iterative pKa Predictor) is a Python tool to monitor and analyse changes in the pKa of ionisable residues during Molecular Dynamics simulations of proteins.

TrIPP uses [PROPKA 3](https://github.com/jensengroup/propka), which is licensed under
the GNU Lesser General Public License (LGPL). PROPKA is not included in this repository
but is required to run the software. See the PROPKA repository for its license terms. 

## Prerequisites

This project has been developed and tested with Python (version 3.9). To make sure you have the right version available on your machine, run the following command: 

```sh
$ python --version
Python 3.9
```

The Visualization class in TrIPP requires a working installation of PyMOL (https://www.pymol.org/) on your machine.

## Table of contents

- [Project Name](#project-name)
  - [Prerequisites](#prerequisites)
  - [Table of contents](#table-of-contents)
  - [Installation](#installation)
  - [Workflow](#workflow)
  - [Development](#development)
  - [Authors](#authors)
  - [License](#license)

## Installation
The recommended way to install TrIPP is via PyPI.
You may want to create a virtual environment before installing the package.
```sh
conda create -n tripp python=3.9
conda activate tripp
pip install tripp
```

(Optional) You can test the installation with the following:
```sh
git clone https://github.com/fornililab/TrIPP.git
cd TrIPP/tests/
pytest -s test_Installation.py
```
Note that you will be prompted for the path of the PyMOL executable when testing the Visualization class.
You may type `skip` to bypass the Visualization class test.

Mac: /Applications/PyMOL.app/Contents/MacOS/MacPyMOL

Linux: which pymol

Windows: where pymol

### Workflow

Please start the conda environement for TrIPP
```sh
conda activate tripp
```
then follow the [tripp_tutorial](tutorial/tripp_tutorial.ipynb) for a comprehensive workflow.

Running the full tutorial on a Macbook Pro (M2 Pro) using 12 cores requires about 6 minutes (2 trajectories, 3087 frames, 1960 atoms).

### Development

Tests for each function are a work in progress.
Users who modifiy the code should pass all tests inside the [tests](tests/) directory.

### Authors

* **Christos Matsingos** - [chmatsingos](https://github.com/chmatsingos)
* **Ka Fu Man** [mkf30](https://github.com/mkf30)
* **Arianna Fornili** [fornililab](https://github.com/fornililab)

### Citation

If you publish results produced with TrIPP or develop methods based on TrIPP, please cite the following paper:

* C. Matsingos, K. Man and A. Fornili: *TrIPP: a Trajectory Iterative pKa Predictor*, ...

### License

TrIPP is licensed under **GPL-3.0**
