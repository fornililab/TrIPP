# TrIPP: Trajectory Iterative pKa Predictor


TrIPP is a Python tool based on PROPKA to calculate the pKa values during Molecular Dynamics trajectories. 

## Prerequisites

This project is tested with Python (version 3.9). To make sure you have the right version available on your machine, try running the following command. 

```sh
$ python --version
Python 3.9
```

The Visualization class in TrIPP requires a working installation of PyMOL (https://www.pymol.org/) on your machine.

Please find the path to your PyMOL executable and provides it in the function.

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
One may want to create a virtual environment before installing the package.
```sh
conda create -n tripp python=3.9
conda activate tripp
pip install tripp
```

(Optional) You may want to test the installation with the following:
```sh
git clone https://github.com/fornililab/TrIPP.git
cd TrIPP/tests/
pytest -s test_Installation.py
```
Note that you will be prompted for the path of PyMOL executable when testing the Visualization class.
You may type `skip` to bypass the Visualization class test.

Mac: /Applications/PyMOL.app/Contents/MacOS/MacPyMOL

Linux: which pymol

Windows: where pymol

### Workflow

Please start the conda environement for TrIPP
```sh
conda activate tripp
```
Then follow [tripp_tutorial](tutorial/tripp_tutorial.ipynb) for a comprehensive workflow.

Running the full tutorial on a Macbook Pro (M2 Pro) using 12 cores requires about 6 minutes (2 trajectories, 3087 frames, 1960 atoms).

### Development

Tests for each function is a work in progress.
Users who modifiy the code should pass all tests inside the [tests](test/) directory.

### Authors

* **Christos Matsingos** - [chmatsingos](https://github.com/chmatsingos)
* **Ka Fu Man** [mkf30](https://github.com/mkf30)
* **Arianna Fornili** [fornililab](https://github.com/fornililab)

### License

The library is licensed by **GPL-3.0**
