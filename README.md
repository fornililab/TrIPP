# TrIPP: Trajectory Iterative pKa Predictor


TrIPP is a Python tool based on PROPKA to calculate the pKa values during Molecular Dynamics trajectories. 

## Prerequisites

This project requires Python (version 3.9 or later). To make sure you have the right version available on your machine, try running the following command. 

```sh
$ python --version
Python 3.9
```

If using the visualization class in TrIPP, we require user to have pre-installed PyMOL (https://www.pymol.org/) in your machine.

Note down the path to your PyMOL executable and then provide it in the function.

## Table of contents

- [Project Name](#project-name)
  - [Prerequisites](#prerequisites)
  - [Table of contents](#table-of-contents)
  - [Installation](#installation)
  - [Workflow](#workflow)
  - [Authors](#authors)
  - [License](#license)

## Installation
The poetry installation requires one to install poetry in their system.

Please follow https://python-poetry.org/docs/#installing-with-the-official-installer to install it if you do not already have one.
```sh
poetry env use PATH/TO/python3.9
poetry shell
poetry install 
```

Note that you may experience the installation of one/some get stuck in 'depending' mode.

This is related to the keyring issue as mentioned in this thread https://github.com/python-poetry/poetry/issues/8623.

A workaround of this is to unlock it with the following code: 
```sh
export PYTHON_KEYRING_BACKEND=keyring.backends.null.Keyring
poetry install
```

Please test the installation by:
```
poetry shell
cd tests
python test_installation.py
```
You will need to provide the path to the PyMOL executable when testing the Visualization class.

Mac: /Applications/PyMOL.app/Contents/MacOS/MacPyMOL

Linux: which pymol

Windows: where pymol

Once all three tests are passed, the installation is complete.

### Workflow

Please start the poetry shell for TrIPP then follow [tripp_tutorial](tutorial/tripp_tutorial.ipynb) for a comprehensive workflow.

Run time for the following cells with Macbook Pro (M2 Pro) using 12 cores are as follow:

Number of trajectories: 2

Length of trajectories: 308.7 ns (3087 frames)

Number of atoms: 1960

Total: 6 minutes 12 seconds

Cell 1.1: 1 minute 13 seconds (2 in total to process)

Cell 2.1+2.2: 1 second

Cell 3.1: 22 seconds

Cell 4.1: 4 minutes 47 seconds (4 mutations each trajectory, 8 in total to process)

Cell 5.1: 1 second

Cell 5.2: 1 second
### Authors

* **Christos Matsingos** - [chmatsingos](https://github.com/chmatsingos)
* **Ka Fu Man** [mkf30](https://github.com/mkf30)
* **Arianna Fornili** [fornililab](https://github.com/fornililab)

### License

The library is licensed by **GPL-3.0**
