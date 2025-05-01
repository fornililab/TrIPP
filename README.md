# TrIPP: Trajectory Iterative pKa Predictor


TrIPP is a Python tool based on PROPKA to calculate the pKa values during Molecular Dynamics trajectories. 

## Prerequisites

This project requires Python (version 3.9 or later). To make sure you have the right version available on your machine, try running the following command. 

```sh
$ python --version
Python 3.9
```

The visualization class in TrIPP requires a working installation of PyMOL (https://www.pymol.org/) on your machine.

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

Running the full tutorial on a Macbook Pro (M2 Pro) using 12 cores requires about 6 minutes (2 trajectories, 3087 frames, 1960 atoms).

### Authors

* **Christos Matsingos** - [chmatsingos](https://github.com/chmatsingos)
* **Ka Fu Man** [mkf30](https://github.com/mkf30)
* **Arianna Fornili** [fornililab](https://github.com/fornililab)

### License

The library is licensed by **GPL-3.0**
