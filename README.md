# TrIPP: Trajectory Iterative pKa Predictor


TrIPP is a Python tool based on PROPKA to calculate the pKa values during Molecular Dynamics trajectories. 

## Prerequisites

This project requires Python (version 3.7 or later). To make sure you have the right version available on your machine, try running the following command. 

```sh
$ python --version
Python 3.9.1
```

## Table of contents

- [Project Name](#project-name)
  - [Prerequisites](#prerequisites)
  - [Table of contents](#table-of-contents)
  - [Getting Started](#getting-started)
  - [Installation](#installation)
  - [Usage](#usage)
    - [Workflow](#workflow)
  - [Authors](#authors)
  - [License](#license)

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for analysis and development purposes. 

## Installation



## Usage 

### Workflow

Input for main class Trajectory:
- Molecular Dynamics trajectory 
- Topology File
- CPU cores 

Input for run method: 
- name of output file
- extract_surface_data boolean 
- chain identifier 
- mutation integer or list of integers 

Output:
- .csv file containing the pKa values 
- .csv file containing surface data 


### Authors

* **Christos Matsingos** - [chmatsingos](https://github.com/chmatsingos)
* **Arianna Fornili** [fornililab](https://github.com/fornililab)

### License

The library is licensed by **GPL-3.0**
