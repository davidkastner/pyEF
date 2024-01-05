![Graphical Summary of README](docs/_static/header.webp)
==============================
# pyEF
*Python electric fields package*

[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/davidkastner/pyef/workflows/CI/badge.svg)](https://github.com/davidkastner/pyef/actions?query=workflow%3ACI)
[![Documentation Status](https://readthedocs.org/projects/pyef/badge/?version=latest)](https://pyef.readthedocs.io/en/latest/?badge=latest)

## Table of Contents
1. **Overview**
2. **Tutorials**
3. **Installation**
    * Download the package from GitHub
    * Creating a python environment
    * Developer install of E3Bind
    * Supporting installations
4. **What is included?**
    * File structure
    * Command Line Interface
5. **Documentation**
    * Update the ReadTheDocs
    * GitHub refresher


## 1. Overview
The purpose of pyEF is to make electric field and electrostatics calculations more accessible.
pyEF is a python package optimized to run using molden files from QM calculations.

## 2. Tutorials
In progress

## 3. Installation
Install the package by running the follow commands inside the downloaded repository. 
This will perform a developmental version install. 
It is good practice to do this inside of a virtual environment.

### Download the package from GitHub
```
git clone git@github.com:davidkastner/pyEF.git
```

### Creating python environment
All the dependencies can be loaded together using the prebuilt environment.yml file.
Compatibility is automatically tested for python versions 3.8 and higher.
Installing all dependencies together has shown to produce more robust installations.

```
cd pyEF
conda env create -f environment.yml
conda activate pyEF
```

### Developer install of pyQMMM
```
cd pyEF
python -m pip install -e .
```

## 4. What's included?
pyEF contains a collection of optimized scripts that work together in pre-built workflows.
If a script is not already included for procedure, many of the functions may be useful in building a procedure.

### File structure

```
.
├── docs                # Readthedocs documentation site
├── pyEF                # Directory containing E3Bind modules
│   ├── cli.py          # Entry point for running E3Bind tools and workflows
│   ├── run             # Runs the workflow
│   ├── analysis        # Analyze the data
│   └── geometry        # Check geometries
└── ...
```


## 5. Documentation
Accurate documentation will always be a high priority for the project.
You can find documentation at the project's [ReadtheDocs](https://pyEF.readthedocs.io/).

### Update the ReadTheDocs

```
make clean
make html
```

### GitHub refresher
#### Push new changes

```
git status
git pull
git add -A .
git commit -m "Change a specific functionality"
git push -u origin main
```

#### Making a pull request
```
git checkout main
git pull

# Before you begin making changes, create a new branch
git checkout -b new-feature-branch
git add -A
git commit -m "Detailed commit message describing the changes"
git push -u origin new-feature-branch

# Visit github.com to add description, submit, merge the pull request

# Once finished on github.com, return to local
git checkout main
git pull

# Delete the remote branch
git branch -d new-feature-branch
```

#### Handle merge conflict

```
git stash push --include-untracked
git stash drop
git pull
```

#### Acknowledgements
Authors: Melissa Manetsch and David W. Kastner

[MolSSi template](https://github.com/molssi/cookiecutter-cms) version 1.6.