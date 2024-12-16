<h1 align="center">Atomistic Simulation Instruments and Database</h1>

<p align="center">
	<a href=https://github.com/ASID-Production/ASID/releases/latest><img src=https://img.shields.io/github/v/release/ASID-Production/ASID?sort=date&style=plastic&color=brightgreen></a>
	<img src=https://img.shields.io/badge/C++-14-blue.svg?style=plastic>
    <img src=https://img.shields.io/badge/Python-3.10-blue.svg?style=plastic>
    <img src=https://img.shields.io/badge/Django-3.2.24-blue.svg?style=plastic>
</p>

__ASID__ is an "easy-to-use" open source program designed for generating and statistically processing large numbers of chemical structures.
The project is a database with a graphical shell, specialized for work in chemical research institutes and industry.

----
#### CPPLIB Build Status:

[![Windows Build](https://github.com/ASID-Production/ASID/actions/workflows/cmake-windows.yml/badge.svg)](https://github.com/ASID-Production/ASID/actions/workflows/cmake-windows.yml)
[![Linux Build](https://github.com/ASID-Production/ASID/actions/workflows/cmake-linux.yml/badge.svg)](https://github.com/ASID-Production/ASID/actions/workflows/cmake-linux.yml)

----

<a name="top"></a>
## Table of Contents: 
click [^](#top) to return here
* [System requirements](#SystemRequirements)
* [Installation on Windows](#InstallationW)
* [Installation on Linux](#InstallationL)
* [Licenses](#Licenses)

<a name="SystemRequirements"></a>
## System Requirements <sup>[^](#top)</sup>
* OS: x64 only, Windows 10/11, Linux (Tested on Ubuntu 22+)
* Processor: 2 or more cores, AMD Piledriver / Intel Haswell (2012) or newer are strongly recommended
* RAM: 4 GB
* Graphics Card: Compatible with OpenGL 4.6
* SSD: 1.4 GB

----
<a name="InstallationW"></a>
## Installation on Windows<sup>[^](#top)</sup>
1. Download Windows-based packege `ASID_v1.0.0_win.zip` or `ASID_v1.0.0_win_installer.exe` of [latest stable version](https://github.com/ASID-Production/ASID/releases/latest).
2. Unpack archive or intrall via installer in any sutable dirrectory.
4. Run the program using the `<ASID_ROOT>/VnE/Start.bat`

----
<a name="InstallationL"></a>
## Installation on Linux<sup>[^](#top)</sup>
* Linux (only x64 version). 
* gcc and g++ (version 8.5 or newer). 
* Python (version 3.10 or newer) with python-dev package.
* CMake(version 3.19 or newer).
* make
* Ninja builder

0. Installing the required dependencies (requares super user rights):
   
   `sudo apt-get -y install python python-dev gcc g++ cmake make ninja-build`

1. Download repository

   Download source code of [latest stable version](https://github.com/ASID-Production/ASID/releases/latest) or
   clone git repository using comand line (`git clone https://github.com/ASID-Production/ASID`).
2. Unpack it in any sutable directory.
3. Run `install-linux.sh` in the ASID directory.
4. Run the program using the `<ASID_ROOT>/VnE/Start` file.

----

<a name="Licenses"></a>
## Licenses <sup>[^](#top)</sup>

Licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0) (the "License").
You can use, redistribute it and/or modify it under the terms of the License.

<details><summary>Installed by PyPI dependencies</summary>

|     Used Library      | Version  | License               |
|:---------------------:|:--------:|:---------------------:|
|       `asgiref`       | 3.5.2    | [BSD 3-Clause License](https://github.com/django/asgiref/blob/main/LICENSE) |    
|       `chardet`       | 5.2.0    | [GNU LGPL v2.1](https://github.com/chardet/chardet/blob/main/LICENSE) |
|       `django`        | 3.2.24   | [BSD 3-Clause License](https://github.com/django/django/blob/main/LICENSE) |
|    `django-filter`    | 22.1     | [BSD 3-Clause License](https://github.com/carltongibson/django-filter/blob/main/LICENSE) |
| `djangorestframework` | 3.14.0   | [BSD 3-Clause License](https://github.com/encode/django-rest-framework/blob/master/LICENSE.md) |
|       `djoser`        | 2.1.0    | [MIT License](https://github.com/sunscrapers/djoser/blob/master/LICENSE) |
|      `drf-yasg`       | 1.21.4   | [BSD 3-Clause License](https://github.com/axnsan12/drf-yasg/blob/master/LICENSE.rst) |
|     `freetype-py`     | latest   | [BSD 3-Clause License](https://github.com/rougier/freetype-py/blob/master/LICENSE.txt) |
|        `gemmi`        | 0.5.8    | [Mozilla Public License 2.0](https://github.com/project-gemmi/gemmi/blob/master/LICENSE.txt) |
|        `numpy`        | 1.26.0   | [BSD 3-Clause License](https://github.com/numpy/numpy/blob/main/LICENSE.txt) |
|      `networkx`       | 2.8.8    | [BSD 3-Clause License](https://github.com/networkx/networkx/blob/main/LICENSE.txt) |
|       `pycifrw`       | 4.4.5    | [PSF License, Version 2](https://github.com/jamesrhester/pycifrw/blob/development/LICENSE) |
|      `pymatgen`       | 2024.3.1 | [MIT License](https://github.com/materialsproject/pymatgen/blob/master/LICENSE) |
|      `pyopengl`       | latest   | [Custom License (based on BSD-3)](https://github.com/Distrotech/PyOpenGL/blob/master/license.txt) |
|       `PySide6`       | latest   | [GNU LGPL v.3](https://doc.qt.io/qt-6/lgpl.html) |
|        `rdkit`        | 2023.9.6 | [BSD 3-Clause License](https://github.com/rdkit/rdkit/blob/master/license.txt) |
|      `requests`       | latest   | [Apache License 2.0](https://github.com/psf/requests/blob/main/LICENSE) |
|      `progress`       | latest   | [ISC License](https://github.com/verigak/progress/blob/master/LICENSE) |
|      `psycopg2`       | latest   | [GNU LGPL v.3](https://github.com/psycopg/psycopg2/blob/master/LICENSE) |
|     `setuptools`      | latest   | [MIT License](https://github.com/pypa/setuptools/blob/main/LICENSE) |

</details>
