[![PyPI version](https://img.shields.io/pypi/v/m2m-postaviz.svg)](https://pypi.org/project/m2m-postaviz/) [![GitHub license](https://img.shields.io/github/license/AuReMe/metage2metabo-postaviz.svg)](https://github.com/AuReMe/metage2metabo-postaviz/blob/main/LICENSE) [![Build 3.12](https://img.shields.io/github/check-runs/AuReMe/metage2metabo-postaviz/main?checkName=build-3.12)](https://github.com/AuReMe/metage2metabo-postaviz/actions/workflows/pythonpackage.yml) [![Build 3.13](https://img.shields.io/github/check-runs/AuReMe/metage2metabo-postaviz/main?checkName=build-3.13)](https://github.com/AuReMe/metage2metabo-postaviz/actions/workflows/pythonpackage.yml) [![Build 3.14](https://img.shields.io/github/check-runs/AuReMe/metage2metabo-postaviz/main?checkName=build-3.14)](https://github.com/AuReMe/metage2metabo-postaviz/actions/workflows/pythonpackage.yml) [![Documentation Status](https://readthedocs.org/projects/metage2metabo-postaviz/badge/?version=latest)](https://metage2metabo-postaviz.readthedocs.io/en/latest/?badge=latest)

# Metage2Metabo-PostAViz (M2M-PostAViz)

M2M-PostAViz (_M2M Post-Analysis and Visualization_) is an interactive platform for exploring metabolic potential predictions from [Metage2Metabo (M2M)](https://github.com/AuReMe/metage2metabo/tree/main).

## Installation

The application is tested on Python versions 3.12, 3.13 and 3.14, on Windows, MacOS and Ubuntu.

Install with pip:

```sh
pip install m2m-postaviz
```

Or from source:

```sh
git clone https://gitlab.inria.fr/postaviz/m2m-postaviz.git
cd m2m-postaviz
pip install -r requirements.txt
pip install .
```

If installation succeeds but running `m2m_postaviz` fails because `polars` is unavailable or incompatible on your machine, install the long-term-support CPU wheel manually:

```sh
pip install polars-lts-cpu
```

## Quickstart

To test the application with example data:

```sh
m2m_postaviz --test
```

For full documentation, usage, and advanced options, see the [online documentation](https://metage2metabo-postaviz.readthedocs.io/).

## License

GNU Lesser General Public License v3 (LGPLv3)

## Authors

Léonard Brindel and [Clémence Frioux](https://cfrioux.github.io)