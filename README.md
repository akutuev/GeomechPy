[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![PyPI version](https://img.shields.io/pypi/v/geomechpy.svg)](https://pypi.org/project/geomechpy/)
[![Python](https://img.shields.io/pypi/pyversions/geomechpy.svg)](https://pypi.org/project/geomechpy/)

# GeomechPy

**GeomechPy** is a Python library providing a comprehensive set of functions to build a Geomechanics workflow for petroleum rock mechanics applications.

Current modules cover:
- **Elastic Properties** — convert any pair of elastic properties into all other notations (Young's modulus, bulk modulus, shear modulus, Poisson's ratio, etc.)
- **Static Elastic Properties** — velocity/slowness to static elastic property conversion
- **Pore Pressure** — gradient-based pore pressure calculations for onshore and offshore settings
- **Overburden Stress** — overburden stress computation from density logs
- **Stress Calculations** — in-situ stress magnitude and orientation calculations
- **Near Wellbore Stresses** — Kirsch solution for borehole wall and radial stress distribution for any borehole orientation
- **Rock Strength** — rock strength property calculations
- **Wellbore Stability** — wellbore stability analysis

## Installation

```bash
pip install geomechpy
```

## Usage

### Elastic Properties

```python
from geomechpy.elastic_properties import ElasticPropertiesConverter

# Convert Young's modulus + Poisson's ratio into all elastic property notations
props = ElasticPropertiesConverter.convert_from_youngs_and_poissons(
    youngs_modulus=1e6,  # psi
    poissons_ratio=0.25
)
print(props.bulk_modulus)
print(props.shear_modulus)
```

### Pore Pressure

```python
from geomechpy.pore_pressure import PorePressureCalculation

# Onshore pore pressure from TVD and gradient
pp = PorePressureCalculation.calculate_pore_pressure_onshore(
    tvd=5000,                             # ft
    formation_pore_pressure_gradient=0.47 # psi/ft
)
print(pp)  # psi
```

### Near Wellbore Stresses

```python
import numpy as np
from geomechpy.near_wellbore_stresses import NearWellboreStressesCalculation

theta = np.linspace(0, 360, 360)  # azimuthal positions around borehole

stresses = NearWellboreStressesCalculation.calculate_kirsch_borehole_wall_stresses(
    shmin=8000, shmax=10000, svert=12000,  # psi
    pore_pressure=4000,                    # psi
    shmax_azimuth=45,                      # deg
    mud_pressure=4500,                     # psi
    theta=theta,                           # deg
    poisson_ratio_static=0.25,
    borehole_deviation=0,                  # deg (vertical well)
    borehole_azimuth=0                     # deg
)
print(stresses.sigma_boreholewall_tt)      # tangential stress
```

## License

This project is licensed under the terms of the **GNU Lesser General Public License v3.0**.  
See the [LICENSE](./LICENSE) file for details.
