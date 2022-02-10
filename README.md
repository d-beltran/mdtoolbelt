# mdtoolbelt

Tools por MD post-processing

**Important**: PIP only installs the package. All the dependencies must be installed separately. To perform a complete installation, please use ANACONDA.


## Conversions

Convert structure and trajectory files from one format to another

### Python

```python
from mdtoolbelt.conversions import convert

convert(input_trajectory_filename='trajectory.xtc', output_trajectory_filename='trajectory.dcd')
```

### Bash

```bash
mdtb convert -it trajectory.xtc -ot trajectory.dcd
```