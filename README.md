# mdtoolbelt

Tools por MD post-processing

These tools may be used both from the terminal or from python.


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