# Set the "standard" file format of every possible file extension
# Note that some formats have different possible extension (e.g. nc, cdf, netcdf)
EXTENSION_FORMATS = {
    # Topologies
    'top': 'top',
    'psf': 'psf',
    'prmtop': 'prmtop',
    'prm7': 'prmtop',
    'txt': 'txt', # charges.txt
    # Structures
    'pdb': 'pdb',
    'gro': 'gro',
    # Trajectories
    'xtc': 'xtc',
    'trr': 'trr',
    'dcd': 'dcd',
    'nc': 'nc',
    'cdf': 'nc',
    'netcdf': 'nc',
    'crd': 'crd',
    'mdcrd': 'crd',
    'trj': 'crd',
    # Other
    'json': 'json',
    'npy': 'npy',
    'in': 'txt'
}

# Topology and trajectory file formats supported by PyTraj
PYTRAJ_SUPPORTED_FORMATS = set([
    # Topologies
    'prmtop', 'top', 'psf', 'pdb'
    # Trajectories
    'nc', 'crd', 'dcd', 'trr', 'xtc'
])

# From GitHub:
# ParmFormatDict = {
#     "AMBERPARM": AMBERPARM,
#     "PDBFILE": PDBFILEPARM,
#     "MOL2FILE": MOL2FILEPARM,
#     "CHARMMPSF": CHARMMPSF,
#     "CIFFILE": CIFFILE,
#     "GMXTOP": GMXTOP,
#     "SDFFILE": SDFFILE,
#     "TINKER": TINKERPARM,
#     "UNKNOWN_PARM": UNKNOWN_PARM,
# }

# Set some flags requeired to write files with pytraj
PYTRAJ_PARM_FORMAT = {
    'prmtop': 'AMBERPARM',
    'psf': 'CHARMMPSF',
    'top': 'GMXTOP',
    'pdb': 'PDBFILE'
}