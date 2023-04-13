from argparse import ArgumentParser, RawTextHelpFormatter

from .conversions import convert
from .filters import filter_atoms
from .subsets import get_trajectory_subset
#from .vmd_spells import chainer
from .mdt_spells import split_merged_trajectories
from .structures import Structure
from .imaging import selective_no_jump, smooth_translation

# Define console call for mdtoolbelt
parser = ArgumentParser(description="Call a tool from mdtoolbelt", formatter_class=RawTextHelpFormatter)
subparsers = parser.add_subparsers(help='Name of the tool to be used', dest="tool")

# The convert command
convert_parser = subparsers.add_parser("convert",
    help="Convert a structure and/or several trajectories to other formats\n" +
        "If several input trajectories are passed they will be merged previously.")
convert_parser.add_argument(
    "-is", "--input_structure",
    help="Path to input structure file")
convert_parser.add_argument(
    "-os", "--output_structure",
    help="Path to output structure file")
convert_parser.add_argument(
    "-it", "--input_trajectories", nargs='*',
    help="Path to input trajectory file(s)")
convert_parser.add_argument(
    "-ot", "--output_trajectory",
    help="Path to output trajectory file")

# The filter command
filter_parser = subparsers.add_parser("filter",
    help="Filter atoms in a structure and/or a trajectory\n")
filter_parser.add_argument(
    "-is", "--input_structure", required=True,
    help="Path to input structure file")
filter_parser.add_argument(
    "-os", "--output_structure",
    help="Path to output structure file")
filter_parser.add_argument(
    "-it", "--input_trajectory",
    help="Path to input trajectory file")
filter_parser.add_argument(
    "-ot", "--output_trajectory",
    help="Path to output trajectory file")
filter_parser.add_argument(
    "-sel", "--selection_string", required=True,
    help="Atom selection formula")
filter_parser.add_argument(
    "-syn", "--selection_syntax", default='vmd',
    help="Atom selection formula (vmd by default)")

# The chainer command
chainer_parser = subparsers.add_parser("chainer",
    help="Set chains on demand in a pdb file.\n" +
        "An atom selection may be provided. Otherwise all the structure is affected.\n" +
        "A chain letter may be passed. Otheriwse a 'chain by fragment' logic is triggered.")
chainer_parser.add_argument(
    "-is", "--input_structure", required=True,
    help="Path to input structure file (pdb)")
chainer_parser.add_argument(
    "-sel", "--atom_selection",
    help="Selection string in VMD syntax for atoms to be chained")
chainer_parser.add_argument(
    "-chn", "--chain_letter",
    help="The chain letter to be set in selected atoms")
chainer_parser.add_argument(
    "-os", "--output_structure",
    help="Path to output structure file (pdb)")

# The split command
split_parser = subparsers.add_parser("split",
    help="Split a trajectory which is the merge of several replicas into the individual replicas.\n" +
        "Splits are guessed by the RMSD profile of the trajectory. A cutoff may be provided.")
split_parser.add_argument(
    "-is", "--input_structure", required=True,
    help="Path to input structure file")
split_parser.add_argument(
    "-it", "--input_trajectory",
    help="Path to input trajectory file")
split_parser.add_argument(
    "-cut", "--cutoff", type=float, default=0.2,
    help="The minimum size of the RMSD jump to consider it is a different trajectory")
split_parser.add_argument(
    "-otp", "--output_trajectory_prefix", default="split",
    help="Prefix for the path to output trajectory files")

# The subset command
subset_parser = subparsers.add_parser("subset",
    help="Get a subset of frames from the current trajectory.")
subset_parser.add_argument(
    "-is", "--input_structure", required=True,
    help="Path to input structure file")
subset_parser.add_argument(
    "-it", "--input_trajectory",
    help="Path to input trajectory file")
subset_parser.add_argument(
    "-ot", "--output_trajectory",
    help="Path to output trajectory file")
subset_parser.add_argument(
    "-start", "--start", type=int, default=0,
    help="Start frame")
subset_parser.add_argument(
    "-end", "--end", type=int, default=None,
    help="End frame")
subset_parser.add_argument(
    "-step", "--step", type=int, default=1,
    help="Frame step")
subset_parser.add_argument(
    "-skip", "--skip", nargs='*', type=int, default=[],
    help="Frames to be skipped")

# The nojump command
nojump_parser = subparsers.add_parser("nojump",
    help="Image a specific atom selection by applying the gromacs 'trjconv -pbc nojump' logic.\n" +
        "WARNING: This may result in periodic boundary artifacts which may not be 'legal'.\n" +
        "WARNING: This logic will not work more than once in a same trajectory.")
nojump_parser.add_argument(
    "-is", "--input_structure", required=True,
    help="Path to input structure file")
nojump_parser.add_argument(
    "-it", "--input_trajectory", required=True,
    help="Path to input trajectory file")
nojump_parser.add_argument(
    "-sel", "--atom_selection", required=True,
    help="Selection string in VMD syntax for atoms to be chained")
nojump_parser.add_argument(
    "-ot", "--output_trajectory",
    help="Path to output trajectory file")

# The smootrans command
smootrans_parser = subparsers.add_parser("smootrans",
    help="Image a trajectory by applying a smooth translation along the trajectory.")
smootrans_parser.add_argument(
    "-is", "--input_structure", required=True,
    help="Path to input structure file")
smootrans_parser.add_argument(
    "-it", "--input_trajectory", required=True,
    help="Path to input trajectory file")
smootrans_parser.add_argument(
    "-trans", "--translation", nargs=3, type=float, required=True,
    help="Translation to be applided (x,y,z). e.g. -trans 1 2.3 -4")
smootrans_parser.add_argument(
    "-ot", "--output_trajectory",
    help="Path to output trajectory file")

# The sumup command
sumup_parser = subparsers.add_parser("sumup",
    help="Get a summary of a pdb")
sumup_parser.add_argument(
    "-is", "--input_structure", required=True,
    help="Path to input structure file")

args = parser.parse_args()

def call():

    # If no command is passed print help
    tool = args.tool
    if not tool:
        parser.print_help()
        return

    # In case the convert tool was called
    if tool == 'convert':
        # If no input arguments are passed print help
        if args.input_structure == None and args.input_trajectories == None:
            convert_parser.print_help()
            return
        if args.input_trajectories == None:
            args.input_trajectories = []
        # Run the convert command
        convert(
            input_structure_filename=args.input_structure,
            output_structure_filename=args.output_structure,
            input_trajectory_filenames=args.input_trajectories,
            output_trajectory_filename=args.output_trajectory,
        )

    # In case the filter tool was called
    if tool == 'filter':
        # Run the convert command
        filter_atoms(
            input_structure_filename=args.input_structure,
            output_structure_filename=args.output_structure,
            input_trajectory_filename=args.input_trajectory,
            output_trajectory_filename=args.output_trajectory,
            selection_string=args.selection_string,
            selection_syntax=args.selection_syntax
        )

    # In case the chainer tool was called
    if tool == 'chainer':
        structure = Structure.from_pdb_file(args.input_structure)
        selection = structure.select(args.atom_selection, syntax='vmd') if args.atom_selection else structure.select_all()
        structure.chainer(selection, args.chain_letter)
        output_structure = args.output_structure if args.output_structure else args.input_structure
        structure.generate_pdb_file(output_structure)

    if tool == 'split':
        split_merged_trajectories(
            input_structure_filename=args.input_structure,
            input_trajectory_filename=args.input_trajectory,
            sudden_jump_cutoff=args.cutoff,
            output_trajectory_prefix=args.output_trajectory_prefix
        )

    if tool == 'subset':
        output_trajectory = args.output_trajectory if args.output_trajectory else args.input_trajectory
        get_trajectory_subset(
            input_structure_filename=args.input_structure,
            input_trajectory_filename=args.input_trajectory,
            output_trajectory_filename=output_trajectory,
            start=args.start,
            end=args.end,
            step=args.step,
            skip=args.skip
        )

    if tool == 'nojump':
        structure = Structure.from_pdb_file(args.input_structure)
        selection = structure.select(args.atom_selection, syntax='vmd')
        output_trajectory = args.output_trajectory if args.output_trajectory else args.input_trajectory
        selective_no_jump(
            input_structure_filename=args.input_structure,
            input_trajectory_filename=args.input_trajectory,
            output_trajectory_filename=output_trajectory,
            selection=selection
        )

    if tool == 'smootrans':
        output_trajectory = args.output_trajectory if args.output_trajectory else args.input_trajectory
        smooth_translation(
            input_structure_filename=args.input_structure,
            input_trajectory_filename=args.input_trajectory,
            output_trajectory_filename=output_trajectory,
            translation=tuple(args.translation)
        )

    if tool == 'sumup':
        structure = Structure.from_pdb_file(args.input_structure)
        structure.display_summary()


    # Tool will always match one of the previous defined options
    # Otherwise argparse returns error itself
    print('Done :)')