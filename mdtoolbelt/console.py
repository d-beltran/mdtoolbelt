from argparse import ArgumentParser, RawTextHelpFormatter

from .conversions import convert
from .subsets import get_trajectory_subset
from .vmd_spells import chainer
from .mdt_spells import split_merged_trajectories

# Define console call for mdtoolbelt
parser = ArgumentParser(description="Call a tool from mdtoolbelt", formatter_class=RawTextHelpFormatter)
subparsers = parser.add_subparsers(help='Name of the tool to be used', dest="tool")

# The convert command
convert_parser = subparsers.add_parser("convert")
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

# The chainer command
chainer_parser = subparsers.add_parser("chainer")
chainer_parser.add_argument(
    "-is", "--input_structure", required=True,
    help="Path to input structure file (pdb)")
chainer_parser.add_argument(
    "-sel", "--atom_selection",
    help="Atoms to be chained")
chainer_parser.add_argument(
    "-chn", "--chain_letter",
    help="The chain letter to be set in selected atoms")
chainer_parser.add_argument(
    "-os", "--output_structure",
    help="Path to output structure file (pdb)")

# The split command
split_parser = subparsers.add_parser("split")
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
subset_parser = subparsers.add_parser("subset")
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
    "-end", "--end", type=int, default=0,
    help="End frame")
subset_parser.add_argument(
    "-step", "--step", type=int, default=1,
    help="Frame step")

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
        # Run the convert command
        convert(
            input_structure_filename=args.input_structure,
            output_structure_filename=args.output_structure,
            input_trajectory_filenames=args.input_trajectories,
            output_trajectory_filename=args.output_trajectory,
        )

    # In case the chainer tool was called
    if tool == 'chainer':
        chainer(
            input_pdb_filename=args.input_structure,
            atom_selection=args.atom_selection,
            chain_letter=args.chain_letter,
            output_pdb_filename=args.output_structure
        )

    if tool == 'split':
        split_merged_trajectories(
            input_structure_filename=args.input_structure,
            input_trajectory_filename=args.input_trajectory,
            sudden_jump_cutoff=args.cutoff,
            output_trajectory_prefix=args.output_trajectory_prefix
        )

    if tool == 'subset':
        get_trajectory_subset(
            input_structure_filename=args.input_structure,
            input_trajectory_filename=args.input_trajectory,
            output_trajectory_filename=args.output_trajectory,
            start=args.start,
            end=args.end,
            step=args.step
        )

    # Tool will always match one of the previous defined options
    # Otherwise argparse returns error itself
    print('Done :)')