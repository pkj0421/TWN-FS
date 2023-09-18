import sys
import logging
import argparse
import subprocess
from pathlib import Path

# ignore warning option
logger = logging.getLogger(__name__)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='TWN automation')
    parser.add_argument('-twn', '--twn', required=True, help='Set your parent path of twn folder')
    parser.add_argument('-region', '--region', required=True, help='Set your region folder path')
    parser.add_argument('-shaep', '--shaep', required=True, help='Set your shaep program path')
    parser.add_argument('-ref', '--ref_form', required=True, help='Set your reference form path')
    parser.add_argument('-out', '--output', required=True, help='Set your output folder path')
    args = parser.parse_args()

    path = Path(rf"{args.twn}")
    twn_folders = [f for f in path.glob('*/')]

    region = Path(rf"{args.region}")
    shape = Path(rf"{args.shaep}")
    ref = Path(rf"{args.ref_form}")
    out = Path(rf"{args.output}")

    for idx, twn in enumerate(twn_folders):
        run = f"-twn {twn.as_posix()} -region {region.as_posix()} -shaep {shape.as_posix()} -ref {ref.as_posix()} -out {out.as_posix()}"
        subprocess.run(args=[sys.executable, 'TWN_analysis.py'] + run.split(' '))
        logger.info(f'Finished : {twn.stem} folder ({idx + 1}/{len(twn_folders)})')


