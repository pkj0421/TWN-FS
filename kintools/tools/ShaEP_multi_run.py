import sys
import logging
import argparse
import subprocess

from pathlib import Path

# ignore warning option
logger = logging.getLogger(__name__)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='BANG grid box')
    parser.add_argument('-tf', '--twn_folders', required=True, help='Set your twn folders')
    parser.add_argument('-cmol', '--compare_molecule', required=True, help='Set your compare molecule')
    parser.add_argument('-shaep', '--shaep', required=True, help='Set your ShaEP program')
    parser.add_argument('-out', '--output', required=True, help='Set your output folder')
    args = parser.parse_args()

    tf = Path(rf"{args.twn_folders}")
    out = Path(rf"{args.output}")
    shaep = Path(rf"{args.shaep}")
    cmol = Path(rf"{args.compare_molecule}")

    tfs = [fd for fd in tf.glob('./*')]
    for fd in tfs:
        fd_out = out / fd.stem
        fd_out.mkdir(parents=True, exist_ok=True)

        shaep_run = f"-twn {fd.as_posix()} -shaep {shaep.as_posix()} -out {fd_out.as_posix()} -cmol {cmol.as_posix()}"
        subprocess.run(args=[sys.executable, 'ShaEP_run.py'] + shaep_run.split(' '))
        logger.info(f"Finished : {fd.stem}")
    logger.info(f"Job 'ShaEP multi run' Done.")


