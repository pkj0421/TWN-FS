import sys
import logging
import subprocess
from pathlib import Path

# ignore warning option
logger = logging.getLogger(__name__)

if __name__ == "__main__":

    path = Path(r"parent path of twn folder")
    twn_folders = [f for f in path.glob('*/')]

    region = Path(f"region folder path")
    shape = Path(f"shaep program path")
    ref = Path(f"reference form file path")
    out = Path(f"output folder path")

    for idx, twn in enumerate(twn_folders):
        run = f"-tw {twn.as_posix()} -region {region.as_posix()} -shaep {shape.as_posix()} -ref {ref.as_posix()} -out {out.as_posix()}"
        subprocess.run(args=[sys.executable, 'TWN_analysis.py'] + run.split(' '))
        logger.info(f'Finished : {twn.stem} folder ({idx + 1}/{len(twn_folders)})')


