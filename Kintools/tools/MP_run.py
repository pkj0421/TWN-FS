import math
import logging
import argparse
import numpy as np
import pandas as pd

from tqdm import tqdm
from pathlib import Path

# ignore warning option
logger = logging.getLogger(__name__)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='BANG grid box')
    parser.add_argument('-twn', '--twn_water', required=True, help='Set your twn folder')
    parser.add_argument('-cmol', '--compare_molecule', required=True, help='Set your compare molecule')
    parser.add_argument('-out', '--output', required=True, help='Set your output folder')
    args = parser.parse_args()

    # set output
    output_path = Path(rf"{args.output}")

    # cmol
    cmol = Path(rf"{args.compare_molecule}")
    # cmol = Path(r"D:\PARK\Lab\HRY\BANG_gridbox\data\3nus.mol2")
    with open(cmol, 'r') as c:
        push = False
        c_crds = {f'{cmol.stem}': []}
        c_lines = c.readlines()
        for c_line in c_lines:
            if c_line.startswith('@<TRIPOS>ATOM'):
                push = True
                continue

            if c_line.startswith('@<TRIPOS>BOND'):
                push = False

            if push:
                c_ats = [l for l in c_line.rstrip().split(' ') if l != ''][2:5]
                c_crds[cmol.stem] += [list(map(float, c_ats))]
    logger.info(f'Read {cmol.stem}')

    # Read TWNs
    TWN = Path(rf'{args.twn_water}')
    MP_out = output_path / 'MP_results'
    MP_out.mkdir(parents=True, exist_ok=True)

    # read pdb
    TWN_fds = [fd for fd in TWN.glob('./*')]
    for fd in TWN_fds:
        fd_name = fd.stem
        fd_out = MP_out / f"{fd_name}_MP.tsv"
        fd_out.mkdir(parents=True, exist_ok=True)
        logger.info(f'Processing {fd_name}...')

        MPs = []
        TWNs = [file for file in fd.glob('./*.pdb')]
        for file in tqdm(TWNs, desc='Calculating MP with TWN'):
            mp_cnt = 0
            twn_box = []
            twn_name = file.stem
            with open(file, 'r') as t:
                lines = t.readlines()
                for line in lines:
                    ats = [l for l in line.rstrip().split(' ') if l != '']
                    if ats[2] == 'OW':
                        twn_crd = list(map(float, ats[5:8]))
                        twn_box += [np.array(twn_crd)]
                        for crd in c_crds[cmol.stem]:
                            if math.dist(crd, twn_crd) <= 1:
                                mp_cnt += 1
            logger.info('Calculating TCP|DA, TCP|DM...')
            twn_cp = np.mean(twn_box, axis=0)
            twn_da = round(np.mean([math.dist(twn_cp, rc) for rc in c_crds[cmol.stem]]), 2)
            twn_dm = round(np.median([math.dist(twn_cp, rc) for rc in c_crds[cmol.stem]]), 2)
            MPs += [{'Molecule': cmol.stem, 'TWN': twn_name, 'MP': mp_cnt, 'TCP|DA': twn_da, 'TCP|DM': twn_dm}]
        MP_result = pd.DataFrame(MPs)
        MP_result.to_csv(fd_out, sep='\t')
        logger.info(f'Saved {fd_out.stem}')





