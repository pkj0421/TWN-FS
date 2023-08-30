import math
import logging
import argparse
import itertools
import numpy as np
import pandas as pd
import multiprocessing


from tqdm import tqdm
from pathlib import Path
from functools import partial

# ignore warning option
logger = logging.getLogger(__name__)


# calculate 3d distance
def cal_distance(PV, w, com):
    a, b, c = PV
    x, y, z = w
    e, f, g = com
    d = abs(a * (x - e) + b * (y - f) + c * (z - g)) / ((a ** 2 + b ** 2 + c ** 2) ** (1 / 2))
    return d


# Calculate molecule's Cross Product
def Region_CP(region_mols, region_name):
    with open(region_mols, 'r') as sr:
        lines = sr.readlines()
        read_part = False
        result = []
        for l_idx, line in enumerate(lines):
            if line.startswith('@<TRIPOS>MOLECULE'):
                mol = lines[l_idx + 1].rstrip()
                comp_x = []
                comp_y = []
                comp_z = []
                MCP = []
                continue

            if line.startswith('@<TRIPOS>ATOM'):
                read_part = True
                continue

            if line.startswith('@<TRIPOS>BOND'):
                # pick center point & calculate euclidean distance
                center = [np.mean(comp_x), np.mean(comp_y), np.mean(comp_z)]
                MCP = pd.DataFrame(MCP)
                MCP['Dist'] = MCP[['X', 'Y', 'Z']].apply(lambda c: math.dist(c, center), axis=1)

                # select three far atoms (FA : Far Atom)
                MCP.sort_values(by=['Dist'], axis=0, ascending=False, inplace=True)
                FA1, FA2, FA3 = MCP[MCP.columns[3:6]].iloc[:3].values

                # Calculate Cross Product
                V1 = FA1 - FA2
                V2 = FA1 - FA3
                CP = np.cross(V1.tolist(), V2.tolist())
                result += [{'Molecule': mol, 'Center': center, 'X': CP[0], 'Y': CP[1], 'Z': CP[2]}]
                read_part = False
                continue

            if read_part:
                atoms = line.split()
                dic_atom = {'Molecule': mol, 'Atom_number': atoms[0], 'Atom_name': atoms[1],
                            'X': float(atoms[2]), 'Y': float(atoms[3]), 'Z': float(atoms[4])}
                MCP += [dic_atom]

                comp_x += [dic_atom['X']]
                comp_y += [dic_atom['Y']]
                comp_z += [dic_atom['Z']]
        result = pd.DataFrame(result)
        result['Region'] = region_name
        return result


def DistAVG(gp, rcp):
    result = []
    for _, R_row in rcp.iterrows():
        R_name = R_row['Molecule']
        R_region = R_row['Region']
        R_CP = tuple(R_row[2:5])
        R_center = tuple(R_row[1])
        gp['Dist'] = gp[['X', 'Y', 'Z']].apply(lambda c: cal_distance(R_CP, c, R_center), axis=1)
        gp['F_Dist'] = gp['Dist'].apply(lambda d: d - gp['Dist'].min())
        row = {'Region': R_region, 'Molecule': R_name, 'TWN_Group': gp['MODEL'].values[0]}
        row['Dist_AVG'] = round(gp['F_Dist'].mean(), 2)
        result += [row]
    return result



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='BANG grid box')
    parser.add_argument('-twn', '--twn_water', required=True, help='Set your twn folder')
    parser.add_argument('-cmol', '--compare_molecule', required=True, help='Set your compare molecule')
    parser.add_argument('-gp', '--group', action='store_true', help='If you want calculate twn group')
    parser.add_argument('-out', '--output', required=True, help='Set your output folder')
    args = parser.parse_args()

    # set output
    output_path = Path(rf"{args.output}")

    # cmol
    cmol = Path(rf"{args.compare_molecule}")
    # cmol = Path(rf"D:\PARK\Lab\HRY\BANG_gridbox\data\Kintools\Region")

    RCP_dic = []
    Regions = [rg for rg in cmol.glob('./*.mol2')]
    for rg_file in tqdm(Regions, desc='Reading subregion mols...'):
        RCP_dic += [Region_CP(rg_file, rg_file.stem)]
    RCP = pd.concat(RCP_dic).reset_index(drop=True)

    R_CP = RCP[['X', 'Y', 'Z']].values[0]
    R_Center = RCP['Center'].values[0]

    # Read TWNs
    TWN = Path(rf'{args.twn_water}')
    DA_out = output_path / 'DA_results'
    DA_out.mkdir(parents=True, exist_ok=True)

    # read pdb
    TWN_fds = [fd for fd in TWN.glob('./*')]
    for fd in TWN_fds:
        fd_name = fd.stem
        fd_out = DA_out / f"{fd_name}_DA.tsv"
        logger.info(f'Processing {fd_name}...')

        TWNs = [file for file in fd.glob('./*.pdb')]
        if not args.group:
            twn_df = []
            for file in tqdm(TWNs, desc='Calculating DA with TWN'):
                dist_box = []
                twn_box = []
                twn_name = file.stem
                with open(file, 'r') as t:
                    lines = t.readlines()
                    for line in lines:
                        ats = [l for l in line.rstrip().split(' ') if l != '']
                        if ats[2] == 'OW':
                            twn_crd = list(map(float, ats[5:8]))
                            twn_box += [{'MODEL': twn_name, 'X': twn_crd[0], 'Y': twn_crd[1], 'Z': twn_crd[2]}]
                twn_df += [pd.DataFrame(twn_box)]
            twn_merge = [tm[1] for tm in pd.concat(twn_df).reset_index(drop=True).groupby('MODEL')]

            func = partial(DistAVG, rcp=RCP)
            pool = multiprocessing.Pool(multiprocessing.cpu_count())
            result = pool.map(func, tqdm(twn_merge))
            result = pd.DataFrame(sum(result, []))
            pool.close()
            pool.join()

            result.to_csv(fd_out, sep='\t', index=False)
            logger.info(f'Saved {fd_out.stem}')

        else:
            for file in tqdm(TWNs, desc='Calculating DA with TWN Group'):
                file = Path(r"D:\PARK\Lab\HRY\BANG_gridbox\data\BANG_results\TWN_freq1.pdb")
                with open(file, 'r') as t:
                    lines = t.readlines()
                    model_box = []
                    for line in lines:
                        if line.startswith('MODEL'):
                            mdn = [l for l in line.rstrip().split(' ') if l != ''][1]

                        if line.startswith('HETATM'):
                            ats = [l for l in line.rstrip().split(' ') if l != '']
                            ats_crd = list(map(float, ats[6:9]))
                            model_box += [{'MODEL': mdn, 'X': ats_crd[0], 'Y': ats_crd[1], 'Z': ats_crd[2]}]

                model_gps = [gp[1] for gp in pd.DataFrame(model_box).groupby('MODEL')]

                func = partial(DistAVG, rcp=RCP)
                pool = multiprocessing.Pool(multiprocessing.cpu_count())
                result = pool.map(func, tqdm(model_gps))
                result = pd.DataFrame(sum(result, []))
                pool.close()
                pool.join()

                result.to_csv(fd_out, sep='\t', index=False)
                logger.info(f'Saved {fd_out.stem}')



