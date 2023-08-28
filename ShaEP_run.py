import time
import logging
import argparse
import subprocess
import numpy as np
import pandas as pd

from pathlib import Path

# ignore warning option
logger = logging.getLogger(__name__)
pd.set_option('mode.chained_assignment', None)


# set logger
def set_log(path_output, log_message):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(path_output / log_message),
            logging.StreamHandler()
        ]
    )



# named molecules
def fix_sdf(file):
    with open(file, 'r') as sdf:
        fix_sdf = []
        inform_push = True
        lines = sdf.readlines()
        for line in lines:
            if 'pdb' in line:
                mol_name = line.rstrip().split('/')[-1].split('.')[0]
                fix_sdf += [mol_name + '\n']
                continue

            if inform_push:
                fix_sdf += [line]

    with open(file, 'w') as fmol2:
        for fm in fix_sdf:
            fmol2.write(fm)


# find twn
def find_twn(file, sch_lst, fto):
    with open(file, 'r') as sdf:
        f_sdf = []
        inform_push = False
        lines = sdf.readlines()
        for line in lines:
            line = line.rstrip()
            if line in sch_lst:
                inform_push = True
                continue

            if inform_push:
                f_sdf += [line + '\n']
                if line == '$$$$':
                    inform_push = False

    with open(fto, 'w') as fmol2:
        for fm in f_sdf:
            fmol2.write(fm)


# Read ShaEP result
def ShaEP_reader(data, out):
    ShaEP_data = pd.read_csv(data, sep='\t').drop(columns=['best_similarity',
                                                           'shape_similarity',
                                                           'ESP_similarity',
                                                           'avg_similarity',
                                                           'N_observations',
                                                           'best_query'])
    ESP_cols = [col for col in ShaEP_data.columns if 'ESP' in col]
    ShaEP_data.drop(columns=ESP_cols, inplace=True)

    # round float
    ShaEP_fcols = [c for c in ShaEP_data.keys() if c != 'molecule' and c != 'sort']
    ShaEP_data[ShaEP_fcols] = round(ShaEP_data[ShaEP_fcols], 2)

    # ShaEP final results
    ShaEP_final = Path(out / f"{data.stem}_final.tsv")
    with open(ShaEP_final, 'w') as f:
        f.write('TWN\tMolecule\tValue\n')
        for i, row in ShaEP_data.iterrows():
            for sfc in ShaEP_fcols:
                if not np.isnan(row[sfc]):
                    f.write(f"{row['molecule']}\t{sfc.split('_')[0]}\t{row[sfc]}\n")
    return ShaEP_final


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='BANG grid box')
    parser.add_argument('-twn', '--twn_water', required=True, help='Set your twn folder')
    parser.add_argument('-cmol', '--compare_molecule', required=True, help='Set your compare molecule')
    parser.add_argument('-shaep', '--shaep', required=True, help='Set your ShaEP program')
    parser.add_argument('-out', '--output', required=True, help='Set your output folder')
    args = parser.parse_args()

    # time
    start = time.time()

    # set output
    output_path = Path(rf"{args.output}")

    # Read TWNs
    TWN = Path(rf'{args.twn_water}')
    shaep_out = output_path / 'ShaEP_results'
    shaep_out.mkdir(parents=True, exist_ok=True)
    ob_out = shaep_out / 'OBMol'
    ob_out.mkdir(parents=True, exist_ok=True)

    # # ShaEP & MP
    ShaEP = Path(rf'{args.shaep}')
    ShaEP_out = shaep_out / 'ShaEP'
    ShaEP_out.mkdir(parents=True, exist_ok=True)

    # initial setting
    set_log(shaep_out, "ShaEP.log")

    # pdb to sdf
    TWNs = [file for file in TWN.glob('./*')]
    ob_exe = subprocess.check_output(['where', 'obabel']).decode().rstrip().split('\r\n')[0]
    ob_run = f'{TWN.as_posix()}/*.pdb -O {ob_out.as_posix()}/*.sdf'
    subprocess.run(args=[ob_exe] + ob_run.split(' '))

    ob_name = Path(f"{ShaEP_out.as_posix()}/{TWN.stem}.sdf")
    ob_c_run = f'{ob_out.as_posix()}/*.sdf -O {ob_name.as_posix()}'
    subprocess.run(args=[ob_exe] + ob_c_run.split(' '))
    fix_sdf(ob_name)

    # # compare molecule
    cmol = Path(rf'{args.compare_molecule}')
    ShaEP_sim = ShaEP_out / f'{ob_name.stem}_{cmol.stem}_sim.txt'
    ShaEP_run = f'--input-file {ob_name} --q {cmol.as_posix()} ' \
                f'--output-file {ShaEP_sim.as_posix()} --noOptimization'
    subprocess.run(args=[str(ShaEP)] + ShaEP_run.split(' '), shell=True)

    # # Calculate ShaEP
    ShaEP_file = ShaEP_reader(ShaEP_sim, ShaEP_out)
    logger.info('Saved ShaEP file')

    ShaEP_final = pd.read_csv(ShaEP_file, sep='\t')
    ShaEP_final.sort_values(by=['Value'], ascending=False, inplace=True)

    # # find twns
    ShaEP_sch = ShaEP_final[ShaEP_final['Value'] >= 0.7]
    print(ShaEP_sch)
    logging.info(f'Search 0.7 upper : {len(ShaEP_sch)}')
    ft_out = ShaEP_out / 'find_twn_by_sim.sdf'
    find_twn(ob_name, ShaEP_sch['TWN'].tolist(), ft_out)



