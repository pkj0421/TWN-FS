import sys
import time
import math
import logging
import argparse
import subprocess
import numpy as np
import pandas as pd
import multiprocessing

from tqdm import tqdm
from rdkit import Chem
from pathlib import Path
from functools import partial

'''
TWN_gridbox.py
parameters :
-twn        TWN folder path (.pdb)
-region     Subregion folder path (.mol2)
-shaep      shaep program path
-ref        Reference form file path (initial columns in Summary, .xlsx)
-out        Output path
'''

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


# parallel translation
def fix_coord(coordinate):
    result = []
    for i in coordinate:
        origin = float(i)
        stem = int(i)
        if round(origin - stem, 2) >= 0.5:
            result += [str(stem + 0.5)]
        else:
            result += [str(stem)]
    return result


# nomination of molecules
def fix_sdf(file):
    with open(file, 'r') as sdf:
        fix_sdf = []
        inform_push = True
        lines = sdf.readlines()
        for line in lines:
            if 'TWN' in line:
                mol_name = line.rstrip().split('/')[-1].split(' ')
                fmol_name = f"{mol_name[0].split('.')[0]}_Group_{mol_name[1]}\n"
                fix_sdf += [fmol_name]
                continue

            if inform_push:
                fix_sdf += [line]

    with open(file, 'w') as fmol2:
        for fm in fix_sdf:
            fmol2.write(fm)


# Store fragment cordinates
def Fragment_cordinates(region_mols):
    # region_mols = Path(r"D:\PARK\Lab\HRY\TWN_gridbox\data\kintools\region\SE.mol2")
    with open(region_mols, 'r') as sr:
        lines = sr.readlines()
        read_part = False
        result = {}
        for l_idx, line in enumerate(lines):
            if line.startswith('@<TRIPOS>MOLECULE'):
                mol = lines[l_idx + 1].rstrip()
                crds = []
                continue

            if line.startswith('@<TRIPOS>ATOM'):
                read_part = True
                continue

            if line.startswith('@<TRIPOS>BOND'):
                crds = pd.DataFrame(crds)
                result[f'{mol}'] = crds[['X', 'Y', 'Z']].astype(float)
                read_part = False
                continue

            if read_part:
                atoms = line.split()
                dic_atom = {'Molecule': mol, 'Atom_number': atoms[0], 'Atom_name': atoms[1],
                            'X': float(atoms[2]), 'Y': float(atoms[3]), 'Z': float(atoms[4])}
                crds += [dic_atom]
        return result


# Read shaep result
def ShaEP_reader(data, out):
    ShaEP_data = pd.read_csv(data, sep='\t').drop(columns=['best_similarity',
                                                           'shape_similarity',
                                                           'ESP_similarity',
                                                           'avg_similarity',
                                                           'N_observations',
                                                           'best_query'])
    ESP_cols = [col for col in ShaEP_data.columns if 'ESP' in col]
    ShaEP_data.drop(columns=ESP_cols, inplace=True)

    # sorted 'molecule' columns
    ShaEP_data['sort'] = ShaEP_data['molecule'].apply(lambda x: int(x.split('_')[-1]))
    ShaEP_data.sort_values(by='sort', inplace=True)

    # round float
    ShaEP_fcols = [c for c in ShaEP_data.keys() if c != 'molecule' and c != 'sort']
    ShaEP_data[ShaEP_fcols] = round(ShaEP_data[ShaEP_fcols], 2)

    # shaep final results
    ShaEP_final = out / f"{ShaEP_sim.stem}_final.tsv"
    with open(ShaEP_final, 'w') as f:
        f.write('Group\tMolecule\tshaep\n')
        for i, row in ShaEP_data.iterrows():
            for sfc in ShaEP_fcols:
                if not np.isnan(row[sfc]):
                    f.write(f"{row['molecule']}\t{sfc.split('_')[0]}\t{row[sfc]}\n")
    return ShaEP_final


# read TWN groups
def TWN_Group_Reader(TG_path, out):
    with open(TG_path, 'r') as wg:
        lines = wg.readlines()
        result = []

        for l_idx, line in enumerate(lines):
            line = [l for l in line.rstrip().split(' ') if l != '']

            if len(line) == 11:
                line = [line[0][:6], line[0][6:]] + line[1:]

            if line[0] == 'MODEL':
                gp_name = f"{TG_path.stem}_Group_{line[1]}"
                continue

            if line[0] == 'HETATM':
                X, Y, Z = list(map(float, line[6:9]))
                gp_prop = {'Group': gp_name, 'TWN FILE': f"{line[3]}{line[5]}", 'X': X, 'Y': Y, 'Z': Z}
                result += [gp_prop]

        result = pd.DataFrame(result)
        ring_size = result[result['Group'] == result['Group'].values[0]]['TWN FILE'].value_counts().tolist()[0]
        result['Number of TWN'] = result['Group'].apply(lambda x : int(result['Group'].value_counts()[x] / ring_size))
        result.to_excel(out, sheet_name='Sheet1', index=False)
    return result.groupby('Group')


def distance_polygon_point(vertices, point):
    # Choose the first three vertices to define the plane
    v0, v1, v2 = vertices[:3]
    normal = np.cross(np.array(v1) - np.array(v0), np.array(v2) - np.array(v0))
    d = -np.dot(normal, np.array(v0))

    # Calculate the distance between the plane and the point
    numerator = abs(np.dot(normal, np.array(point)) + d)
    denominator = math.sqrt(np.dot(normal, normal))
    distance = numerator / denominator
    return distance


# calculate average of distances between ligand centroid and DBSCAN centroid
def Distance_Average(R_box, region, TGs, TGCPs):
    result = []
    R_pdb, R_crds = R_box
    RCP = list(R_crds.mean())  # centroid of fragment

    # Group file
    for TG_name, TG in TGs:
        # DBSCAN Center Points
        DBCPs = [db.split('|') for db in TGCPs[TGCPs['MODEL'] == TG_name]['DBCP'].values[0].split('&')]
        Polygons = [[[x - 100 for x in list(map(float, c.split('-')))] for c in cp] for cp in DBCPs]

        polygon_dists = [distance_polygon_point(Polygon, RCP) for Polygon in Polygons]
        DB_DA = round(np.mean(polygon_dists), 2)

        # DBSCAN distance average
        row = {'region': region, 'Molecule': R_pdb, 'Group': TG_name, 'DBSCAN|Distance Average': DB_DA}
        result += [row]
    return result


# extract top similarity in shaep outputs
def Select_top_value(subrg, ShaEP_output, distbox, sc):
    ShaEP_subrg = ShaEP_output / f"ShaEP_{subrg}"
    ShaEP_files = [sf for sf in ShaEP_subrg.glob('./*.tsv')]

    vals = {}
    pdbs = set()
    for file in ShaEP_files:
        data = pd.read_csv(file, sep='\t', dtype={'Molecule': 'string'})
        groups = data.groupby('Molecule')
        gp_distbox = distbox[subrg].drop(columns='region')
        for gp_key, gp in groups:
            pdbs.add(gp_key)
            gp_dist = pd.merge(gp, gp_distbox, how='left', on=['Group', 'Molecule'])
            max_cols = [c for c in gp_dist.columns if c != 'Molecule']

            # sort the DataFrame by the specified columns and sort order
            gp_sorted = gp_dist.sort_values(by=[col for col, ascending in sc],
                                            ascending=[ascending for col, ascending in sc])
            max_val = gp_sorted.iloc[0]
            vals[f"{max_val['Molecule']}"] = [max_val[mc] for mc in max_cols]

    # generate framework
    result = pd.DataFrame(columns=['pdb'] + max_cols, data={'pdb': list(pdbs)}).set_index('pdb')

    # input data
    for key, val in vals.items():
        for cn, cv in zip(max_cols, val):
            result.at[key, cn] = cv

    result.to_excel(ShaEP_output / f"ShaEP_{subrg}_top_sims.xlsx", sheet_name='Sheet1', header=True)
    top_out = (subrg, result)
    return top_out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='TWN analysis')
    parser.add_argument('-twn', '--twn_water', required=True, help='Set your twn folder')
    parser.add_argument('-region', '--region', required=True, help='Set your region folder')
    parser.add_argument('-shaep', '--shaep', required=True, help='Set your shaep program')
    parser.add_argument('-cond', '--condition', default=False, help='If you want extract condition')
    parser.add_argument('-ref', '--ref_form', required=True, help='Set your reference form')
    parser.add_argument('-out', '--output', required=True, help='Set your output folder')
    args = parser.parse_args()

    # time
    start = time.time()

    # set output
    output_path = Path(rf"{args.output}")

    # set regions
    # you can select specific subregion [AP, FP, GA, SE, X]
    # subregions = ['AP', 'FP', 'GA', 'SE', 'X']
    subregions = ['AP']
    logger.info(f'Subregions : {subregions}')

    # Prepare fragment coordinates for distance calculation with twn
    Fragments = {}
    Region_path = Path(rf'{args.region}')
    Regions = [rg for rg in Region_path.glob('./*.mol2')]
    for rg_file in tqdm(Regions, desc='Reading coordinates of subregion molecules...'):
        Fragments[rg_file.stem] = Fragment_cordinates(rg_file)

    # # Read TWNs

    # initial setting (generate log file and folder)
    TWN = Path(rf'{args.twn_water}')
    twn_out = output_path / f'TWN-FS_{TWN.stem}'
    ob_out = twn_out / 'OBMol'
    twn_out.mkdir(parents=True, exist_ok=True)
    ob_out.mkdir(parents=True, exist_ok=True)

    set_log(twn_out, "TWN_analysis.log")
    logger.info('Start TWN_gridbox.py')

    DA_dic = {}
    twn_run = f'-twn {TWN.as_posix()} -o {twn_out.as_posix()}'
    subprocess.run(args=[sys.executable, 'TWN_gridbox.py'] + twn_run.split(' '))

    # pdb to sdf (change molecular format)
    twn_file = Path(f'{twn_out.as_posix()}/TWN.pdb')
    ob_file = Path(f'{ob_out.as_posix()}/TWN.sdf')
    ob_exe = subprocess.check_output(['where', 'obabel']).decode().rstrip().split('\r\n')[0]
    ob_run = f'{twn_file.as_posix()} -O {ob_file.as_posix()} --errorlevel 1 --addinindex 1'
    subprocess.run(args=[ob_exe] + ob_run.split(' '))
    fix_sdf(ob_file)

    # # shaep & Distance Average
    ShaEP = Path(rf'{args.shaep}')
    ShaEP_out = twn_out / 'ShaEP'
    ShaEP_out.mkdir(parents=True, exist_ok=True)

    DA_out = twn_out / 'Distance_Average'
    DA_out.mkdir(parents=True, exist_ok=True)
    TGs = TWN_Group_Reader(twn_file, twn_file.with_name(f'{twn_file.stem}_Group_inform.xlsx'))
    TGCPs = pd.read_csv(twn_out / f'{twn_file.stem}_Group_Center_Point.tsv', sep='\t')

    # check sdf
    check_file = [m for m in Chem.SDMolSupplier(str(ob_file))]
    lc = len(check_file)
    if lc >= 500:
        sdf_box = []
        na_sp = lc // 500
        logger.info(f"sdf file is too big -> divided lc // 500 = {na_sp} files")
        sdfs = np.array_split(check_file, na_sp)

        d_folder = ob_file.parent / 'D_sdf'
        d_folder.mkdir(parents=True, exist_ok=True)
        for sf_num, sf in enumerate(sdfs):

            d_file = f"{ob_file.stem}_{sf_num}"
            c_file = d_folder / f"{d_file}.sdf"
            with Chem.SDWriter(str(c_file)) as w:
                for m in sf:
                    w.write(m)
            sdf_box += [c_file]
    else:
        sdf_box = [ob_file]

    # start calculations
    for region in subregions:

        ShaEP_sim = ShaEP_out / f'TWN_ShaEP_{region}_Similarity.txt'
        ShaEP_region = Region_path / f'{region}.mol2'
        ShaEP_region_out = ShaEP_out / f'ShaEP_{region}'
        ShaEP_region_out.mkdir(parents=True, exist_ok=True)

        if len(sdf_box) >= 2:
            sim_combine = []
            for d_sdf in sdf_box:
                D_sim = d_sdf.parent / f'{d_sdf.stem}_{region}.txt'
                sim_combine += [D_sim]

                ShaEP_run = f'--input-file {d_sdf.as_posix()} --q {ShaEP_region.as_posix()} ' \
                            f'--output-file {D_sim.as_posix()} --noOptimization'
                subprocess.run(args=[str(ShaEP)] + ShaEP_run.split(' '), shell=True)
                logger.info(f"Processed shaep : {d_sdf} file")

            logger.info('Combining sims...')
            comb_sims = pd.concat([pd.read_csv(f, sep='\t') for f in sim_combine]).reset_index(drop=True)
            comb_sims.to_csv(ShaEP_sim, sep='\t', index=False)

        else:
            ShaEP_run = f'--input-file {ob_file.as_posix()} --q {ShaEP_region.as_posix()} ' \
                        f'--output-file {ShaEP_sim.as_posix()} --noOptimization'
            subprocess.run(args=[str(ShaEP)] + ShaEP_run.split(' '), shell=True)

        # # Calculate shaep
        ShaEP_file = ShaEP_reader(ShaEP_sim, ShaEP_region_out)
        logger.info(f'Saved {region} shaep file')

        # # Calculate distance average (primary grouping & secondary grouping)
        pool_DA = multiprocessing.Pool(multiprocessing.cpu_count())
        DA = partial(Distance_Average, region=region, TGs=TGs, TGCPs=TGCPs)

        # # Calculate distance average
        DA_tqdm = tqdm(Fragments[region].items(),
                         desc=f'Calculating {region} Distance Average Point...')
        DA_result = pd.DataFrame(sum(pool_DA.map(DA, DA_tqdm), []))
        pool_DA.close()
        pool_DA.join()

        DA_result['Molecule'] = DA_result['Molecule'].astype('string')
        DA_dic[f'{region}'] = DA_result
        DA_result.to_csv(DA_out / f'{twn_file.stem}_{region}_Distance_Average.tsv', sep='\t', index=False)
        logger.info(f'Saved {region}_Distance_Average file')

    # # extract top
    logger.info(f'Selecting top values...')

    # prompt the user to enter a column name to sort by
    logger.info(f"Column names : {['shaep'] + DA_result.columns.tolist()}")
    logger.info('Input example : shaep-desc, DBSCAN|Distance Average-asc')

    if args.condition:
        # for log
        time.sleep(3)
        while True:
            sort_cols = input('Enter column names-sort order (asc or desc), separated by commas : ')

            try:
                # convert the input string to a list of tuples with column names and sort order
                sort_cols = [col.strip().split('-') for col in sort_cols.split(',')]
                sort_cols = [(col, True if order.strip().lower() == 'asc' else False) for col, order in sort_cols]
                print(f'\nYour command : {sort_cols}')
            except:
                print('\nWrong input format')
                continue

            # check input
            check = True
            for col, asending in sort_cols:
                if (col == 'shaep') and (asending in [True, False]):
                    continue
                elif (col in DA_result.columns) and (asending in [True, False]):
                    continue
                else:
                    check = False

            if check:
                break
            else:
                print("That's not the correct. Please try again.")
    else:
        sort_cols = [('shaep', False), ('DBSCAN|Distance Average', True)]
    logger.info(f'Sorting columns : {sort_cols}')

    pool_ST = multiprocessing.Pool(len(subregions))
    ST = partial(Select_top_value, ShaEP_output=ShaEP_out, distbox=DA_dic, sc=sort_cols)
    ShaEP_tops = dict(pool_ST.map(ST, subregions))
    pool_ST.close()
    pool_ST.join()
    logger.info('Extracted top similarity.')

    # # Summary
    ref = Path(rf"{args.ref_form}")
    logger.info('Writing Summary...')
    with pd.ExcelWriter(twn_out / 'Summary.xlsx') as writer:
        ref_form = pd.read_excel(ref, sheet_name='Sheet1', dtype={'PDB': 'string'}).drop(columns=['ID'])
        for region in subregions:
            ref_region_form = ref_form[ref_form['region'] == region].drop(columns=['region'])
            ref_region_form['PDB'] = ref_region_form['PDB'].apply(lambda x: x.replace('_', '-'))

            # region Tops
            RT = ShaEP_tops[region]
            RT_merge = pd.merge(ref_region_form, RT, how='left', left_on='PDB', right_on='pdb')

            # Save summary in region sheet
            RT_merge.to_excel(writer, sheet_name=region, index=False)
    logger.info('Saved Summary')

    # finish time
    f_time = time.strftime("%H:%M:%S", time.gmtime(time.time() - start))
    logger.info('All process have been completed')
    logger.info(f'Learning time : {f_time}')

