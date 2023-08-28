import re
import math
import logging
import argparse
import operator
import numpy as np
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from sklearn.cluster import DBSCAN

# module info option
logger = logging.getLogger(__name__)


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


# trans format
def trans_format(form, *vals):
    pdb_form = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}"
    ter_form = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}"
    model_form = '{:6s}     {:>3d}'
    trans = {'pdb': pdb_form, 'ter': ter_form, 'model': model_form}
    return trans[form].format(*vals)


# for max_boundary function
def update_boundary(crd):
    if abs(round(crd - int(crd), 2)) >= 0.5:
        return float(str(int(crd)) + '.5')
    else:
        return int(crd)


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


# for broad_match function
def fix_ac(coords):
    result = []
    for coord in coords:
        if coord - int(coord) == 0:
            result += [str(int(coord))]
        else:
            result += [str(coord)]
    return '-'.join(result)


# twn expand search
def broad_match(t_key, twn_inform):
    space = []
    x, y, z = map(float, t_key.split('-'))
    rg = 0.5
    for ax in [x - rg, x, x + rg]:
        for ay in [y - rg, y, y + rg]:
            for az in [z - rg, z, z + rg]:
                ac = fix_ac([ax, ay, az])
                if ac in twn_inform.keys():
                    space += twn_inform[ac]
    space = set(space)
    return space


# merge twn groups
def merger(ct, cts):
    merge_lst = []
    for ct2_key, ct2 in cts:
        merge = True
        for ct1_crd in ct:
            dist_box = []
            for ct2_crd in ct2:
                cl_dist = math.dist(ct1_crd, ct2_crd)
                dist_box += [cl_dist]
            if min(dist_box) >= 1:
                merge = False
        if merge:
            merge_lst += [ct2_key]
    return merge_lst


# update max boundary range
def max_boundary(bd, bd_rg):
    with open(bd, 'r') as b:
        lines = b.readlines()
        x_crds = []
        y_crds = []
        z_crds = []
        for line in lines:
            line = line.rstrip()
            if not 'HETATM' in line:
                continue

            crd = [float(x) for x in line[32:55].split(' ') if x != '']
            x_crds += [crd[0]]
            y_crds += [crd[1]]
            z_crds += [crd[2]]

        # update x coordinates
        xp_max = update_boundary(max(x_crds) + bd_rg) + 100
        xm_max = update_boundary(min(x_crds) - bd_rg) + 100

        # update y coordinates
        yp_max = update_boundary(max(y_crds) + bd_rg) + 100
        ym_max = update_boundary(min(y_crds) - bd_rg) + 100

        # update x coordinates
        zp_max = update_boundary(max(z_crds) + bd_rg) + 100
        zm_max = update_boundary(min(z_crds) - bd_rg) + 100

    return [xp_max, xm_max, yp_max, ym_max, zp_max, zm_max]


# read protein
def Protein_reader(protein_file):
    with open(protein_file, 'r') as f:
        lines = f.readlines()
        protein_inform = []
        for line in lines:
            line = line.rstrip()
            atoms = [p for p in line.split(' ') if p != '']

            if len(atoms) != 11:
                protein_inform += [line]
                continue

            if atoms[3] == 'SOL' or atoms[3] == 'CL':
                continue

            else:
                pdb_row = trans_format('pdb', atoms[0], int(atoms[1]), atoms[2], '', atoms[3], 'A', int(atoms[4]), '',
                                       float(atoms[5]), float(atoms[6]), float(atoms[7]),
                                       float(atoms[8]), float(atoms[9]), atoms[10], '')
                protein_inform += [pdb_row]

    protein_pdb = protein_inform + ['TER', 'ENDMDL']
    return protein_pdb


# read twn water
def TWN_reader(twn, bd_max=False):
    twn_box = {}
    twn_inform = {}
    twn_atom_num = 1
    for twn_file in tqdm(twn):
        twn_name = twn_file.stem
        file_box = {'twn_crd': [], 'twn_box': [], 'twn_inform': []}
        with open(twn_file, 'r') as f:
            draw = True
            lines = f.readlines()
            for line in lines:
                if not draw:
                    break
                line = line.rstrip()
                atoms = [p for p in line.split(' ') if p != '']

                if atoms[2] != 'OW':
                    continue

                if len(atoms) == 10:
                    atoms += ['O']

                if len(atoms) != 11:
                    continue

                coordinate = [c + 100 for c in list(map(float, atoms[5:8]))]
                fix_coordinate = fix_coord(coordinate)

                if bd_max:
                    if (float(fix_coordinate[0]) > bd_max[0]) or (float(fix_coordinate[0]) < bd_max[1]):
                        draw = False
                    if (float(fix_coordinate[1]) > bd_max[2]) or (float(fix_coordinate[1]) < bd_max[3]):
                        draw = False
                    if (float(fix_coordinate[2]) > bd_max[4]) or (float(fix_coordinate[2]) < bd_max[5]):
                        draw = False

                pdb_row = trans_format('pdb', 'HETATM', twn_atom_num, atoms[2], '',
                                       "{0:-<3}".format(twn_name.split('_')[2]), 'W', int(twn_name.split('_')[1]),
                                       '', float(atoms[5]), float(atoms[6]), float(atoms[7]),
                                       float(atoms[8]), float(atoms[9]), atoms[10], '')

                file_box['twn_crd'] += ['-'.join(fix_coordinate)]
                file_box['twn_box'] += [pdb_row]
            twn_atom_num += 1

            if draw:
                for f_crd in file_box['twn_crd']:
                    try:
                        twn_inform[f_crd] += [twn_name]
                    except:
                        twn_inform[f_crd] = [twn_name]

                twn_box[twn_name] = file_box['twn_box']
    return twn_box, twn_inform


# write twn water
def TWN_writer(twn_box, twn_inform, twn_out):

    # twn_inform = TWN_inform
    cluster = []
    for t_key, twn_pocket in twn_inform.items():

        # search twn coordinates
        for tr in twn_pocket:

            # bring gridbox information include atom of twn
            area = [t[0] for t in twn_inform.items() if tr in t[1]]

            # expand area each atom
            space = [broad_match(ar, twn_inform) for ar in area]

            # find intersection atoms
            twn_group = eval(' & '.join([f"space[{s}]" for s in range(0, len(area))]))

            if not twn_group:
                continue

            if (len(twn_group) == 1) or (twn_group in cluster):
                continue
            else:
                cluster += [twn_group]

    # remove sub-clusters
    remove_subsets = [c for c in cluster if len([s for s in cluster if c.issubset(s)]) == 1]
    logger.info(f'Search TWN clusters : {len(cluster)} -> Remove sub-clusters {len(remove_subsets)}')

    # re-grouping cluster (DBSCAN)
    db_atoms = {}
    db_centers = []
    dbscan = DBSCAN(eps=1, min_samples=2)
    for cl_idx, cl in enumerate(remove_subsets):
        gp_crds = [[c[0] for c in twn_inform.items() if d in c[1]] for d in cl]
        gp_ats = sum([[list(map(float, c.split('-'))) for c in crds] for crds in gp_crds], [])
        gp_df = pd.DataFrame([{'x': gpc[0], 'y': gpc[1], 'z': gpc[2]} for gpc in gp_ats])
        gp_df['cluster'] = dbscan.fit_predict(gp_df)

        # check dbscan
        cls_cnt = len(gp_df['cluster'].unique())
        if cls_cnt != len(gp_crds[0]):
            logger.info(f'Failed clustering (DBSCAN) : {cl} is {len(gp_crds[0])}-ring, but clustered {cls_cnt}...')

        # calculate center point each db_cluster
        db_atoms[cl_idx] = gp_df
        db_centers += [(cl_idx, [np.mean(dbc[1].values, axis=0)[:3] for dbc in gp_df.groupby('cluster')])]

    # calculate distance
    regroups = []
    duple_setnum = []
    for cl_key, cl_center in db_centers:
        merge_lst = merger(cl_center, db_centers)
        m_centers = [m for m in db_centers if m[0] in merge_lst]
        for mc_idx, mc_vals in enumerate(m_centers):
            remerge = set(merger(mc_vals[1], m_centers))
            if len(remerge) != len(merge_lst):
                m_set = set().union(*[remove_subsets[r] for r in remerge])
                if not m_set in regroups:
                    regroups += [m_set]
                    duple_setnum += [remerge]

    remove_gp = set().union(*duple_setnum)
    regroup_sets = regroups + [sb[1] for sb in enumerate(remove_subsets) if not sb[0] in remove_gp]
    regroup_sets_num = duple_setnum + [{sb[0]} for sb in enumerate(db_atoms.keys()) if not sb[0] in remove_gp]

    logger.info(f'{len(remove_gp)} groups merging... generate {len(regroups)} new sets.')
    logger.info(f'final sets : {len(remove_subsets)} - {len(remove_gp)} + {len(regroups)} = {len(regroup_sets)}')

    # regroups db_centers
    rdb_centers = []
    for rcl_idx, rcls in enumerate(regroup_sets_num):
        rgp_df = pd.concat([db_atoms[rn] for rn in rcls])
        db_cp = '&'.join(['|'.join(['-'.join([str(c) for c in list(ct)]) for ct in db_centers[rn][1]]) for rn in rcls])

        # calculate center point each db_cluster
        rdb_centers += [(rcl_idx + 1,
                         (rgp_df[['x', 'y', 'z']].sum().to_numpy() / len(rgp_df)) - 100,
                         len(rgp_df['cluster'].unique()), db_cp)]

    # save rdb_centers
    with open(twn_out.parent / f'{twn_out.stem}_Group_Center_Point.tsv', 'w') as rdb:
        rdb.write('MODEL\tTGCP\tRING\tDBCP\n')
        for rdb_idx, tgcp, ring_num, dbcp in rdb_centers:
            tgcp_crds = '|'.join([str(t) for t in tgcp])
            rdb.write(f"{twn_out.stem}_Group_{rdb_idx}\t{tgcp_crds}\t{ring_num}\t{dbcp}\n")
    logger.info(f"Saved RDBC_inform : {twn_out.parent / f'{twn_out.stem}_Group_Center_Point.tsv'}")

    # naming twn cluster
    twn_each = []
    for st_idx, st in enumerate(regroup_sets):
        model_row = trans_format('model', 'MODEL', st_idx + 1)
        twn_each += [[model_row] + sorted(sum([twn_box[d] for d in st], []), key=lambda x: int(x[22:26])) + ['ENDMDL']]
    twn_each += [['ENDMDL']]
    twn_pdb = sum(twn_each, [])

    # write twn
    with open(twn_out, 'w') as f:
        for line in twn_pdb:
            f.write(f"{line}\n")

    return twn_pdb


# for complex pdb
def trans_name(pdb, key):
    for l_idx, l in enumerate(pdb):
        if l.startswith('HETATM'):
            atom_inform = [ls for ls in l.split(' ') if ls != '']
            atom_inform.insert(3, '')
            atom_inform.insert(7, '')
            atom_inform.insert(14, '')
            for a_idx, a in enumerate(atom_inform):
                if a_idx in [1, 6]:
                    atom_inform[a_idx] = int(a)
                if a_idx == 5:
                    atom_inform[a_idx] = key
                if a_idx in [8, 9, 10, 11, 12]:
                    atom_inform[a_idx] = float(a)
            pdb[l_idx] = trans_format('pdb', *atom_inform)
    return pdb


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='BANG grid box')
    parser.add_argument('-pt', '--protein', default=False, help='Set your protein file')
    parser.add_argument('-bd', '--boundary', default=False, help='Set your boundary file')
    parser.add_argument('-tw', '--twn_water', required=True, help='Set your twn folder')
    parser.add_argument('-o', '--output', required=True, help='Set your output folder')
    parser.add_argument('-c', '--complex', action='store_true', help='If you hope to gain complex.pdb')
    args = parser.parse_args()

    # set output
    output_path = Path(rf"{args.output}")
    set_log(output_path, "BANG_gridbox.log")

    # set protein & boundary
    if args.protein:
        Protein = Path(rf"{args.protein}")
        Protein_pdb = Protein_reader(Protein)
        logger.info(f'Set protein : {Protein}')

    if args.boundary:
        boundary_file = Path(rf"{args.boundary}")
        boundary_range = 3
        bd_max = max_boundary(boundary_file, boundary_range)
        logger.info(rf'Limited boundary: Use {boundary_file}')

    # start twn water
    TWN_path = Path(rf"{args.twn_water}")
    TWN = [twn for twn in TWN_path.glob('./*.pdb')]
    TWN_out = output_path / f"TWN.pdb"
    logger.info(f'Set TWN water : {TWN_path} | {len(TWN)} files in folder')

    logger.info(f'Start TWN analysis...')
    if args.boundary:
        TWN_box, TWN_inform = TWN_reader(TWN, bd_max)
    else:
        TWN_box, TWN_inform = TWN_reader(TWN)
    logger.info(f'Generated {len(TWN_inform)} gridbox')

    TWN_pdb = TWN_writer(TWN_box, TWN_inform, TWN_out)
    logger.info(f'Saved {TWN_out}')

    # write complex pdb
    if args.complex:
        TWN_pdb = trans_name(TWN_pdb, 'T')
        Complex_pdb = Protein_pdb + TWN_pdb
        Complex_out = output_path / f"Complex.pdb"
        with open(Complex_out, 'w') as c:
            for line in Complex_pdb:
                c.write(f"{line}\n")
        logger.info(f'Saved Complex pdb : {Complex_out}')

logger.info('Complete.')
