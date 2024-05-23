import pandas as pd
import sys
import numpy as np
import tqdm
import json
from scipy.stats import binom, norm, chi2
from mass2chem.formula import calculate_formula_mass
from itertools import product
from default_search_config import *


m0_reagent = float(sys.argv[1])
m13c_reagent = float(sys.argv[2])
reagent_mass_delta = float(sys.argv[3])
ft = pd.read_csv(sys.argv[4], sep="\t")
ratio = m0_reagent / (m0_reagent + m13c_reagent)

possible_reagents = [{"num_C": x, "prob": binom.pmf(x, REAGENT_C, C13_PROB), "ratio": ratio} for x in range(REAGENT_C)]
possible_reagents += [{"num_C": x, "prob": binom.pmf(x-2, REAGENT_C-2, C13_PROB), "ratio": 1-ratio} for x in range(2, REAGENT_C)]

reagents_for_groups = {num_groups: [
    {
        "num_C": np.sum([x["num_C"] for x in combo]),
        "prob": np.product([x["prob"] for x in combo]) * np.product([x["ratio"] for x in combo])
    }
    for combo in product(possible_reagents, repeat=num_groups)
    if np.product([x["prob"] for x in combo]) * np.product([x["ratio"] for x in combo]) > PROB_LIMIT
] for num_groups in range(1, NUM_GROUPS + 1)}

prxu = [[binom.pmf(z, REAGENT_C * i, C13_PROB) for z in range(REAGENT_C * i)] for i in range(1, NUM_GROUPS + 1)]
ctrl_pool = [x for x in ft.columns[11:] if "Pool_QC_cellpellets_12C" in x and "Pool_QC_cellpellets_12C13C" not in x]
labl_pool = [x for x in ft.columns[11:] if "Pool_QC_cellpellets_12C13C" in x]

ctrl_medians = ft[ctrl_pool].median(axis=1, skipna=True)
labl_medians = ft[labl_pool].median(axis=1, skipna=True)
ctrl_nonzero = ctrl_medians.replace(0, np.nan)
labl_nonzero = labl_medians.replace(0, np.nan)

#working_ft = ft.copy()

#ft = ft[ft['mz'].between(670, 680)]
possibly_isotoped = ft.loc[labl_medians > ctrl_medians, 'id_number'].tolist()
possibly_mono = ft.loc[ctrl_medians > labl_medians * (1 / ratio * 0.75), 'id_number'].tolist()
ignore = ft.loc[(labl_medians <= ctrl_medians) & (ctrl_medians <= labl_medians * (1 / ratio * 0.75)), 'id_number'].tolist()


results = []
ft_dict = {}
ctrl_vals = {None: [0 for _ in ctrl_pool]}
labl_vals = {None: [0 for _ in labl_pool]}
for x in ft.to_dict(orient='records'):
    ft_dict[x['id_number']] = x
    ctrl_vals[x['id_number']] = [ft_dict.get(x['id_number'], {}).get(z, 0) for z in ctrl_pool]
    labl_vals[x['id_number']] = [ft_dict.get(x['id_number'], {}).get(z, 0) for z in labl_pool]

for f_root in tqdm.tqdm(ft[ft['id_number'].isin(possibly_mono)].to_dict(orient='records')):
    root_mz = f_root['mz']
    max_groups = max(int(np.ceil(root_mz / reagent_mass_delta)), NUM_GROUPS)

    possible = 0
    while possible * 12 < root_mz and possible < 32 - NUM_GROUPS * 2 and possible < MAX_ISO:
        possible += 1

    isotopologues = [[f_root['id_number']]] + [[] for _ in list(range(possible))]
    results_for_root = []


    for i in range(1, possible):
        iso_mz = root_mz + i * C13_MASS_DELTA
        iso_mz_err = iso_mz / 1e6 * 10
        for f_iso in ft[ft['mz'].between(iso_mz - iso_mz_err, iso_mz + iso_mz_err) & ft['rtime'].between(f_root['rtime_left_base'], f_root['rtime_right_base'])].to_dict(orient='records'):
            isotopologues[i].append(f_iso['id_number'])
    
    if sum([len(x) for x in isotopologues]) > 15:
        continue
    

    for index in np.ndindex(tuple([len(x) + 1 for x in isotopologues])):
        if index[0] == 1:
            working_isotopologues = [isotopologues[n][i-1] if i != 0 else None for n, i in enumerate(index)]
            if np.any([x in possibly_isotoped for x in working_isotopologues]) and len(working_isotopologues) >= MIN_CLIQUE_SIZE:
                ctrl_intensity_matrix = [ctrl_vals[f_id] for f_id in working_isotopologues]
                labl_intensity_matrix = [labl_vals[f_id] for f_id in working_isotopologues]
                norm_labl = labl_intensity_matrix / np.sum(labl_intensity_matrix, axis=0)
                norm_ctrl = ctrl_intensity_matrix / np.sum(ctrl_intensity_matrix, axis=0)
                ratio_labl_raw = np.median(norm_labl, axis=1)
                ratio_ctrl_raw = np.median(norm_ctrl, axis=1)
                std_labl = np.std(norm_labl, axis=1)
                std_ctrl = np.std(norm_ctrl, axis=1)
                for num_groups in range(1, min(NUM_GROUPS + 1, max_groups)):
                    if ratio_labl_raw[num_groups * 2] > 0:
                        corrected = []
                        fill_missing = False
                        for i, uncorrected in enumerate(ratio_ctrl_raw):
                            if LABELLED is False:
                                if uncorrected == 0 or np.isnan(uncorrected):
                                    fill_missing = True
                            if fill_missing:
                                corrected.append(0)
                            else:
                                working = uncorrected
                                for j in range(i):
                                    correction_value = corrected[j] * prxu[num_groups-1][i - j]
                                    working -= correction_value
                                working /= prxu[num_groups-1][0]
                                corrected.append(working)

                        corrected = list(corrected / np.sum(corrected))
                        predicted = {}
                        for reagent in reagents_for_groups[num_groups]:
                            for C13_from_X, X_prob in enumerate(corrected):
                                total_C13 = int(reagent["num_C"] + C13_from_X)
                                total_prob = reagent["prob"] * X_prob
                                if total_C13 not in predicted:
                                    predicted[total_C13] = 0
                                predicted[total_C13] += total_prob
                        for x in predicted:
                            predicted[x] = float(predicted[x])

                        theoretical = []
                        delta = []
                        found_prob = sum([predicted[i] for i, z in enumerate(ratio_labl_raw) if z > 0])
                        for i, z in enumerate(ratio_labl_raw):
                            if predicted[i] < DETECTION_LIMIT and z == 0:
                                pass
                                #delta.append(0)
                            elif predicted[i] > DETECTION_LIMIT and z == 0:
                                theoretical.append(predicted[i])
                                var_to_use = np.mean([x for x in std_labl if x > 0]) + np.mean([x for x in std_ctrl if x > 0])
                                delta.append(((predicted[i]/found_prob - z) ** 2) / var_to_use)
                            else:
                                theoretical.append(predicted[i])
                                delta.append(((predicted[i]/found_prob - z) ** 2) / (std_labl[i] + std_ctrl[i]))
                        chi2_value = np.sum(delta)
                        p_val = 1-chi2.sf(chi2_value, df=len(delta))
                        if p_val <= THRESHOLD and len(delta) >= MIN_CLIQUE_SIZE:
                            inferred_masses = [ft_dict[x]['mz'] - (i * C13_MASS_DELTA) - (num_groups * reagent_mass_delta) for i,x in enumerate(working_isotopologues) if x],
                            results_for_root.append({
                                "clique": [x for x in working_isotopologues if x],
                                "features": [
                                    {'mz': z['mz'], 
                                    'rtime': z['rtime'],
                                    'rtime_left_base': z['rtime_left_base'],
                                    'rtime_right_base': z['rtime_right_base']} for z in [ft_dict[x] for x in working_isotopologues if x]],
                                #"reagent_mass_delta": reagent_mass_delta,
                                "num_groups": num_groups,
                                #"ratio_labl": ratio_labl_raw.tolist(),
                                "ratio_ctrl": corrected,
                                "expected_dist": predicted,
                                "std_labl": std_labl.tolist(),
                                "inferred_mass": np.mean(inferred_masses),
                                "inferred_mass_range": [np.mean(inferred_masses) - np.std(inferred_masses) * 2.5, np.mean(inferred_masses) + np.std(inferred_masses)],
                                "chi2_p": p_val,
                                "chi2_stats": delta,
                                "theory": theoretical,
                                #"raw": ratio_labl_raw.tolist()
                            })
    if results_for_root:
        if GREEDY:
            results.append(sorted(results_for_root, key=lambda x: (-len(x["clique"]), x["chi2_p"]))[0])
        else:
            results.extend(results_for_root)

            
import json
print(json.dumps(results, indent=4))