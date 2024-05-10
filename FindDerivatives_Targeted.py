import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mass2chem.formula import calculate_formula_mass

ft = pd.read_csv(sys.argv[1], sep="\t")
targets = pd.read_csv(sys.argv[2])
reagents = pd.read_csv(sys.argv[3])
sample_metadata = pd.read_excel(sys.argv[4])

PROTON = 1.00727647

def generate_derivatized_targets(targets, reagents):
    derivatized_targets = []
    reagent_mass_vector = reagents['MassDelta'].values.tolist()
    reagent_name_vector = reagents['Name'].values.tolist()
    for t in targets.to_dict(orient='records'):
        t['Underiv_Mass'] = calculate_formula_mass(t['Formula']) + PROTON
        max_reagent_vector = []
        for rname in reagent_name_vector:
            max_reagent_vector.append(t.get(rname, 0) + 1)
        for reagent_combination in np.ndindex(tuple(max_reagent_vector)):
            dt = {k:v for k,v in t.items()}
            mass_delta = np.dot(reagent_combination, reagent_mass_vector)
            dt['Deriv_Mass'] = dt['Underiv_Mass'] + mass_delta
            dt['Num_Rxns'] = sum(reagent_combination)
            dt['Deriv_Name'] = dt['TargetName'] + "_" + ','.join([str(x) + y for x,y in zip(reagent_combination, reagent_name_vector)])
            for rname, rcount in zip(reagent_name_vector, reagent_combination):
                dt[rname] = rcount
            derivatized_targets.append(dt)
    return derivatized_targets

def color_groups(sample_metadata):
    group_colors = {}
    i = 0
    for sample in sample_metadata.to_dict(orient='records'):
        if sample["GROUP"] not in group_colors:
            c = list(mcolors.BASE_COLORS.keys())[i]
            if c == "w":
                c = "grey"
            group_colors[sample["GROUP"]] = c
            i += 1
    group_colors['LBL_STANDARDS'] = "orange"
    group_colors['STANDARDS'] = 'cyan'
    return group_colors

def map_samples_to_groups(sample_metadata, table):
    sname_to_group = {}
    for sample in table.columns[11:]:
        for meta in sample_metadata.to_dict(orient='records'):
            if meta["Sample_name"] in sample:
                sname_to_group[sample] = meta["GROUP"]
    for sample in table.columns[11:]:
        if sample not in sname_to_group:
            if "Standard" in sample:
                if "12C13C" in sample:
                    sname_to_group[sample] = "LBL_STANDARDS"
                else:
                    sname_to_group[sample] = "STANDARDS"
    return sname_to_group

def map_samples_to_colors(sample_group_map, group_color_map):
    sample_color_map = {}
    for k in sample_group_map:
        sample_color_map[k] = group_color_map[sample_group_map[k]]
    return sample_color_map

def search_for_targets(derivatized_targets, table, mz_tol=5):
    results = []
    for dt in derivatized_targets:
        mz_err = dt['Deriv_Mass'] * mz_tol / 1e6
        mz_matches = table[table['mz'].between(dt['Deriv_Mass'] - mz_err, dt['Deriv_Mass'] + mz_err)]
        for mz_match in mz_matches.to_dict(orient='records'):
            rtime = mz_match['rtime']
            mass_delta_13C = dt['Num_Rxns'] * 2 * 1.003355
            iso_matches_r = table[table['rtime'].between(rtime - 5, rtime + 5)]
            iso_matches_r_mz = iso_matches_r[iso_matches_r['mz'].between(dt['Deriv_Mass'] + mass_delta_13C - mz_err, dt['Deriv_Mass'] + mass_delta_13C + mz_err)]
            if not iso_matches_r_mz.empty:
                for C13_match in iso_matches_r_mz.to_dict(orient='records'):
                    results.append({
                        "Target": dt["Deriv_Name"],
                        "m0": mz_match["id_number"],
                        "m13C": C13_match["id_number"]
                    })
            else:
                results.append({
                    "Target": dt["Deriv_Name"],
                    "m0": mz_match["id_number"],
                    "m13C": None
                })
    return results

def plot_pairs(result_pairs, table, s_g_map, s_c_map):
    samples = []
    samples_colors = []
    standards = []
    standards_colors = []
    for (s, g), (_, c) in zip(s_g_map.items(), s_c_map.items()):
        if g != "LBL_STANDARDS" and g != "STANDARDS":
            samples.append(s)
            samples_colors.append(c)
        else:
            standards.append(s)
            standards_colors.append(c)


    for pair in result_pairs:
        m0_feature = pair["m0"] 
        m0_samples_values = table[table['id_number'] == m0_feature][samples].values[0].tolist()
        m0_standards_values = table[table['id_number'] == m0_feature][standards].values[0].tolist()
        if sum(m0_samples_values):
            m13C_feature = pair["m13C"]
            if m13C_feature is not None:
                m13C_samples_values = table[table['id_number'] == m13C_feature][samples].values[0].tolist()
                m13C_standards_values = table[table['id_number'] == m13C_feature][standards].values[0].tolist()

            if m13C_feature is not None:
                plt.bar(list(range(len(m0_samples_values + m0_standards_values))), m0_samples_values + m0_standards_values, color=samples_colors + standards_colors)
                plt.suptitle(pair["Target"])
                plt.show()
            else:
                _, (ax1, ax2, ax3) = plt.subplot(1,3)
                ax1.bar(list(range(len(m0_samples_values + m0_standards_values))), np.log10(m0_samples_values + m0_standards_values), color=samples_colors + standards_colors)
                ax2.bar(list(range(len(m13C_samples_values + m13C_standards_values))), np.log10(m13C_samples_values + m13C_standards_values), color=samples_colors + standards_colors)
                ratio = []
                for a, b in zip(m0_samples_values + m0_standards_values, m13C_samples_values + m13C_standards_values):
                    ratio.append(a/(a+b))
                ax3.bar(list(range(len(ratio))), ratio)
                plt.suptitle(pair["Target"])
                plt.show()


dtargets = generate_derivatized_targets(targets, reagents)
g_c_map = color_groups(sample_metadata)
s_g_map = map_samples_to_groups(sample_metadata, ft)
s_c_map = map_samples_to_colors(s_g_map, g_c_map)
deriv_pairs = search_for_targets(dtargets, ft)
plot_pairs(deriv_pairs, ft, s_g_map, s_c_map)