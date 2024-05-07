import pandas as pd
import sys
import matplotlib.pyplot as plt
import numpy as np
from mass2chem.formula import calculate_formula_mass

ft = pd.read_csv(sys.argv[1], sep="\t")
meta = pd.read_excel("~/Analyses/STS_sample_list.xlsx")
targets = [
    ("Estrone", "C18H22O2", 1),
    ("17β-Estradiol", "C18H24O2", 1),
    ("Estriol", "C18H24O3", 1)
]

PROTON = 1.00727647
dncl_mass = 504.22031 - 271.16926

deriv_targets = []
for t in targets:
    base_mass = calculate_formula_mass(t[1]) + PROTON
    deriv_targets.append([t[0], t[1], base_mass, 0])
    print(t[0], base_mass)

    for i in range(1, t[2] + 1):
        deriv_mass = base_mass + dncl_mass * i 
        deriv_targets.append([t[0] +  "_DnCl_" + str(i), t[1], deriv_mass, i])
print(deriv_targets)

samples = ft.columns[11:]

colors = []
samples = []
standards = []
standards_colors = []
for sample in ft.columns[11:]:
    genotype = None
    media = None
    if "Standard" in sample:
        standards.append(sample)
        if "12C13C" in sample:
            standards_colors.append('m')
        else:
            standards_colors.append('orange')
        
    else:
        for m_d in meta.to_dict(orient='records'):
            if m_d["Sample_name"] in sample:
                genotype = m_d["Genotype"]
                media = m_d["Media"]
                break
        if genotype and media:
            samples.append(sample)
            if genotype == "WT" and media == "Basal+DHEAS":
                colors.append('g')
            elif genotype == "WT" and media == "Basal":
                colors.append('b')
            elif genotype != "WT" and media == "Basal+DHEAS":
                colors.append('r')
            elif genotype != "WT" and media == "Basal":
                colors.append('k')
print("\n".join(standards))
for dt in deriv_targets:
    print(dt)
    mz_err = dt[2] * 5 / 1e6
    matches = ft[ft["mz"].between(dt[2] - mz_err, dt[2] + mz_err)]
    if matches.shape[0] == 1:
        rt = matches['rtime'].values[0]
        print(dt[2])
        C13_matches = ft[ft["mz"].between((dt[2] + 1.003355 * dt[3] * 2) - mz_err, (dt[2] + 1.003355 * dt[3] * 2) + mz_err)]
        C13_matches = C13_matches[C13_matches['rtime'].between(rt - 1, rt + 1)]
        print(C13_matches)
        if C13_matches.shape[0] == 1:
            values = matches[samples].values[0]
            s_values = matches[standards].values[0]
            C13_values = C13_matches[samples].values[0]
            C13_s_values = C13_matches[standards].values[0]

            if sum(values):
                if dt[3] > 0:
                    fig, (ax1, ax2, ax3) = plt.subplots(1,3)
                    comb_colors = [*colors, *standards_colors]
                    all_samples = [*samples, *standards]
                    all_values = np.log10([x+1 for x in [*values, *s_values]])
                    plt.suptitle(dt[0] + " rtime: " + str(rt) + "(s)")
                    ax1.set_ylim([0, 10])
                    ax1.bar(list(range(len(all_samples))), all_values, color=comb_colors)
                    C13_all_values = np.log10([x+1 for x in [*C13_values, *C13_s_values]])
                    ax2.set_ylim([0, 10])
                    ax2.bar(list(range(len(all_samples))), C13_all_values, color=comb_colors)
                    ratio = []
                    for a, b in zip([*values, *s_values], [*C13_values, *C13_s_values]):
                        ratio.append(a/(a+b) if b > 0 else 0)
                    ax3.set_ylim([0,1])
                    ax3.bar([x[0] for x in enumerate(ratio)], ratio)
                    ax3.set_ylabel("12C/13C Ratio")
                    plt.show()
                else:
                    plt.suptitle(dt[0] + " rtime: " + str(rt) + "(s)")
                    comb_colors = [*colors, *standards_colors]
                    all_samples = [*samples, *standards]
                    all_values = np.log10([x+1 for x in [*values, *s_values]])
                    plt.bar(list(range(len(all_samples))), all_values, color=comb_colors)
                plt.show()