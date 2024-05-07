import pandas as pd
import sys
import matplotlib.pyplot as plt
import numpy as np
from mass2chem.formula import calculate_formula_mass

ft = pd.read_csv(sys.argv[1], sep="\t")
meta = pd.read_excel("~/Analyses/STS_sample_list.xlsx")
targets = [
    ("Desoxycorticosterone", "C21H30O3", 2),
    ("17alpha-hydroxyprogesterone", "C21H30O3", 2),
    ("Aldosterone", "C21H28O5", 3),
    ("Androstenedione", "C19H26O2", 2),
    ("Androsterone", "C19H30O2", 1),
    ("Corticosterone", "C21H30O4", 2),
    ("Cortisone", "C21H28O5", 3),
    ("Prasterone", "C19H28O2", 1),
    ("DHEA", "C19H28O2", 1),
    ("Stanolone", "C19H30O2", 1),
    ("Eticholanone", "C19H30O2", 1),
    ("Hydrocortisone", "C21H30O5", 2),
    ("Progesterone", "C21H30O2", 2),
    ("Testosterone", "C19H28O2", 1),
    ("11-Deoxycortisol", "C21H30O4", 2),
    ("DHEA-S", "C19H28O5S", 0)

]

PROTON = 1.00727647
dnhz_mass = 578.3047-331.22677

deriv_targets = []
for t in targets:
    base_mass = calculate_formula_mass(t[1]) + PROTON
    deriv_targets.append([t[0], t[1], base_mass, 0])
    for i in range(1, t[2] + 1):
        deriv_mass = base_mass + dnhz_mass * i
        deriv_targets.append([t[0] +  "_DnHz_" + str(i), t[1], deriv_mass, i])
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

for dt in deriv_targets:
    mz_err = dt[2] * 5 / 1e6
    matches = ft[ft["mz"].between(dt[2] - mz_err, dt[2] + mz_err)]
    for match in matches.to_dict(orient='records'):
        rt = match['rtime']
        C13_matches = ft[ft["mz"].between((dt[2] + 1.003355 * dt[3] * 2) - mz_err, (dt[2] + 1.003355 * dt[3] * 2) + mz_err)]
        C13_matches = C13_matches[C13_matches['rtime'].between(rt - 1, rt + 1)]
        if not C13_matches.empty:
            for C13_match in C13_matches.to_dict(orient='records'):
                values = [match[x] for x in samples]
                s_values = [match[x] for x in standards]
                C13_values = [C13_match[x] for x in samples]
                C13_s_values = [C13_match[x] for x in standards]
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
        else:
            values = [match[x] for x in samples]
            if sum(values):
                plt.suptitle(dt[0] + " rtime: " + str(rt) + "(s)")
                plt.bar(list(range(len(samples))), np.log10([x+1 for x in values]), color=colors)
                plt.show()            