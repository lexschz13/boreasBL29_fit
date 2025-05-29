import boreas_fit as bf


bf.data_to_csv(["Tue27", "Wed28", "Thu29"], "boreasBL29-gherranz-May2025")


T, gap, errgap, rel_errgap, cost, rel_cost = bf.compute_gap("boreasBL29-gherranz-May2025", 'LH', 'S1', 'Off', 0,
                        mode="lorentz", emin=528.5, emax=532, ylim_gaps=(0.7,1.0))